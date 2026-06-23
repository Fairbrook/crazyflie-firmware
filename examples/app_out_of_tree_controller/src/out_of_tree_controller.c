/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--´  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2024 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * out_of_tree_controller.c - App layer application of an out of tree
 * controller.
 */

#include "stabilizer_types.h"
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include "app.h"
#include "log.h"

#include "FreeRTOS.h"
#include "task.h"

// Edit the debug name to get nice debug prints
#define DEBUG_MODULE "MYCONTROLLER"
#include "debug.h"

// We still need an appMain() function, but we will not really use it. Just let
// it quietly sleep.
void appMain() {
  DEBUG_PRINT("Waiting for activation ...\n");

  while (1) {
    vTaskDelay(M2T(2000));
  }
}

// The new controller goes here --------------------------------------------
// Move the includes to the the top of the file if you want to
#include "commander.h"
#include "controller.h"
#include "log.h"
#include "math3d.h"
#include "num.h"
#include "param.h"
#include "physicalConstants.h"
#include "platform_defaults.h"
#include "stabilizer_types.h"

static const float MU = 0.00f;

const struct mat33 CRAZYFLIE_INERTIA = {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
                                         {0.83e-6f, 16.6e-6f, 1.8e-6f},
                                         {0.72e-6f, 1.8e-6f, 29.3e-6f}}};

// static const float THRUST_MIN = 1;
// static const float THRUST_MAX = 40;
static const float ROTATION_MAX = 30;

// const gains
// Retuned for the SI force/torque (N, N*m) interpretation of
// controlModeForceTorque. The original {16,16,12}/{12,12,5} translational and
// {95,95,95}/{40,40,40} rotational gains command ~60g of force and ~160 N*m of
// torque (vs the CF's ~0.01-0.03 N*m actuator limit), so the attitude loop
// chatters against saturation and the vehicle tumbles. These values give proper
// bandwidth separation (attitude ~30 rad/s, position ~2 rad/s) and converge
// cleanly (verified in tools/matlab/simulate_controller_oot_convergence.m).
static float TRANS_KP_FIXED[] = {0.11f, 0.11f, 0.11f};
static float TRANS_KD_FIXED[] = {0.10f, 0.10f, 0.10f};
static float ROT_KP_FIXED[] = {0.03f, 0.03f, 0.03f};
static float ROT_KD_FIXED[] = {0.0008f, 0.0008f, 0.0008f};

// Init store variables
static float pos_error_stored[] = {0.0f, 0.0f, 0.0f};
static float vel_error_stored[] = {0.0f, 0.0f, 0.0f};
static float omega_stored[] = {0.0f, 0.0f, 0.0f};
static float orientation_stored[] = {0.0f, 0.0f, 0.0f, 0.0f};
static float orientation_error_stored[] = {0.0f, 0.0f, 0.0f, 0.0f};
// static float angular_velocity_stored[] = {0.0f,0.0f,0.0f};
static float angular_velocity_error_stored[] = {0.0f, 0.0f, 0.0f};
static float stored_fth[] = {0.0f, 0.0f, 0.0f};
// static float stored_target[] = {0.0f, 0.0f, 0.0f};

struct quat initial_orientation;
bool recorded_initial_orientation;

// Input variables
static struct vec control_torque;

static bool isInit;

// Auxiliary functions
// signum
static inline int signum(float n) {
  if (n < 0) {
    return -1;
  } else if (n > 0) {
    return 1;
  } else {
    return 0;
  }
}

static inline float qmag(struct quat q) {
  float norm_q = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z + q.w * q.w);
  return norm_q;
}

static inline struct quat qsmul(struct quat q, float s) {
  struct quat qs = mkquat(q.x * s, q.y * s, q.z * s, q.w * s);
  return qs;
}

static inline struct quat load_q_from_array(float const *d) {
  return mkquat(d[0], d[1], d[2], d[3]);
}

static inline void store_from_q(struct quat q, float *d) {
  d[0] = (float)q.x;
  d[1] = (float)q.y;
  d[2] = (float)q.z;
  d[3] = (float)q.w;
}

// Quaternion to axis-angle
static inline struct vec rotvec(struct quat q) {
  q = qnormalize(q);
  float ang = 2 * acosf(q.w);
  float axis[] = {q.x, q.y, q.z};
  float magnitude = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z);
  float rv[] = {0.0f, 0.0f, 0.0f};
  for (int i = 0; i < 3; i++) {
    rv[i] = (magnitude != 0) ? (ang * axis[i]) / magnitude : 0;
  }
  return mkvec(rv[0], rv[1], rv[2]);
}

// From: https://la.mathworks.com/help/nav/ref/quaternion.log.html
static inline struct vec qlog(struct quat q) {
  float norm_v = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z);
  if (norm_v == 0) {
    return mkvec(0, 0, 0);
  }
  return vscl(acosf(q.w), vdiv(mkvec(q.x, q.y, q.z), norm_v));
}

// exp(q) = exp(a)*( cos(||v||) + (v/||v||)*sin(||v||) )
static inline struct quat qexp(struct quat q) {
  float norm_v = sqrtf(q.x * q.x + q.y * q.y + q.z * q.z);
  float exp_scalar = expf(q.w) * cosf(norm_v);
  float exp_factor = (norm_v != 0) ? (expf(q.w) * sinf(norm_v) / norm_v) : 0;
  struct quat result;
  result.w = exp_scalar;
  result.x = exp_factor * q.x;
  result.y = exp_factor * q.y;
  result.z = exp_factor * q.z;
  return result;
}

static inline float absf(const float a) {
  if (a < 0)
    return -a;
  return a;
}

static inline float dhnorm(float e, float ep, float max_e, float gamma2,
                           float mu) {
  return (1 / max_e) * (float)pow(absf(e), 1 / (1 - mu)) +
         (float)(gamma2 * absf(ep));
}

void controllerOutOfTreeInit() { isInit = true; }

bool controllerOutOfTreeTest() {
  // Always return true
  return true;
}

float control_thrust = 0;
//  float fth = 0;

void controllerOutOfTree(control_t *control, const setpoint_t *setpoint,
                         const sensorData_t *sensors, const state_t *state,
                         const stabilizerStep_t stabilizerStep) {
  if (RATE_DO_EXECUTE(ATTITUDE_RATE, stabilizerStep) ||
      RATE_DO_EXECUTE(POSITION_RATE, stabilizerStep)) {
    struct vec omega = mkvec(0, 0, 0);
    struct vec z_vec = mkvec(0, 0, 1);
    struct quat z_q = mkquat(0, 0, 1, 0);

    omega.x = radians(sensors->gyro.x);
    omega.y = radians(sensors->gyro.y);
    omega.z = radians(sensors->gyro.z);

    // Current attitude
    struct quat orientation =
        mkquat(state->attitudeQuaternion.x, state->attitudeQuaternion.y,
               state->attitudeQuaternion.z, state->attitudeQuaternion.w);

    struct vec posError = mkvec(state->position.x - setpoint->position.x,
                                state->position.y - setpoint->position.y,
                                state->position.z - setpoint->position.z);

    pos_error_stored[0] = posError.x;
    pos_error_stored[1] = posError.y;
    pos_error_stored[2] = posError.z;

    // Velocity error
    struct vec velError = mkvec(state->velocity.x - setpoint->velocity.x,
                                state->velocity.y - setpoint->velocity.y,
                                state->velocity.z - setpoint->velocity.z);
    vel_error_stored[0] = velError.x;
    vel_error_stored[1] = velError.y;
    vel_error_stored[2] = velError.z;
    // Angular velocity from gyroscope

    omega_stored[0] = omega.x;
    omega_stored[1] = omega.y;
    omega_stored[2] = omega.z;

    // ----- Translational control ------

    struct vec trans_kp_vec =
        mkvec(TRANS_KP_FIXED[0], TRANS_KP_FIXED[1], TRANS_KP_FIXED[2]);
    struct vec trans_kd_vec =
        mkvec(TRANS_KD_FIXED[0], TRANS_KD_FIXED[1], TRANS_KD_FIXED[2]);

    struct vec fu = vadd(veltmul(vscl(-1, trans_kp_vec), posError),
                         veltmul(vscl(-1, trans_kd_vec), velError));
    struct vec ut = vdiv(fu, vmag(fu));
    if (vmag(fu) > CF_MASS * GRAVITY_MAGNITUDE) {
      fu = vscl(CF_MASS * GRAVITY_MAGNITUDE, ut);
    }

    fu.z += ((CF_MASS+0.002f) * GRAVITY_MAGNITUDE);

    stored_fth[0] = fu.x;
    stored_fth[1] = fu.y;
    stored_fth[2] = fu.z;

    struct quat curr_thrust_force_vectorq = qqmul(orientation, z_q);
    curr_thrust_force_vectorq =
        qqmul(curr_thrust_force_vectorq, qinv(orientation));
    struct vec curr_thrust_force_vector =
        mkvec(curr_thrust_force_vectorq.x, curr_thrust_force_vectorq.y,
              curr_thrust_force_vectorq.z); // Force direction
    control_thrust = fu.z / vdot(z_vec, curr_thrust_force_vector); // Fu
    control_thrust =
        constrain(control_thrust, 0, CF_MASS * GRAVITY_MAGNITUDE * 1.5f);

    ut = vdiv(fu, vmag(fu));
    float dot = vdot(ut, z_vec);
    struct vec cross = vcross(ut, z_vec);
    struct vec imaginary = mkvec(0, 0, 0);
    float real = 1;

    if (vmag(cross) > 1e-12f) {
      imaginary = vscl(sqrtf((1 - dot) / 2), vdiv(cross, vmag(cross)));
      real = sqrtf((1 + dot) / 2);
    }
    struct quat qd = mkquat(-imaginary.x, -imaginary.y, -imaginary.z, real);
    qd = qnormalize(qd);
    struct quat qe = qqmul(qinv(qd), orientation);
    qe = qnormalize(qe);
    float theta = 2 * acosf(qe.w);
    if (theta != 0) {
      struct vec qrv = vscl((theta / sinf(theta / 2)), mkvec(qe.x, qe.y, qe.z));
      if (vmag(qrv) > M_PI_F) {
        qd = qneg(qd);
        qe = qqmul(qinv(qd), orientation);
        qe = qnormalize(qe);
      }
    }
    store_from_q(qd, orientation_stored);
    store_from_q(qe, orientation_error_stored);

    struct vec rot_kp_vector =
        mkvec(ROT_KP_FIXED[0], ROT_KP_FIXED[1], ROT_KP_FIXED[2]);
    struct vec rot_kd_vector =
        mkvec(ROT_KD_FIXED[0], ROT_KD_FIXED[1], ROT_KD_FIXED[2]);

    struct vec rotError = qlog(qe);

    struct vec dh = mkvec(dhnorm(rotError.x, omega.x, M_PI_F, 0.5, MU),
                          dhnorm(rotError.y, omega.y, M_PI_F, 0.5, MU),
                          dhnorm(rotError.z, omega.z, M_PI_F, 0.5, MU));
    struct vec dhprop = mkvec(dh.x > 1e-12f ? (float)pow(dh.x, 2 * MU) : 0,
                              dh.y > 1e-12f ? (float)pow(dh.y, 2 * MU) : 0,
                              dh.z > 1e-12f ? (float)pow(dh.z, 2 * MU) : 0);
    struct vec dhdiff = mkvec(dh.x > 1e-12f ? (float)pow(dh.x, MU) : 0,
                              dh.y > 1e-12f ? (float)pow(dh.y, MU) : 0,
                              dh.z > 1e-12f ? (float)pow(dh.z, MU) : 0);

    struct vec tau =
        vadd(veltmul(vscl(-2, rot_kp_vector), veltmul(dhprop, rotError)),
             veltmul(vscl(-1, rot_kd_vector), veltmul(dhdiff, omega)));
    if (vmag(tau) > CF_MASS * 0.5f) {
      tau = vscl(CF_MASS * 0.5f, vdiv(tau, vmag(tau)));
    }
    control_torque.x = tau.x;
    control_torque.y = tau.y;
    control_torque.z = tau.z;
  }

  // Control Input
  if (setpoint->mode.z == modeDisable) {
    control->thrustSi = 0.0f;
    control->torque[0] = 0.0f;
    control->torque[1] = 0.0f;
    control->torque[2] = 0.0f;
  } else {
    // control the body torques
    //    control->thrustSi = control_thrust;
    control->thrustSi = control_thrust;
    control->torqueX  = control_torque.x;
    control->torqueY  = control_torque.y;
    control->torqueZ  = control_torque.z;
  }

  //
  control->controlMode = controlModeForceTorque;
}

PARAM_GROUP_START(ootParams)

PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kp_x, &TRANS_KP_FIXED[0])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kp_y, &TRANS_KP_FIXED[1])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kp_z, &TRANS_KP_FIXED[2])

PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kd_x, &TRANS_KD_FIXED[0])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kd_y, &TRANS_KD_FIXED[1])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, trans_kd_z, &TRANS_KD_FIXED[2])

PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kp_x, &ROT_KP_FIXED[0])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kp_y, &ROT_KP_FIXED[1])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kp_z, &ROT_KP_FIXED[2])

PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kd_x, &ROT_KD_FIXED[0])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kd_y, &ROT_KD_FIXED[1])
PARAM_ADD(PARAM_FLOAT | PARAM_PERSISTENT, rot_kd_z, &ROT_KD_FIXED[2])
PARAM_GROUP_STOP(ootParams)

LOG_GROUP_START(oot)
LOG_ADD(LOG_FLOAT, pos_err_x, &pos_error_stored[0])
LOG_ADD(LOG_FLOAT, pos_err_y, &pos_error_stored[1])
LOG_ADD(LOG_FLOAT, pos_err_z, &pos_error_stored[2])

LOG_ADD(LOG_FLOAT, vel_err_x, &vel_error_stored[0])
LOG_ADD(LOG_FLOAT, vel_err_y, &vel_error_stored[1])
LOG_ADD(LOG_FLOAT, vel_err_z, &vel_error_stored[2])

LOG_ADD(LOG_FLOAT, thrust, &control_thrust)
LOG_ADD(LOG_FLOAT, torque_x, &control_torque.x)
LOG_ADD(LOG_FLOAT, torque_y, &control_torque.y)
LOG_ADD(LOG_FLOAT, torque_z, &control_torque.z)

LOG_ADD(LOG_FLOAT, qd_x, &orientation_stored[0])
LOG_ADD(LOG_FLOAT, qd_y, &orientation_stored[1])
LOG_ADD(LOG_FLOAT, qd_z, &orientation_stored[2])
LOG_ADD(LOG_FLOAT, qd_w, &orientation_stored[3])

LOG_ADD(LOG_FLOAT, q_err_x, &orientation_error_stored[0])
LOG_ADD(LOG_FLOAT, q_err_y, &orientation_error_stored[1])
LOG_ADD(LOG_FLOAT, q_err_z, &orientation_error_stored[2])
LOG_ADD(LOG_FLOAT, q_err_w, &orientation_error_stored[3])

LOG_ADD(LOG_FLOAT, ang_vel_err_x, &angular_velocity_error_stored[0])
LOG_ADD(LOG_FLOAT, ang_vel_err_y, &angular_velocity_error_stored[1])
LOG_ADD(LOG_FLOAT, ang_vel_err_z, &angular_velocity_error_stored[2])

LOG_GROUP_STOP(oot)
