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

const struct mat33 CRAZYFLIE_INERTIA = {{{16.6e-6f, 0.83e-6f, 0.72e-6f},
                                         {0.83e-6f, 16.6e-6f, 1.8e-6f},
                                         {0.72e-6f, 1.8e-6f, 29.3e-6f}}};

// static const float THRUST_MIN = 1;
// static const float THRUST_MAX = 40;
static const float ROTATION_MAX = 30;

// const gains
static float TRANS_KP_FIXED[] = {16.0f, 16.0f, 12.0f};
static float TRANS_KD_FIXED[] = {12.0f, 12.0f, 5.0f};
static float ROT_KP_FIXED[] = {95.0f, 95.0f, 95.0f};
static float ROT_KD_FIXED[] = {40.0f, 40.0f, 40.0f};

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

void controllerOutOfTreeInit() { isInit = true; }

bool controllerOutOfTreeTest() {
  // Always return true
  return true;
}

void controllerOutOfTree(control_t *control, const setpoint_t *setpoint,
                         const sensorData_t *sensors, const state_t *state,
                         const stabilizerStep_t stabilizerStep) {
  float omega[3] = {0};
  //  float fth = 0;
  float control_thrust = 0;
  struct vec z_vec = mkvec(0, 0, 1);
  //  struct vec unitz = mkvec(0, 0, 1);
  struct quat z_q = mkquat(0, 0, 0, 1);
  omega[0] = radians(sensors->gyro.x);
  omega[1] = radians(sensors->gyro.y);
  omega[2] = radians(sensors->gyro.z);

  if (RATE_DO_EXECUTE(ATTITUDE_RATE, stabilizerStep) ||
      RATE_DO_EXECUTE(POSITION_RATE, stabilizerStep)) {

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

    omega_stored[0] = omega[0];
    omega_stored[1] = omega[1];
    omega_stored[2] = omega[2];

    // ----- Translational control ------

    struct vec trans_kp_vec =
        mkvec(TRANS_KP_FIXED[0], TRANS_KP_FIXED[1], TRANS_KP_FIXED[2]);
    struct vec trans_kd_vec =
        mkvec(TRANS_KD_FIXED[0], TRANS_KD_FIXED[1], TRANS_KD_FIXED[2]);

    struct vec trans_control =
        vadd(veltmul(trans_kp_vec, posError), veltmul(trans_kd_vec, velError));
    trans_control.z -= GRAVITY_MAGNITUDE;
    trans_control = vscl(-CF_MASS, trans_control);
    float norm_trans_control = vmag(trans_control);
    struct vec control_direction = vzero();

    if (norm_trans_control != 0) {
      // control_direction = vdiv(trans_control,norm_trans_control);
      // trans_control = vscl(THRUST_MAX*tanhf(norm_trans_control/THRUST_MAX),
      // control_direction);
      control_direction = vdiv(trans_control, norm_trans_control);
    }

    struct quat curr_thrust_force_vectorq = qqmul(orientation, z_q);
    curr_thrust_force_vectorq =
        qqmul(curr_thrust_force_vectorq, qinv(orientation));
    struct vec curr_thrust_force_vector =
        mkvec(curr_thrust_force_vectorq.x, curr_thrust_force_vectorq.y,
              curr_thrust_force_vectorq.z); // Fth
    control_thrust =
        trans_control.z / vdot(z_vec, curr_thrust_force_vector); // Fu

    /*
    // Current attitude
    struct quat orientation =
        mkquat(state->attitudeQuaternion.x, state->attitudeQuaternion.y,
               state->attitudeQuaternion.z, state->attitudeQuaternion.w);
    if (!recorded_initial_orientation) {
      initial_orientation = orientation;
      recorded_initial_orientation = true;
    }

    struct vec posError = mkvec(state->position.x - setpoint->position.x,
                                state->position.y - setpoint->position.y,
                                state->position.z - setpoint->position.z);
    //    struct vec pos =
    //        mkvec(state->position.x, state->position.y, state->position.z);
    // pos = qvrot(initial_orientation, pos);
    // struct vec posError = mkvec(pos.x - 4, pos.y - 0, pos.z - 2);

    stored_target[0] = setpoint->position.x;
    stored_target[1] = setpoint->position.y;
    stored_target[2] = setpoint->position.z;

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

    omega_stored[0] = omega[0];
    omega_stored[1] = omega[1];
    omega_stored[2] = omega[2];

    kp_x = TRANS_KP_FIXED[0];
    kp_y = TRANS_KP_FIXED[1];
    kp_z = TRANS_KP_FIXED[2];

    struct vec trans_kp_vec =
        mkvec(TRANS_KP_FIXED[0], TRANS_KP_FIXED[1], TRANS_KP_FIXED[2]);
    struct vec trans_kd_vec =
        mkvec(TRANS_KD_FIXED[0], TRANS_KD_FIXED[1], TRANS_KD_FIXED[2]);
    struct vec ut = vadd(veltmul(vneg(trans_kp_vec), posError),
                         veltmul(vneg(trans_kd_vec), velError));
    ut.z += GRAVITY_MAGNITUDE * CF_MASS;

    fth = vmag(ut);
    struct vec fu = vzero();

    if (fth != 0) {
      fu = vdiv(ut, fth);
    }

    stored_fth[0] = ut.x;
    stored_fth[1] = ut.y;
    stored_fth[2] = ut.z;

    struct quat curr_thrust_force_vectorq = qqmul(orientation, z_q);
    curr_thrust_force_vectorq =
        qqmul(curr_thrust_force_vectorq, qinv(orientation));
    struct vec curr_thrust_force_vector =
        mkvec(curr_thrust_force_vectorq.x, curr_thrust_force_vectorq.y,
              curr_thrust_force_vectorq.z);                        // Fth
    control_thrust = ut.z / vdot(unitz, curr_thrust_force_vector); // Fu

    float dot_temp = vdot(fu, unitz);
    struct vec vcross_temp = vcross(fu, unitz);

    struct vec imaginary = mkvec(0, 0, 0);
    float cross_mag = vmag(vcross_temp);
    if (cross_mag != 0) {
      imaginary = vscl(sqrtf((1 - dot_temp) / 2), vdiv(vcross_temp, cross_mag));
    }

    float real = sqrtf((1 + dot_temp) / 2);

    struct quat qd = mkquat(-imaginary.x, -imaginary.y, -imaginary.z, real);
    struct quat qz = mkquat(0, 0, 0, 1);
    struct quat qe = qnormalize(qqmul(qqmul(qinv(qz), qd), orientation));
    float theta = 2 * acosf(qe.w);
    struct vec qrv = vscl(theta, mkvec(qe.x, qe.y, qe.z));
    if (vmag(qrv) > M_PI_F) {
      qd = qneg(qd);
      qe = qnormalize(qqmul(qqmul(qinv(qz), qd), orientation));
    }

    struct vec orientationErrorVector = qlog(qd);

    struct vec angVelocityErrorVector =
        mkvec(omega[0] - setpoint->attitudeRate.pitch,
              omega[1] - setpoint->attitudeRate.roll,
              omega[2] - setpoint->attitudeRate.yaw);

    struct vec angVelocityVector = mkvec(omega[0], omega[1], omega[2]);

    store_from_q(qd, orientation_stored);
    store_from_q(qe, orientation_error_stored);

    struct vec rot_kp_vec =
        mkvec(ROT_KP_FIXED[0], ROT_KP_FIXED[1], ROT_KP_FIXED[2]);
    struct vec rot_kd_vec =
        mkvec(ROT_KD_FIXED[0], ROT_KD_FIXED[1], ROT_KD_FIXED[2]);

    angular_velocity_error_stored[0] = angVelocityErrorVector.x;
    angular_velocity_error_stored[1] = angVelocityErrorVector.y;
    angular_velocity_error_stored[2] = angVelocityErrorVector.z;

    struct vec orientation_control =
        vadd(veltmul(vneg(rot_kp_vec), orientationErrorVector),
             veltmul(vneg(rot_kd_vec), angVelocityErrorVector));
    // struct vec control_torque =
    control_torque = vadd(
        orientation_control,
        vcross(angVelocityVector, mvmul(CRAZYFLIE_INERTIA, angVelocityVector)));
    */
  }

  // Control Input
  if (setpoint->mode.z == modeDisable) {
    control->thrustSi = 0.0f;
    control->torque[0] = 0.0f;
    control->torque[1] = 0.0f;
    control->torque[2] = 0.0f;
  } else {
    // control the body torques
    control->thrustSi = control_thrust;
    control->torqueX = control_torque.x;
    control->torqueY = control_torque.y;
    control->torqueZ = control_torque.z;
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

LOG_ADD(LOG_FLOAT, fth_x, &stored_fth[0])
LOG_ADD(LOG_FLOAT, fth_y, &stored_fth[1])
LOG_ADD(LOG_FLOAT, fth_z, &stored_fth[2])

// LOG_ADD(LOG_FLOAT, target_x, &stored_target[0])
// LOG_ADD(LOG_FLOAT, target_y, &stored_target[1])
// LOG_ADD(LOG_FLOAT, target_z, &stored_target[2])

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
