// Assume that the IMU observation is 3-dimensional acceleration
_float_t z_imu[3] = {accel_x, accel_y, accel_z};  //  Observation value
_float_t hx_imu[3];  // Predicted observed value (needs to be calculated based on the model)
_float_t H_imu[3*OEKF_N];  //  3x14 observation matrix
_float_t R_imu[3*3] = {0.1,0,0, 0,0.1,0, 0,0,0.1};  //  diagonal noise matrix

//  Call lightweight update (m=3)
oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);

/ octonion_rotate(&ekf->state.q, accel_body, hx_imu);

// Observation prediction after the coupling term:
_float_t hx_imu[3];
octonion_rotate(&ekf->state.q, accel_body, hx_imu);
// Add the correction of i[5] (sequential deviation) and i[6] (disturbance coupling) to the observation
hx_imu[0] += ekf->state.q.i[5] * 0.1;  // Correct the X-axis for sequential bias
hx_imu[1] += ekf->state.q.i[5] * 0.1;  // Correct the Y-axis for sequential bias
hx_imu[2] += ekf->state.q.i[6] * 0.05;  // Perturbation coupling correction Z-axis

// Construct the H matrix (supplement the partial derivatives of the coupling terms):
_float_t H_imu[3*OEKF_N] = {0};
// Original attitude and velocity related partial derivatives (omitted)...
// Add the partial derivatives of i[5] and i[6] (corresponding to x[6] and x[7])
H_imu[0*OEKF_N + 6] = 0.1;  // Sensitivity of X-axis observation to i[5]
H_imu[1*OEKF_N + 6] = 0.1;  // Sensitivity of Y-axis observation to i[5]
H_imu[2*OEKF_N + 7] = 0.05;  // Sensitivity of Z-axis observation to i[6]

// GPS observation model (6 dimensions: 3 positions + 3 velocities)
_float_t z_gps[6] = {gps_p_x, gps_p_y, gps_p_z, gps_v_x, gps_v_y, gps_v_z}; // Observed values (position x/y/z, velocity x/y/z)
_float_t hx_gps[6]; // Predicted observed value
_float_t H_gps[6*OEKF_N] = {0}; // 6x14 observation matrix
_float_t R_gps[6*6] = { // Observation noise matrix (diagonal matrix, set according to sensor accuracy)
    1.0, 0, 0, 0, 0, 0,
    0, 1.0, 0, 0, 0, 0,
    0, 0, 1.5, 0, 0, 0,  // High-level noise can be slightly larger
    0, 0, 0, 0.1, 0, 0,
    0, 0, 0, 0, 0.1, 0,
    0, 0, 0, 0, 0, 0.15
};

// Predicted observed value hx_gps: directly extract position and velocity from the state vector
hx_gps[0] = ekf->state.p[0]; // Position x (corresponding to x[11])
hx_gps[1] = ekf->state.p[1]; // Position y (corresponding to x[12])
hx_gps[2] = ekf->state.p[2]; // Position z (corresponding to x[13])
hx_gps[3] = ekf->state.v[0]; // Velocity x (corresponding to x[8])
hx_gps[4] = ekf->state.v[1]; // Velocity y (corresponding to x[9])
hx_gps[5] = ekf->state.v[2]; // Velocity z (corresponding to x[10])

// Construct the H matrix (the partial derivative of the observation with respect to the state)
// The partial derivative of position x with respect to state x[11] is 1
H_gps[0*OEKF_N + 11] = 1.0;
// The partial derivative of position y with respect to state x[12] is 1
H_gps[1*OEKF_N + 12] = 1.0;
// The partial derivative of position z with respect to state x[13] is 1
H_gps[2*OEKF_N + 13] = 1.0;
// The partial derivative of velocity x with respect to state x[8] is 1
H_gps[3*OEKF_N + 8] = 1.0;
// The partial derivative of velocity y with respect to state x[9] is 1
H_gps[4*OEKF_N + 9] = 1.0;
// The partial derivative of velocity z with respect to state x[10] is 1
H_gps[5*OEKF_N + 10] = 1.0;

// Call the update function (m=6)
oekf_update(&ekf, z_gps, hx_gps, H_gps, R_gps, 6);

// Barometer Observation Model (1D: Altitude)
_float_t z_baro[1] = {baro_altitude}; // Observed value (barometer height, unit consistent with the state)
_float_t hx_baro[1]; // Predict the observed value
_float_t H_baro[1*OEKF_N] = {0}; // 1x14 observation matrix
_float_t R_baro[1*1] = {2.0}; // Altitude observation noise (set according to sensor accuracy)

// Predicted observed value hx_baro: Extract position z from the state vector
hx_baro[0] = ekf->state.p[2]; // 对应x[13]

// Constructing the H matrix: the partial derivative of height with respect to position z is 1
H_baro[0*OEKF_N + 13] = 1.0;

// Call the update function (m=1)
oekf_update(&ekf, z_baro, hx_baro, H_baro, R_baro, 1);

// Optical flow observation model (2D: horizontal velocity x/y)
_float_t z_flow[2] = {flow_vx, flow_vy}; // Observation value (horizontal velocity x/y calculated by optical flow)
_float_t hx_flow[2]; // Predict the observed value
_float_t H_flow[2*OEKF_N] = {0}; // 2x14 observation matrix
_float_t R_flow[2*2] = { // Optical flow velocity noise matrix
    0.05, 0,
    0, 0.05
};

// Predicted observation value hx_flow: Extract horizontal velocity from the state vector
hx_flow[0] = ekf->state.v[0]; // Velocity x (corresponding to x[8])
hx_flow[1] = ekf->state.v[1]; // Velocity y (corresponding to x[9])

// Constructing the H matrix: The partial derivative of the horizontal velocity with respect to the state velocity x/y is 1
H_flow[0*OEKF_N + 8] = 1.0;  // Partial derivative of velocity x with respect to x[8]
H_flow[1*OEKF_N + 9] = 1.0;  // The partial derivative of velocity y with respect to x[9]

// Call the update function (m=2)
oekf_update(&ekf, z_flow, hx_flow, H_flow, R_flow, 2);


Initialize the 14-dimensional state (octonion attitude, velocity, position).
Process multi-modal sensor data such as IMU (acceleration, angular velocity), GPS (position, velocity), and barometer (altitude).
Call oekf_predict_async to handle timestamp misalignment (for example, IMU frequency is 200Hz and GPS frequency is 10Hz).
Design observation models hx and Jacobian matrices H for different sensors (for example, when GPS observes positions, H has non-zero 
values only in the position dimension).

//各传感器的 oekf_update 调用需根据实际数据触发时机执行（例如 GPS 数据频率较低，可通过时间戳判断是否有新数据再更新）。
//The oekf_update calls of each sensor need to be executed according to the actual data triggering timing (for example, the GPS
data frequency is low, and it can be judged whether there is new data through the timestamp before updating).
