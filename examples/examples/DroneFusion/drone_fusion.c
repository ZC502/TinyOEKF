// Assume that the IMU observation is 3-dimensional acceleration假设IMU观测为3维加速度
_float_t z_imu[3] = {accel_x, accel_y, accel_z};  // 观测值 Observation value
_float_t hx_imu[3];  // 预测观测值（需根据模型计算）Predicted observed value (needs to be calculated based on the model)
_float_t H_imu[3*OEKF_N];  // 3x14观测矩阵 3x14 observation matrix
_float_t R_imu[3*3] = {0.1,0,0, 0,0.1,0, 0,0,0.1};  // 对角噪声矩阵 diagonal noise matrix

// 调用轻量化更新（m=3） Call lightweight update (m=3)
oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);
