// Assume that the IMU observation is 3-dimensional acceleration假设IMU观测为3维加速度
_float_t z_imu[3] = {accel_x, accel_y, accel_z};  // 观测值 Observation value
_float_t hx_imu[3];  // 预测观测值（需根据模型计算）Predicted observed value (needs to be calculated based on the model)
_float_t H_imu[3*OEKF_N];  // 3x14观测矩阵 3x14 observation matrix
_float_t R_imu[3*3] = {0.1,0,0, 0,0.1,0, 0,0,0.1};  // 对角噪声矩阵 diagonal noise matrix

// 调用轻量化更新（m=3） Call lightweight update (m=3)
oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);

/ octonion_rotate(&ekf->state.q, accel_body, hx_imu);

// 新增耦合项后的观测预测：
_float_t hx_imu[3];
octonion_rotate(&ekf->state.q, accel_body, hx_imu);
// 加入i[5]（顺序偏差量）和i[6]（扰动耦合）对观测的修正
hx_imu[0] += ekf->state.q.i[5] * 0.1;  // 顺序偏差修正X轴
hx_imu[1] += ekf->state.q.i[5] * 0.1;  // 顺序偏差修正Y轴
hx_imu[2] += ekf->state.q.i[6] * 0.05;  // 扰动耦合修正Z轴

// 构建H矩阵（补充耦合项的偏导数）：
_float_t H_imu[3*OEKF_N] = {0};
// 原有姿态和速度相关偏导数（略）...
// 新增i[5]和i[6]的偏导数（对应x[6]和x[7]）
H_imu[0*OEKF_N + 6] = 0.1;  // X轴观测对i[5]的敏感度
H_imu[1*OEKF_N + 6] = 0.1;  // Y轴观测对i[5]的敏感度
H_imu[2*OEKF_N + 7] = 0.05;  // Z轴观测对i[6]的敏感度

初始化 14 维状态（八元数姿态、速度、位置）。
处理 IMU（加速度、角速度）、GPS（位置、速度）、气压计（高度）等多模态传感器数据。
调用 oekf_predict_async 处理时间戳错位（如 IMU 频率 200Hz，GPS 频率 10Hz）。
针对不同传感器设计观测模型 hx 和雅各比矩阵 H（如 GPS 观测位置时，H 仅在位置维度有非零值）。
