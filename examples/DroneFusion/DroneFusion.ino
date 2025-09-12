#include <Wire.h>
#include "tinyoekf.h"
#include "tinyoekf_custom.h"

// Sensor data structure
typedef struct {
    // IMU data (acceleration m/s², angular velocity rad/s)
    _float_t accel[3] = {0};
    _float_t gyro[3] = {0};
    // GPS data (position m, speed m/s)
    _float_t gps_pos[3] = {0};
    _float_t gps_vel[3] = {0};
    // Barometer data (altitude m)
    _float_t baro_alt = 0;
    // Timestamp (microseconds)
    uint64_t timestamp = 0;
} SensorData;

// Global variable
oekf_t ekf;                  // EKF example
SensorData sensor_data;      // Sensor data
_float_t Q[OEKF_N*OEKF_N] = {0};  // Process noise matrix
_float_t R_gps[6*6] = {0};   // GPS observation noise matrix
_float_t R_baro[1*1] = {0};  // Barometer observation noise matrix

// Initialize EKF parameters
void init_ekf() {
    // Initialize the state covariance matrix (the diagonal elements are the initial uncertainties)
    _float_t pdiag[OEKF_N] = {
        0.1f,  // Real part of octonion
        0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f, 0.1f,  // Imaginary part of octonion
        0.5f, 0.5f, 0.5f,  // speed
        1.0f, 1.0f, 1.0f   // position
    };
    oekf_initialize(&ekf, pdiag);

    // Initialization process noise (motion uncertainty)
    for (int i = 0; i < OEKF_N; i++) {
        Q[i*OEKF_N + i] = (i < 8) ? 1e-4f : 1e-3f;  // Attitude noise < velocity/position noise
    }

    // Initialize GPS noise (position noise > velocity noise)
    for (int i = 0; i < 3; i++) {
        R_gps[i*6 + i] = 1.0f;       // Position noise
        R_gps[(i+3)*6 + (i+3)] = 0.5f;  // speed noise
    }

    // Initialize barometer noise
    R_baro[0] = 0.2f;
}

// Simulated sensor data (to be replaced with hardware reading in practical applications)
void read_sensors() {
    // Generate simulated data (sine motion example)
    static uint32_t t = 0;
    t++;
    sensor_data.timestamp = micros();

    // Simulated acceleration (including noise)
    sensor_data.accel[0] = sin(t*0.01f) * 2.0f + (random()%100 - 50)*0.01f;
    sensor_data.accel[1] = cos(t*0.01f) * 2.0f + (random()%100 - 50)*0.01f;
    sensor_data.accel[2] = 9.81f + (random()%100 - 50)*0.02f;

    // Simulated angular velocity (including noise)
    sensor_data.gyro[0] = (random()%100 - 50)*0.001f;
    sensor_data.gyro[1] = (random()%100 - 50)*0.001f;
    sensor_data.gyro[2] = (random()%100 - 50)*0.001f;

    // Simulated GPS data (including noise)
    sensor_data.gps_pos[0] = t*0.01f + (random()%100 - 50)*0.05f;
    sensor_data.gps_pos[1] = sin(t*0.01f) + (random()%100 - 50)*0.05f;
    sensor_data.gps_pos[2] = cos(t*0.01f) + (random()%100 - 50)*0.05f;
    sensor_data.gps_vel[0] = 0.01f + (random()%100 - 50)*0.01f;
    sensor_data.gps_vel[1] = cos(t*0.01f)*0.01f + (random()%100 - 50)*0.01f;
    sensor_data.gps_vel[2] = -sin(t*0.01f)*0.01f + (random()%100 - 50)*0.01f;

    // Simulated barometer height (including noise)
    sensor_data.baro_alt = sensor_data.gps_pos[2] + (random()%100 - 50)*0.02f;
}

// IMU observation update
void update_imu() {
    _float_t z_imu[3] = {sensor_data.accel[0], sensor_data.accel[1], sensor_data.accel[2]};
    _float_t hx_imu[3];
    _float_t H_imu[3*OEKF_N] = {0};

    // Predicted observed value (attitude rotation + coupling term correction)
    octonion_rotate(&ekf.state.q, sensor_data.accel, hx_imu);
    hx_imu[0] += ekf.state.q.i[5] * 0.1f;  // Sequential deviation correction of the X-axis
    hx_imu[1] += ekf.state.q.i[5] * 0.1f;  // Correct the Y-axis for sequential deviation
    hx_imu[2] += ekf.state.q.i[6] * 0.05f; // Perturbation coupling correction Z-axis

    // Construct the observation matrix H (supplement the partial derivatives of the coupling terms)
    H_imu[0*OEKF_N + 6] = 0.1f;  // Sensitivity of the X-axis to i[5]
    H_imu[1*OEKF_N + 6] = 0.1f;  // Sensitivity of the Y-axis to i[5]
    H_imu[2*OEKF_N + 7] = 0.05f; // Sensitivity of the Z-axis to i[6]

    // Call lightweight update
    oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);
}

// GPS observation update
void update_gps() {
    _float_t z_gps[6] = {
        sensor_data.gps_pos[0], sensor_data.gps_pos[1], sensor_data.gps_pos[2],
        sensor_data.gps_vel[0], sensor_data.gps_vel[1], sensor_data.gps_vel[2]
    };
    _float_t hx_gps[6];
    _float_t H_gps[6*OEKF_N] = {0};

    // Predicted observed values (directly taken from the state)
    memcpy(hx_gps, &ekf.state.p, 3*sizeof(_float_t));  // position
    memcpy(&hx_gps[3], &ekf.state.v, 3*sizeof(_float_t));  // speed

    // Construct the observation matrix H (only the position and velocity dimensions have non-zero values)
    for (int i = 0; i < 3; i++) {
        H_gps[i*OEKF_N + (11 + i)] = 1.0f;  // Position corresponding to states 11-13
        H_gps[(i+3)*OEKF_N + (8 + i)] = 1.0f;  // Speed corresponds to states 8-10
    }

    // Call update
    oekf_update(&ekf, z_gps, hx_gps, H_gps, R_gps, 6);
}

// Barometer observation update
void update_baro() {
    _float_t z_baro[1] = {sensor_data.baro_alt};
    _float_t hx_baro[1] = {ekf.state.p[2]};  // The predicted height is the Z-axis position
    _float_t H_baro[1*OEKF_N] = {0};
    H_baro[0*OEKF_N + 13] = 1.0f;  // Height corresponds to state 13 (Z-axis position)

    // Call update
    oekf_update(&ekf, z_baro, hx_baro, H_baro, R_baro, 1);
}

void setup() {
    Serial.begin(115200);
    Wire.begin();
    init_ekf();
    randomSeed(analogRead(0));  // Initialize the random number seed
}

void loop() {
    // Read sensor data
    read_sensors();

    // Asynchronous prediction (handling timestamp misalignment)
    oekf_predict_async(&ekf, sensor_data.timestamp, Q, sensor_data.accel, sensor_data.gyro);

    // Sensor update (scheduled by frequency, simplified here as updating every cycle)
    update_imu();         // IMU high-frequency update (200Hz simulation)
    if (millis() % 100 < 10) update_gps();  // GPS low-frequency update (10Hz simulation)
    if (millis() % 50 < 10) update_baro();  // Barometer medium-frequency update (20Hz simulation)

    // Output the fusion result
    Serial.print("State: ");
    Serial.print("q=[");
    Serial.print(ekf.state.q.r, 4);
    for (int i = 0; i < 7; i++) {
        Serial.print(",");
        Serial.print(ekf.state.q.i[i], 4);
    }
    Serial.print("], v=[");
    for (int i = 0; i < 3; i++) {
        if (i > 0) Serial.print(",");
        Serial.print(ekf.state.v[i], 4);
    }
    Serial.print("], p=[");
    for (int i = 0; i < 3; i++) {
        if (i > 0) Serial.print(",");
        Serial.print(ekf.state.p[i], 4);
    }
    Serial.println("]");

    delay(5);  // Control the loop frequency (approximately 200Hz)
}
