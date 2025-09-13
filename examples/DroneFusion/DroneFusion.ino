/*
 * Example of UAV Multi-Sensor Fusion (based on TinyOEKF)
 * Function: Fuse IMU, GPS, and barometer data to estimate the position, velocity, and attitude of the UAV
 * Hardware Requirements:
 *   - Arduino-compatible boards (such as Teensy 4.1, ESP32)
 *   - MPU9250 (IMU, acceleration + angular velocity)
 *   - GPS module (supports NMEA protocol, such as NEO-M8N)
 *   - BMP280 (barometer)
 * Wiring Instructions:
 *   - MPU9250: SDA -> IMU_SDA, SCL -> IMU_SCL, VCC -> 3.3V, GND -> GND
 *   - GPS: TX -> GPS_RX, RX -> GPS_TX (UART communication), VCC -> 3.3V, GND -> GND
 *   - BMP280: SDA -> IMU_SDA (Can share I2C with IMU), SCL -> IMU_SCL, VCC -> 3.3V, GND -> GND
 */

#include <Wire.h>
#include <TinyOEKF.h>
#include <MPU9250_WE.h>  // It is recommended to use the MPU9250 library with calibration function（https://github.com/wollewald/MPU9250_WE）
#include <Adafruit_BMP280.h>  // BMP280 library
#include <SoftwareSerial.h>   // If GPS uses a software serial port

// -------------------------- Hardware Pin Definition --------------------------
#define IMU_SDA         21    // I2C SDA pin (ESP32 example)
#define IMU_SCL         22    // I2C SCL pin
#define GPS_RX          16    // GPS receiving pin (UART)
#define GPS_TX          17    // GPS transmit pin
#define GPS_BAUDRATE    9600  // GPS baud rate

// -------------------------- Sensor Object --------------------------
MPU9250_WE imu;                  // IMU object
Adafruit_BMP280 baro;            // Barometer object
SoftwareSerial gpsSerial(GPS_RX, GPS_TX);  // GPS software serial port (hardware serial port can be replaced with Serial2, etc.)

// -------------------------- OEKF Configuration --------------------------
oekf_t ekf;                       // Example of OEKF
// Initial covariance matrix (diagonal elements, set according to sensor accuracy)
_float_t P_diag[OEKF_N] = {
  // Octonions (8-dimensional: real part + 7 imaginary parts) - Initial uncertainty is small
  1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3,
  // Velocity (3D) - Initial uncertainty is medium
  1e-1, 1e-1, 1e-1,
  // Position (3D) - Initial uncertainty is large
  1e1, 1e1, 1e1
};
// Process noise covariance matrix Q (a diagonal matrix, reflecting model uncertainty)
_float_t Q[OEKF_N * OEKF_N] = {0};  // Initialize all to 0

// -------------------------- Sensor data cache --------------------------
_float_t accel_body[3] = {0};  // Body coordinate system acceleration (m/s²)
_float_t omega[3] = {0};       // Angular velocity in body coordinate system (rad/s)
_float_t gps_p[3] = {0};       // GPS position (m, northeast-up coordinate system)
_float_t gps_v[3] = {0};       // GPS speed (m/s, Northeast-East-Sky coordinate system)
_float_t baro_h = 0;           // Barometer height (m, relative to sea level)
bool gps_new_data = false;     // GPS new data marker
bool baro_new_data = false;    // Barometer new data flag

// -------------------------- Function declaration --------------------------
void initSensors();                  // Sensor initialization
void updateIMU();                    // Update IMU data
void updateGPS();                    // Update GPS data
void updateBaro();                   // Update barometer data
void gps_hx(const OEKF_State *state, _float_t hx[6]);  // GPS observation model prediction
void gps_H(_float_t H[6*OEKF_N]);    // Construction of GPS observation matrix

void setup() {
  Serial.begin(115200);       // Debug serial port
  Wire.begin(IMU_SDA, IMU_SCL);  // Initialize I2C
  gpsSerial.begin(GPS_BAUDRATE);  // Initialize the GPS serial port

  // Initialize process noise Q (diagonal element setting: adjust according to motion characteristics)
  for (int i = 0; i < OEKF_N; i++) {
    if (i < 8) Q[i*OEKF_N + i] = 1e-6;  // Octonion process noise
    else if (i < 11) Q[i*OEKF_N + i] = 1e-4;  // velocity process noise
    else Q[i*OEKF_N + i] = 1e-2;  // position process noise
  }

  // Initialize the sensor
  initSensors();

  // Initialize OEKF
  oekf_initialize(&ekf, P_diag);
  ekf.last_timestamp = micros();  // Initial timestamp
  Serial.println("System initialization completed");
}

void loop() {
  uint64_t now = micros();  // Current time (microseconds)

  // 1. Update IMU data (high frequency: ~200Hz)
  updateIMU();
  // Asynchronous prediction (calculating dt based on timestamps to handle inconsistent sensor frequencies)
  oekf_predict_async(&ekf, now, Q, accel_body, omega);

  // 2. Update GPS data (low frequency: ~10Hz)
  updateGPS();
  if (gps_new_data) {
    gps_new_data = false;  // Clear the flag
    _float_t z_gps[6] = {gps_p[0], gps_p[1], gps_p[2], gps_v[0], gps_v[1], gps_v[2]};
    _float_t hx_gps[6];
    _float_t H_gps[6*OEKF_N] = {0};  // Initialize the observation matrix

    // Calculate GPS predicted observation values and observation matrix
    gps_hx(&ekf.state, hx_gps);
    gps_H(H_gps);

    // GPS observation noise matrix R (diagonal elements, set according to the GPS datasheet)
    _float_t R_gps[6*6] = {
      1.0, 0, 0, 0, 0, 0,    // Position x noise
      0, 1.0, 0, 0, 0, 0,    // Position y noise
      0, 0, 1.5, 0, 0, 0,    // Position z noise (lower precision)
      0, 0, 0, 0.1, 0, 0,    // Speed x Noise
      0, 0, 0, 0, 0.1, 0,    // Speed y noise
      0, 0, 0, 0, 0, 0.15    // Speed z noise
    };

    // Perform EKF update
    oekf_update(&ekf, z_gps, hx_gps, H_gps, R_gps, 6);
    Serial.println("GPS data update completed");
  }

  // 3. Update barometer data (medium frequency: ~50Hz)
  updateBaro();
  if (baro_new_data) {
    baro_new_data = false;  // Clear the flag
    _float_t z_baro[1] = {baro_h};
    _float_t hx_baro[1];
    _float_t H_baro[1*OEKF_N] = {0};  // Initialize the observation matrix

    // Barometer observation model: Predicted height = position z in the state
    hx_baro[0] = ekf.state.p[2];
    // Observation matrix: The partial derivative of height with respect to position z (x[13]) is 1
    H_baro[0*OEKF_N + 13] = 1.0;

    // Barometer observation noise (set according to sensor accuracy)
    _float_t R_baro[1*1] = {2.0};  // Highly noisy

    // Perform EKF update
    oekf_update(&ekf, z_baro, hx_baro, H_baro, R_baro, 1);
  }

  // 4. Output status (output once every 100ms to avoid serial port blocking)
  static uint64_t last_print = 0;
  if (now - last_print > 100000) {  // 100ms
    last_print = now;
    Serial.printf("t:%llu, p:%.2f,%.2f,%.2f, v:%.2f,%.2f,%.2f, q_i5:%.4f, q_i6:%.4f\n",
                  now / 1000,  // Convert to milliseconds
                  ekf.state.p[0], ekf.state.p[1], ekf.state.p[2],
                  ekf.state.v[0], ekf.state.v[1], ekf.state.v[2],
                  ekf.state.q.i[5], ekf.state.q.i[6]);
  }
}

// -------------------------- Sensor initialization function --------------------------
void initSensors() {
  // 初始化MPU9250
  if (!imu.init()) {
    Serial.println("MPU9250 initialization failed! Please check the wiring.");
    while (1);  // Stop here
  }
  imu.setAccRange(MPU9250_ACC_RANGE_8G);  // Set the acceleration range
  imu.setGyroRange(MPU9250_GYRO_RANGE_500);  // Set the angular velocity range
  imu.calibrateImu();  // Set the angular velocity range
  Serial.println("MPU9250 initialization completed");

  // Initialize BMP280
  if (!baro.begin(0x76)) {  // 0x76 is a common I2C address (another possible one is 0x77)
    Serial.println("BMP280 initialization failed! Please check the wiring");
    while (1);
  }
  baro.setSampling(Adafruit_BMP280::MODE_NORMAL,  // Normal mode
                   Adafruit_BMP280::SAMPLING_X2,  // Temperature sampling
                   Adafruit_BMP280::SAMPLING_X16, // Pressure sampling
                   Adafruit_BMP280::FILTER_X16,   // filtering
                   Adafruit_BMP280::STANDBY_MS_500);  // standby time
  Serial.println("BMP280 initialization completed");

  Serial.println("All sensors have been initialized, waiting for data...");
}

// -------------------------- Sensor data update function --------------------------
void updateIMU() {
  // Read IMU data (converted to the body coordinate system)
  imuData_t data = imu.getImuData();
  // Acceleration: Convert to m/s² (the original data unit is g, 1g ≈ 9.81m/s²)
  accel_body[0] = data.ax * 9.81;
  accel_body[1] = data.ay * 9.81;
  accel_body[2] = data.az * 9.81;
  // Angular velocity: Convert to rad/s (the unit of the original data is °/s, conversion formula: °/s * π/180)
  omega[0] = data.gx * DEG_TO_RAD;
  omega[1] = data.gy * DEG_TO_RAD;
  omega[2] = data.gz * DEG_TO_RAD;
}

void updateGPS() {
  // Simple NMEA parsing (it is recommended to use a dedicated GPS library such as TinyGPSPlus in actual projects)
  while (gpsSerial.available() > 0) {
    char c = gpsSerial.read();
    // Assume that the position and velocity are updated after parsing the GGA or RMC frame (this is an example logic here)
    // In practical applications, it needs to be replaced with real parsing code (such as extracting latitude and longitude and converting them to meters)
    static bool parsed = false;
    if (parsed) {  // Simulation parsing successful
      // This is just an example; in practice, it needs to be converted based on GPS data (such as from WGS84 to the northeast-up coordinate system)
      gps_p[0] += 0.1;  // Simulate the position increment in the x-direction
      gps_p[1] += 0.05; // Simulate the position increment in the y direction
      gps_p[2] += 0.02; // Simulate z-direction position increment
      gps_v[0] = 1.2;   // Simulate the velocity in the x-direction
      gps_v[1] = 0.8;   // Simulate the velocity in the y-direction
      gps_v[2] = 0.3;   // Simulate the velocity in the z-direction
      gps_new_data = true;
      parsed = false;
    }
    // Mark as to be parsed when the GGA frame header is detected
    if (c == '$' && gpsSerial.peek() == 'G' && gpsSerial.peek() == 'G') {
      parsed = true;
    }
  }
}

void updateBaro() {
  static uint32_t last_baro = 0;
  if (millis() - last_baro > 20) {  // 50Hz sampling
    last_baro = millis();
    // Read the air pressure and convert it to altitude (m)
    baro_h = baro.readAltitude(1013.25);  // 1013.25 is the standard atmospheric pressure at sea level
    baro_new_data = true;
  }
}

// -------------------------- GPS observation model function --------------------------
void gps_hx(const OEKF_State *state, _float_t hx[6]) {
  // Predicted observed value = position and velocity in the state (directly extracted)
  hx[0] = state->p[0];  // Position x
  hx[1] = state->p[1];  // Position y
  hx[2] = state->p[2];  // Position z
  hx[3] = state->v[0];  // speed x
  hx[4] = state->v[1];  // speed y
  hx[5] = state->v[2];  // speed z
}

void gps_H(_float_t H[6*OEKF_N]) {
  // Initialize all to 0
  memset(H, 0, 6*OEKF_N * sizeof(_float_t));
  // The partial derivative of position x with respect to state x[11] (p[0]) is 1
  H[0*OEKF_N + 11] = 1.0;
  // The partial derivative of position y with respect to state x[12] (p[1]) is 1
  H[1*OEKF_N + 12] = 1.0;
  // The partial derivative of position z with respect to state x[13] (p[2]) is 1.
  H[2*OEKF_N + 13] = 1.0;
  // The partial derivative of velocity x with respect to state x[8] (v[0]) is 1
  H[3*OEKF_N + 8] = 1.0;
  // The partial derivative of velocity y with respect to state x[9] (v[1]) is 1
  H[4*OEKF_N + 9] = 1.0;
  // The partial derivative of velocity z with respect to state x[10] (v[2]) is 1
  H[5*OEKF_N + 10] = 1.0;
}
