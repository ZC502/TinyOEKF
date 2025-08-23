/* SensorFusion: Sensor fusion on Arduino using TinyOEKF.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 */


// These must be defined before including TinyEKF.h
#define EKF_N 2     // pressure, temperature
#define EKF_M 3     // baro pressure, baro temperature, LM35 temperature

static const uint8_t LM35_PIN = 0;

#include <tinyekf.h>
#include <SFE_BMP180.h>
#include <Wire.h>
#include <tinyoekf.h>

// Original sensor object
SFE_BMP180 bmp;
const int LM35_PIN = A0;

// Added: OEKF object and asynchronous processing parameters
oekf_t oekf;
unsigned long lm35_last_time = 0;  // LM35 (high-frequency temperature) timestamp
unsigned long bmp_last_time = 0;   // BMP180 (low-frequency air pressure/temperature) timestamp
const float LM35_RATE = 0.01;      // LM35 sampling rate: 100Hz (interval of 0.01 seconds)
const float BMP_RATE = 0.1;        // BMP180 sampling rate: 10Hz (interval of 0.1 seconds)

// Added: Noise matrix (adjusted according to sensor characteristics)
_float_t Q_lm35[OEKF_N*OEKF_N] = {0};  // Process noise corresponding to LM35
_float_t Q_bmp[OEKF_N*OEKF_N] = {0};   // The process noise corresponding to BMP180
_float_t R_lm35[1*1] = {0.1};          // LM35 observation noise (temperature)
_float_t R_bmp[2*2] = {0.5, 0, 0, 0.5};// BMP180 observation noise (air pressure, temperature)

static const float EPS = 1e-4;

static const float Q[EKF_N*EKF_N] = {

    EPS, 0,   
    0,   EPS
};

static const float R[EKF_M*EKF_M] = {

    EPS, 0,   0,
    0,   EPS, 0,
    0,   0,   EPS
};

// So process model Jacobian is identity matrix
static const float F[EKF_N*EKF_N] = {
    1, 0,
    0, 1
};

static const float H[EKF_M*EKF_N] = {

    1, 0,
    0, 1,
    0, 1
};

static ekf_t _ekf;

static SFE_BMP180 baro;

// Adapted from https://github.com/sparkfun/BMP180_Breakout
static void getBaroReadings(double & T, double & P)
{
    char status = baro.startTemperature();

    if (status != 0) {
        delay(status);
        status = baro.getTemperature(T);
        if (status != 0) {
            status = baro.startPressure(3);
            if (status != 0) {
                delay(status);
                status = baro.getPressure(P,T);
                if (status == 0)
                    Serial.println("error retrieving pressure measurement");
            }
            else Serial.println("error starting pressure measurement");
        }
        else Serial.println("error retrieving temperature measurement");
    }
    else Serial.println("error starting temperature measurement");
}

void setup() 
{
    // Use identity matrix as initiali covariance matrix
    const float Pdiag[EKF_N] = {1, 1};

    ekf_initialize(&_ekf, Pdiag);

    Serial.begin(115200);
     if (!bmp.begin()) {
        Serial.println("Could not find BMP180!");
        while (1);
    }

    // Start reading from baro
    baro.begin();

    // Set up to read from LM35
    analogReference(INTERNAL);
     // Added: Initialize OEKF
    _float_t pdiag[OEKF_N] = {
        1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,  // Octonion (attitude)
        0.1,0.1,0.1,                               // Speed (can be simplified to 0 here)
        1.0,1.0,1.0                                // Position/Altitude (Core Fusion Quantity)
    };
    oekf_initialize(&oekf, pdiag);

    // Added: Initialize the noise matrix (set values on the diagonal, set non-diagonal elements to 0)
    for (int i=0; i<OEKF_N; i++) {
        Q_lm35[i*OEKF_N + i] = 1e-4;  // Process noise during high-frequency updates of LM35
        Q_bmp[i*OEKF_N + i] = 1e-3;   // The process noise of BMP180 during low-frequency updates (slightly large)
    }
}

void loop() {
     unsigned long now = millis() / 1000.0;  // Current time (seconds)

    // 1. Handling LM35 (high-frequency temperature sensor): prediction only, no update
    if (now - lm35_last_time >= LM35_RATE) {
        // Read the LM35 temperature (analog quantity)
        float lm35_temp = (analogRead(LM35_PIN) * 0.0048828125) * 100.0;

        // Asynchronous prediction: Update the state according to the time interval dt (using Q corresponding to LM35)
        oekf_predict_async(&oekf, LM35_RATE, Q_lm35);

        lm35_last_time = now;
    }

    // 2. Processing BMP180 (low-frequency pressure/temperature sensor): prediction + update
    if (now - bmp_last_time >= BMP_RATE) {
        // Read BMP180 data
        float bmp_pressure = bmp.readPressure() / 100.0;  // hectopascal
        float bmp_temp = bmp.readTemperature();

        // Step 1: Predict up to the current time (using Q corresponding to BMP)
        oekf_predict_async(&oekf, BMP_RATE, Q_bmp);

        // Step 2: Construct the observed value z (air pressure, temperature)
        _float_t z[2] = {bmp_pressure, bmp_temp};

        // Step 3: Calculate the predicted observation hx (derived from the current state)
        _float_t hx[2];
        hx[0] = oekf.state.p[2];  // Assume that altitude corresponds to air pressure (actual mapping is required, simplified here)
        hx[1] = oekf.x[9];        // Assume that state x[9] corresponds to temperature (adjustments need to be made according to the actual definition)

        // Step 4: Construct the observation matrix H (2x14), associating only the relevant states
        _float_t H[2*OEKF_N] = {0};
        H[0*OEKF_N + 13] = 1;  // Pressure-related position z (x[13])
        H[1*OEKF_N + 9] = 1;   // Temperature-related state x[9]

        // Step 5: Detect disturbances (residual thresholds: air pressure 5 hectopascals, temperature 2°C)
        bool perturbed = oekf_detect_perturbation(z, hx, 5.0);  // The threshold takes the maximum value
        oekf_adapt_Q(Q_bmp, perturbed, 2.0);  // When disturbed, Q is amplified by a factor of 2

        // Step 6: Execute the update
        oekf_update(&oekf, z, hx, H, R_bmp);

        bmp_last_time = now;
        
    double baroTemperature, baroPressure;
    getBaroReadings(baroTemperature, baroPressure);

    // Read temperature from LM35
    const float lm35Temperature = analogRead(LM35_PIN) / 9.31;

    // Set the observation vector z
    const float z[EKF_M] = {baroPressure, baroTemperature, lm35Temperature};

    // Process model is f(x) = x
    const float fx[EKF_N] = { _ekf.x[0], _ekf.x[1] };

    // Run the prediction step of the DKF
    ekf_predict(&_ekf, fx, F, Q);

    // Measurement function simplifies the relationship between state
    // and sensor readings for convenience.  A more realistic
    // measurement function would distinguish between state value and
    // measured value; e.g.:
    //   hx[0] = pow(this->x[0], 1.03);
    //   hx[1] = 1.005 * this->x[1];
    //   hx[2] = .9987 * this->x[1] + .001;
    const float hx[EKF_M] = {_ekf.x[0], _ekf.x[1], _ekf.x[1] };

    // Run the update step
    ekf_update(&_ekf, z, hx, H, R);

    // Report measured and predicte/fused values
    Serial.print("BMP180Press:");
    Serial.print(z[0]);
    Serial.print(" ");
    Serial.print(" BMP180Temp:");
    Serial.print(z[1]);
    Serial.print(" LM35Temp:");
    Serial.print(z[2]);
    Serial.print(" EKFPress:");
    Serial.print(_ekf.x[0]);
    Serial.print(" EKFTemp:");
    Serial.println(_ekf.x[1]);
    // Output the fusion result (retain the original logic and modify it to the OEKF state)
    Serial.print("Fused Altitude: ");
    Serial.print(oekf.state.p[2]);
    Serial.print(" m, Fused Temp: ");
    Serial.println(oekf.x[9]);
    }
}

// Extended explanation: If adding MPU6050 (IMU, 200Hz) and GPS (1Hz), the following logic can be referred to:
// - When the IMU is triggered, call oekf_predict_async(dt_imu, Q_imu)
// - When the GPS is triggered, call oekf_predict_async(dt_gps, Q_gps) and then perform position update

// 扩展说明：若添加MPU6050（IMU，200Hz）和GPS（1Hz），可参照以下逻辑：
// - IMU触发时调用oekf_predict_async(dt_imu, Q_imu)
// - GPS触发时调用oekf_predict_async(dt_gps, Q_gps)后执行位置更新
