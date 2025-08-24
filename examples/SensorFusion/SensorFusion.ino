/* SensorFusion: Sensor fusion on Arduino using TinyOEKF.  
 *
 * Copyright (C) 2015 Simon D. Levy
 *Modifications Copyright (C) 2025 ZuoCen Liu
 * MIT License
 */


#include <tinyoekf.h>
#include <SFE_BMP180.h>
#include <Wire.h>

static const uint8_t LM35_PIN = A0;

// Original sensor object
SFE_BMP180 bmp;

// OEKF object and asynchronous parameters
oekf_t oekf;
unsigned long lm35_last_time = 0;  // LM35 timestamp (milliseconds)
unsigned long bmp_last_time = 0;   // BMP180 timestamp (milliseconds)
const float LM35_RATE = 0.01f;     // 100Hz (interval of 0.01 seconds)
const float BMP_RATE = 0.1f;       // 10Hz (interval of 0.1 seconds)

// noise matrix
_float_t Q_lm35[OEKF_N*OEKF_N] = {0};  // LM35 process noise
_float_t Q_bmp[OEKF_N*OEKF_N] = {0};   // BMP process noise
_float_t R_lm35[1*1] = {0.1f};         // LM35 observation noise (temperature)
_float_t R_bmp[2*2] = {0.5f, 0, 0, 0.5f};// BMP observation noise (air pressure, temperature)

// Initialization process noise Q (diagonal matrix)
void init_Q(_float_t Q[OEKF_N*OEKF_N], _float_t val) {
    memset(Q, 0, OEKF_N*OEKF_N*sizeof(_float_t));
    for (int i=0; i<OEKF_N; i++) Q[i*OEKF_N +i] = val;
}

// State constraint function
void apply_constraints(oekf_t *oekf) {
    // Height is non-negative
    if (oekf->state.p[2] < 0) oekf->state.p[2] = 0;
    // Temperature range constraints (such as -40~85°C)
    if (oekf->x[7] < -40) oekf->x[7] = -40;
    else if (oekf->x[7] > 85) oekf->x[7] = 85;
}

static const float EPS = 1e-4;

// So process model Jacobian is identity matrix

// Adapted from https://github.com/sparkfun/BMP180_Breakout
static void getBaroReadings(double & T, double & P) {
    char status = bmp.startTemperature();
    if (status != 0) {
        delay(status);
        status = bmp.getTemperature(T);
        if (status != 0) {
            status = bmp.startPressure(3);
            if (status != 0) {
                delay(status);
                status = bmp.getPressure(P, T);
            }
        }
    }
}

void setup() 
{
    // Use identity matrix as initiali covariance matrix
    const float Pdiag[OEKF_N] = {1, 1};

    Serial.begin(115200);
     if (!bmp.begin()) {
        Serial.println("Could not find BMP180!");
        while (1);
    }

    // Start reading from baro
    baro.begin();

    // Set up to read from LM35
    analogReference(INTERNAL);
    
     // Initialize the OEKF state covariance
    _float_t pdiag[OEKF_N] = {
        1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,  // Octonion (attitude)
        0.1f,0.1f,0.1f,                            // speed
        1.0f,1.0f,1.0f                             // Position (mainly height)
    };
    oekf_initialize(&oekf, pdiag);
    
    // Initialize the noise matrix
    init_Q(Q_lm35, 1e-4f);  // The process noise of LM35 is relatively small.
    init_Q(Q_bmp, 1e-3f);   // The BMP process noise is slightly larger

    // Initialize the temperature state (reuse x[7])
    oekf.x[7] = 25.0f;  // The initial temperature is assumed to be 25°C
}

void loop() {
     float now = (float)millis() / 1000.0f;  // Current time (in seconds, floating-point type)

    // Processing LM35 (high-frequency temperature sensor): prediction + update
  if (now - lm35_last_time >= LM35_RATE) {
    // Read the LM35 temperature
    float lm35_temp = (analogRead(LM35_PIN) * 0.0048828125) * 100.0;

    // 1. Asynchronous prediction
    oekf_predict_async(&oekf, LM35_RATE, Q_lm35);

    // 2. Construct observations (1-dimensional: temperature)
      // Reuse the imaginary part i6 (x[7]) of the octonion to store the temperature
    oekf.x[7] = 25.0f;  // Reuse the initial temperature in setup()
      _float_t z_lm35[1] = {lm35_temp};
    _float_t hx_lm35[1] = {oekf.x[7]};  // Predict temperature (x[7])
    
      // 3. Observation matrix H (1x14): only associated with the temperature state
    _float_t H_lm35[1*OEKF_N] = {0};
    H_lm35[0*OEKF_N +7] = 1;

    // 4. Update (1-dimensional observation)
    ekf_custom_scalar_update(&oekf, z_lm35[0], hx_lm35[0], H_lm35, R_lm35[0]);

     //Application Constraints
        apply_constraints(&oekf);

      lm35_last_time = now;
}

    // Processing BMP180 (low-frequency air pressure + temperature): prediction + update
    if (now - bmp_last_time >= BMP_RATE) {
        
       double T_bmp, P_bmp;
      getBaroReadings(T_bmp, P_bmp);  
      float bmp_pressure = P_bmp / 100.0f;  
       float bmp_temp = (float)T_bmp;

        // 1. Asynchronous prediction
        oekf_predict_async(&oekf, BMP_RATE, Q_bmp);

        // 2. Barometric pressure update (associated with position z)
        _float_t h_p[1] = {oekf.state.p[2]};
        _float_t H_p[OEKF_N] = {0};
       H_p[13] = 1;  // Position z corresponds to x[13]
       ekf_custom_scalar_update(&oekf, z[0], h_p[0], H_p, R_bmp[0]);

       // 3. Temperature update (associated with x[7])
       _float_t h_t[1] = {oekf.x[7]};
       _float_t H_t[OEKF_N] = {0};
       H_t[7] = 1;
       ekf_custom_scalar_update(&oekf, z[1], h_t[0], H_t, R_bmp[3]);  // R_bmp[3] is temperature noise

        // 4. Observation matrix H (2x14)
        _float_t H[2*OEKF_N] = {0};
        H[0*OEKF_N +13] = 1;  // Position z (x[13])
        H[1*OEKF_N +7] = 1;   // Temperature (x[7])

        // 5. Disturbance detection (Threshold: air pressure 3hPa, temperature 2°C)
        _float_t res_p = z[0] - hx[0];
        _float_t res_t = z[1] - hx[1];
        bool perturbed = (fabs(res_p) > 3.0f) || (fabs(res_t) > 2.0f);
         oekf_adapt_Q(Q_bmp, perturbed, 2.0f);  // Amplify Q during disturbance
        
        // 6. Update （Note: It is necessary to ensure that OEKF_M >= 2, or modify it to a custom update.）
        oekf_update(&oekf, z, hx, H, R_bmp);

        // Application Constraints
        apply_constraints(&oekf);

        // Output the fusion result
        Serial.print("Fused Altitude: ");
        Serial.print(oekf.state.p[2]);
        Serial.print(" m, Fused Temp: ");
        Serial.println(oekf.x[7]);
        Serial.print(" BMP180Press:");
        Serial.print(bmp_pressure);
        Serial.print(" BMP180Temp:");
        Serial.print(bmp_temp);
        Serial.print(" LM35Temp:");
        Serial.println(lm35_temp);  // It is necessary to ensure that lm35_temp is within the scope (it can be defined as a global variable)
        
        bmp_last_time = now;
    }
    }
        
    // Measurement function simplifies the relationship between state
    // and sensor readings for convenience.  A more realistic
    // measurement function would distinguish between state value and
    // measured value; e.g.:
    //   hx[0] = pow(this->x[0], 1.03);
    //   hx[1] = 1.005 * this->x[1];
    //   hx[2] = .9987 * this->x[1] + .001;


// Extended explanation: If adding MPU6050 (IMU, 200Hz) and GPS (1Hz), the following logic can be referred to:
// - When the IMU is triggered, call oekf_predict_async(dt_imu, Q_imu)
// - When the GPS is triggered, call oekf_predict_async(dt_gps, Q_gps) and then perform position update

// 扩展说明：若添加MPU6050（IMU，200Hz）和GPS（1Hz），可参照以下逻辑：
// - IMU触发时调用oekf_predict_async(dt_imu, Q_imu)
// - GPS触发时调用oekf_predict_async(dt_gps, Q_gps)后执行位置更新
