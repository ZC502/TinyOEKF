# DroneFusion Example explanation

## Functional Overview
This example demonstrates how to use the TinyOEKF library to implement multi-sensor fusion navigation for unmanned aerial vehicles (UAVs). It fuses data from IMU, GPS, barometer, and optical flow sensors through the Quaternion Extended Kalman Filter (OEKF) to achieve high-precision estimation of attitude, velocity, and position, which is suitable for high-speed maneuvering scenarios of UAVs.


## Description of Core Code Functions

### 1. Sensor Observation Model and Update Function
The code implements data fusion through the construction of observation models for different sensors and the calling of oekf_update, with the core logic as follows:

#### IMU (accelerometer) fusion
```c
// Calculation of IMU observations (3-axis acceleration) and predicted values
_float_t z_imu[3] = {accel_x, accel_y, accel_z};  // Original acceleration observation
_float_t hx_imu[3];  // Predicted observation value based on the current state
octonion_rotate(&ekf->state.q, accel_body, hx_imu);  // Attitude rotation converts the body acceleration to the navigation frame
hx_imu[0] += ekf->state.q.i[5] * 0.1;  // Add sequential bias correction
hx_imu[1] += ekf->state.q.i[5] * 0.1;
hx_imu[2] += ekf->state.q.i[6] * 0.05;  // Add disturbance coupling correction

// Construct the observation matrix H_imu (describing the relationship between observations and the derivative of the state)
_float_t H_imu[3*OEKF_N] = {0};
H_imu[0*OEKF_N + 6] = 0.1;  // Sensitivity of X-axis observation to i[5]
H_imu[1*OEKF_N + 6] = 0.1;  // Sensitivity of Y-axis observation to i[5]
H_imu[2*OEKF_N + 7] = 0.05;  // Sensitivity of Z-axis observation to i[6]

// Call EKF update (3D observation)
oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);
```
- Function: Convert IMU acceleration data to the navigation coordinate system, and correct it through the imaginary part of the quaternion (i[5] sequence deviation, i[6] disturbance coupling) to improve the accuracy of acceleration estimation in dynamic scenarios.


#### GPS fusion
```c
// Mapping between GPS observations (3-axis position + 3-axis velocity) and predicted values
_float_t z_gps[6] = {gps_p_x, gps_p_y, gps_p_z, gps_v_x, gps_v_y, gps_v_z};
_float_t hx_gps[6];
hx_gps[0] = ekf->state.p[0];  // Position x (corresponding to state x[11])
hx_gps[1] = ekf->state.p[1];  // Position x (corresponding to state x[11])
hx_gps[2] = ekf->state.p[2];  // Position z (corresponding to state x[13])
hx_gps[3] = ekf->state.v[0];  // Velocity x (corresponding to state x[8])
hx_gps[4] = ekf->state.v[1];  // Velocity y (corresponding to state x[9])
hx_gps[5] = ekf->state.v[2];  // Velocity z (corresponding to state x[10])

// Construct the observation matrix H_gps (only the position and velocity dimensions in the state have non-zero values)
_float_t H_gps[6*OEKF_N] = {0};
H_gps[0*OEKF_N + 11] = 1.0;  // The derivative of position x with respect to x[11] is 1
H_gps[1*OEKF_N + 12] = 1.0;  // The derivative of position y with respect to x[12] is 1.
H_gps[2*OEKF_N + 13] = 1.0;  // The derivative of position z with respect to x[13] is 1
H_gps[3*OEKF_N + 8] = 1.0;   // The derivative of velocity x with respect to x[8] is 1
H_gps[4*OEKF_N + 9] = 1.0;   // The derivative of velocity y with respect to x[9] is 1
H_gps[5*OEKF_N + 10] = 1.0;  // The derivative of velocity z with respect to x[10] is 1

// Call EKF update (6-dimensional observation)
oekf_update(&ekf, z_gps, hx_gps, H_gps, R_gps, 6);
```
- **Function**: Directly map the GPS position and velocity observations to the position (`p[0]-p[2]`) and velocity (`v[0]-v[2]`) dimensions in the EKF state, and clarify the linear relationship between the state and observations through the observation matrix.


#### Barometer fusion
```c
// Mapping between barometer observations (altitude) and predicted values
_float_t z_baro[1] = {baro_altitude};  // Barometer altitude observation
_float_t hx_baro[1] = {ekf->state.p[2]};  // Directly use the position z in the state (x[13]) as the predicted value

// Construct the observation matrix H_baro
_float_t H_baro[1*OEKF_N] = {0};
H_baro[0*OEKF_N + 13] = 1.0;  // The derivative of height with respect to position z (x[13]) is 1

// Call EKF update (1-dimensional observation)
oekf_update(&ekf, z_baro, hx_baro, H_baro, R_baro, 1);
```
- **Function**: Use the altitude data from the barometer to correct the Z-axis position in the state, making up for the insufficient accuracy of GPS in the vertical direction.


#### Optical flow fusion
```c
// Optical flow observations (horizontal velocity) and predicted value mapping
_float_t z_flow[2] = {flow_vx, flow_vy};  // Horizontal velocity of optical flow calculation
_float_t hx_flow[2];
hx_flow[0] = ekf->state.v[0];  // Velocity x (x[8])
hx_flow[1] = ekf->state.v[1];  // Velocity y (x[9])

// Construct the observation matrix H_flow
_float_t H_flow[2*OEKF_N] = {0};
H_flow[0*OEKF_N + 8] = 1.0;  // The derivative of velocity x with respect to x[8] is 1.
H_flow[1*OEKF_N + 9] = 1.0;  // The derivative of velocity y with respect to x[9] is 1.

// Call EKF update (2D observation)
oekf_update(&ekf, z_flow, hx_flow, H_flow, R_flow, 2);
```
- **Function**: In indoor scenarios without GPS, correct the X/Y axis velocities in the state through the horizontal velocity observation of optical flow to improve the low-speed hovering accuracy.


### 2. Asynchronous time alignment function
```c
// Handling inconsistent sensor timestamps (such as IMU at 200Hz and GPS at 10Hz)
oekf_predict_async(&ekf, dt, accel_body, omega);
```
- **Function**: Dynamically adjust the prediction step size (`dt`) based on the timestamp differences of sensor data to ensure that sensor data with different frequencies are aligned in the time dimension and avoid fusion delays.

  
## Sensor data parsing dependency library

### 1. IMU（MPU9250）
- **Dependent Libraries**: `MPU9250_WE` (recommended) or `Adafruit_MPU6050`
- **Installation method**：
  1. Open Arduino IDE and go to `Tools > Manage Libraries`
  2. Search for `MPU9250_WE` (Author: Wolfgang Ewald) and click Install
  3. If using the Adafruit library, search for `Adafruit MPU6050` and install it. Additionally, the dependent `Adafruit BusIO` library needs to be installed.


### 2. GPS（NEO-M8N）
- **Dependency library**：`TinyGPSPlus`
- **Installation method**：
1. Search for `TinyGPSPlus` (Author: Mikal Hart) in the Arduino Library Manager.
2. Click Install. It supports NMEA protocol parsing and is lightweight, suitable for embedded scenarios.


### 3. Barometer（BMP280）
- **dependent library**：`Adafruit BMP280 Library`
- **Installation method**：
  1. Search for `Adafruit BMP280` in the library manager.
  2. The dependent `Adafruit Unified Sensor` library will be installed automatically during the installation.


### 4. Optical Flow (PX4Flow)
- **Dependency Library**: `PX4Flow` (custom or community library)
- **Installation Method**:
1. Clone the repository from GitHub: `git clone https://github.com/PX4/PX4Flow.git`
2. Copy the `PX4Flow` folder to the Arduino libraries directory


## Host computer visualization tool

Using Python's `pyserial` and `matplotlib` to real-time plot the UAV status (position, speed, attitude), the steps are as follows:

### 1. Install dependencies
```bash
pip install pyserial matplotlib numpy
```

### 2. Example code (serial_plotter.py)
```python
import serial
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Configure the serial port (modify according to the actual port)
ser = serial.Serial('COM3', 115200, timeout=0.1)

# Initialize data storage
t = []
pos_x, pos_y, pos_z = [], [], []
vel_x, vel_y, vel_z = [], []

# Initialize the graphics
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
line1, = ax1.plot([], [], 'r-', label='X')
line2, = ax1.plot([], [], 'g-', label='Y')
line3, = ax1.plot([], [], 'b-', label='Z')
ax1.set_ylabel('Position (m)')
ax1.legend()

line4, = ax2.plot([], [], 'r-', label='Vx')
line5, = ax2.plot([], [], 'g-', label='Vy')
line6, = ax2.plot([], [], 'b-', label='Vz')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Velocity (m/s)')
ax2.legend()

def init():
    ax1.set_xlim(0, 10)
    ax1.set_ylim(-5, 5)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(-2, 2)
    return line1, line2, line3, line4, line5, line6

def update(frame):
    # Read serial port data (assuming the format: t,px,py,pz,vx,vy,vz)
    data = ser.readline().decode().strip()
    if data:
        parts = data.split(',')
        if len(parts) == 7:
            try:
                t_val = float(parts[0])
                px = float(parts[1])
                py = float(parts[2])
                pz = float(parts[3])
                vx = float(parts[4])
                vy = float(parts[5])
                vz = float(parts[6])

                t.append(t_val)
                pos_x.append(px)
                pos_y.append(py)
                pos_z.append(pz)
                vel_x.append(vx)
                vel_y.append(vy)
                vel_z.append(vz)

                # Update the drawing data
                line1.set_data(t, pos_x)
                line2.set_data(t, pos_y)
                line3.set_data(t, pos_z)
                line4.set_data(t, vel_x)
                line5.set_data(t, vel_y)
                line6.set_data(t, vel_z)

                # Dynamically adjust the X-axis range
                if t_val > ax1.get_xlim()[1]:
                    ax1.set_xlim(t_val - 10, t_val)
                    ax2.set_xlim(t_val - 10, t_val)
            except:
                pass
    return line1, line2, line3, line4, line5, line6

ani = FuncAnimation(fig, update, init_func=init, interval=50, blit=True)
plt.show()
ser.close()
```

### 3. Usage Method
1. Ensure that the drone flight controller outputs status data through the serial port (format: `timestamp, px, py, pz, vx, vy, vz`; it is necessary to add serial port printing logic in `DroneFusion.ino`)
2. Modify the serial port in the code (such as `COM3` or `/dev/ttyUSB0`)
3. Run the script: `python serial_plotter.py` to view the position and velocity curves in real-time

## Precautions
- The sensor data must be consistent with the coordinate system in the code (recommended: ENU coordinate system, with the X-axis eastward, Y-axis northward, and Z-axis upward)
- If the noise characteristics of the sensor are different, it is necessary to adjust the observation noise matrices such as `R_imu` and `R_gps`
- In high-speed maneuvering scenarios, it is recommended to increase the noise of the attitude coupling terms (`i[5]`, `i[6]`) in the `Q` matrix to improve dynamic response
