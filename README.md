# TinyOEKF: Lightweight Octonion Extended Kalman Filter for UAVs

# TinyEKF: Lightweight C/C++ Extended Kalman Filter with Python for prototyping

<a href="examples/SensorFusion/SensorFusion.ino"><img src="media/barotemp2.jpg" width=750></a>
<a href="python/altitude_fuser.py"><img src="media/altitude.png" width=1000></a>

TinyEKF is a simple, header-only C/C++ implementation of the
[Extended Kalman Filter](https://simondlevy.github.io/ekf-tutorial/) 
that is general enough to use on different projects.  It supports both single- and
double-precision floating-point computation.  In order to make it practical for
running on Arduino, STM32, and other microcontrollers, it uses static
(compile-time) memory allocation (no "new" or "malloc").  The **examples**
folder includes both a "pure C" example from the literature, as well as an
Arduino example of sensor fusion.  The **python** folder includes a Python
class that you can use to prototype your EKF before implementing it in C or C++.

Arduino users can simply install or drag the whole TinyEKF folder into their Arduino libraries folder. 
The **examples/SensorFusion** folder contains a little sensor fusion example using a 
[BMP180 barometer](https://www.sparkfun.com/products/11824) and 
[LM35 temperature sensor](http://www.robotshop.com/en/dfrobot-lm35-linear-temperature-sensor.html).
I have run this example on an Arduino Uno and a Teensy 3.2. The BMP180, being an I^2C sensor, should be connected
to pins 4 (SDA) and 5 (SCL) of the Uno, or pins 18 (SDA) and 19 (SCL) of the Teensy.  For other Arduino boards,
consult the [documentation](https://www.arduino.cc/en/Reference/Wire) on the Wire library. The analog output
from the LM35 should go to the A0 pin of your Arduino or Teensy.

In addition to the class definition, the **python** folder has an example of mouse tracking, using OpenCV. 
So you will have to install OpenCV to run this example. There is also a sensor-fusion example in this folder.


# TinyOEKF: Lightweight Octonion Extended Kalman Filter for UAVs

TinyOEKF is an enhanced lightweight C/C++ implementation of the Octonion Extended Kalman Filter, evolved from TinyEKF. It is specifically optimized for high-speed UAV maneuvering scenarios, supporting non-commutative attitude modeling, coupled motion prediction, and asynchronous multi-sensor fusion. Designed for embedded systems like Arduino and STM32, it uses static memory allocation (no dynamic allocation) to ensure stability on resource-constrained devices.
1. Installation
1.1 Arduino/Embedded Platforms
Download or clone the repository
Drag the entire TinyOEKF folder into your Arduino libraries directory (typically Documents/Arduino/libraries/)
Verify installation by opening File > Examples > TinyOEKF > SensorFusion in the Arduino IDE
1.2 Python Prototyping
For rapid algorithm validation, the Python implementation requires:
bash
pip install numpy opencv-python  # Dependencies
cd python
python setup.py install          # Install TinyOEKF Python package
2. Core Features
Lightweight: Pure C implementation with no dynamic memory allocation, compatible with Arduino Uno/Teensy/STM32 microcontrollers
Coupled Modeling: 7 octonion imaginary parts encode rotation-translation sequence, disturbance coupling, and other physical correlations, improving accuracy in dynamic scenarios
Multi-Sensor Fusion: Supports asynchronous sensors (IMU, GPS, barometer, optical flow) with automatic timestamp alignment (e.g., 200Hz IMU + 10Hz GPS)
Dual Precision: Supports both single- and double-precision floating-point computation
Python Prototype: Includes Python class for rapid algorithm iteration before embedded deployment
3. Technical Principle: Octonion & 4D Spacetime Algebra
The connection between octonion non-associativity, time dimension, and 4D spacetime is essentially realized through "implicit encoding of physical processes via algebraic structure" — even without explicit time variables, octonion operations capture time irreversibility through non-associativity, while their high-dimensional structure maps 4D spacetime dynamic coupling.
3.1 Non-Associativity: "Algebraic Mirror" of Time Order
Octonions' core property is non-associativity ((a*b)*c ≠ a*(b*c)), which naturally matches the "irreversibility" of time (physical process order cannot be swapped).
3.1.1 Physical Significance
In UAV motion, "rotate A-axis → rotate B-axis → translate" produces different results from "rotate B-axis → rotate A-axis → translate". This order-induced difference is quantifiable via octonion multiplication non-associativity. For example, the difference between (q1*q2)*q3 and q1*(q2*q3) (where q1,q2,q3 are octonions for three operations) directly encodes the physical deviation from step order changes.
3.1.2 Code Implementation
In octonion_mult (defined in src/tinyoekf_custom.h), imaginary part cross terms (e.g., i[3] includes (a->i[0]*b->i[5] - a->i[5]*b->i[0]) inherently contain "operation order" information. Cumulative non-associativity in multiple multiplications (e.g., incremental octonion dq multiplied with current attitude q in prediction steps) essentially models "time-sequential operation order" algebraically.
3.2 Implicit Time Embedding: Coupled via Dynamic Strength
Though no explicit time variable is defined, time is implicitly embedded in octonion updates through physical process "rate" and "timing":
3.2.1 Time Rate Influence
In octonion_update_coupling_strength, angular velocity magnitude omega_norm (essentially rotation rate dθ/dt) scales the real part r — higher motion intensity (faster time rate) reduces r, increasing imaginary part coupling weights. This indirectly encodes "time scale" in the real part: time's impact on coupling is amplified during rapid motion.
3.2.2 Time Interval Accumulation
In the prediction step oekf_predict, time interval dt directly participates in attitude updates via incremental octonion dq = cos(θ/2) + i*sin(θ/2) (where θ = omega * dt). Cumulative dt across multiple predictions superimposes octonion non-associativity effects, essentially converting "time passage" into "cumulative algebraic deviation".
3.3 Mapping to 4D Spacetime: High-Dimensional Imaginary Parts as "Bridges"
4D spacetime (3D space + 1D time) is characterized by "spacetime inseparability" (e.g., relativistic spacetime co-evolution). Octonions' 8-dimensional structure (1 real + 7 imaginary parts) maps this coupling as follows:
3.3.1 Real Part & Time Correlation
The real part r acts as both a coupling strength benchmark and an implicit "time basis vector". For example:
When the UAV is stationary (time passes but space motion is zero), r≈1 and imaginary parts are near zero, reflecting "uniform time flow"
During motion, r decreases and imaginary parts (space coupling terms) strengthen, corresponding to enhanced "spacetime interaction"
3.3.2 Imaginary Parts & Spacetime Cross Terms
Low-dimensional imaginary parts i[0]-i[2] correspond to 3D space rotation axes, directly encoding spatial attitude
Mid-to-high-dimensional imaginary parts i[3]-i[6] encode "space-time" cross-coupling:
i[5] (sequential deviation) → cumulative effect of spatial operation order over time
i[6] (disturbance coupling) → impact of sudden time-domain disturbances on spatial state
This mapping, while not a strict 4D decomposition, unifies "spatial motion" and "temporal evolution" through 7 imaginary parts — analogous to "coordinated spacetime event description" in 4D spacetime.
3.4 Summary: "Isomorphism" Between Algebra & Physical Spacetime
Octonion non-associativity abstracts "time order irreversibility", while their 8D structure simulates 4D spacetime "co-evolution" via dynamic adjustment of real (time scale) and imaginary (space/spacetime coupling) parts. Even without explicit time variables, physical time is embedded in octonion operations through "order", "rate", and "interval", making high-dimensional imaginary parts a natural bridge between algebra and spacetime physics.
4. Quick Start Examples
4.1 Arduino Sensor Fusion (Basic)
The examples/SensorFusion folder demonstrates fusion of:
BMP180 barometer (pressure/altitude)
LM35 temperature sensor
Wiring Guide:
BMP180 (I²C): Connect SDA to pin 4 (Uno) / 18 (Teensy), SCL to pin 5 (Uno) / 19 (Teensy)
LM35: Analog output to A0 pin
Usage:
Open examples/SensorFusion/SensorFusion.ino in Arduino IDE
Select your board (e.g., Arduino Uno) and port
Upload and monitor serial output (115200 baud)
4.2 UAV Multi-Sensor Fusion (Advanced)
The examples/DroneFusion folder includes a 14-dimensional state estimation example fusing:
IMU (3-axis acceleration + 3-axis angular velocity)
GPS (3D position + 3D velocity)
Barometer (altitude)
Optical flow (2D horizontal velocity)
Key code snippet for multi-sensor update:
c
运行
// IMU update (3D acceleration)
_float_t z_imu[3] = {accel_x, accel_y, accel_z};
_float_t hx_imu[3];
octonion_rotate(&ekf->state.q, accel_body, hx_imu); // Attitude coupling
oekf_update(&ekf, z_imu, hx_imu, H_imu, R_imu, 3);

// GPS update (6D: position + velocity)
_float_t z_gps[6] = {gps_p_x, gps_p_y, gps_p_z, gps_v_x, gps_v_y, gps_v_z};
oekf_update(&ekf, z_gps, hx_gps, H_gps, R_gps, 6);

// Barometer update (1D altitude)
_float_t z_baro[1] = {baro_altitude};
oekf_update(&ekf, z_baro, hx_baro, H_baro, R_baro, 1);
4.3 Python Prototyping
python/mousetracker.py: Demonstrates real-time mouse tracking using OpenCV
python/altitude_fuser.py: Simulates altitude fusion from barometer and temperature data
Run with:
bash
python python/mousetracker.py  # Requires OpenCV
5. Hardware Requirements
Sensor	Interface	Recommended Model	Purpose
IMU (Accel + Gyro)	I²C	MPU9250	Attitude/velocity sensing
GPS	UART	NEO-M8N (10Hz)	Absolute position/velocity
Barometer	I²C	BMP280	Altitude correction
(Optional) Optical Flow	I²C/UART	PX4Flow	Low-altitude velocity
6. API Reference
Core Structures
oekf_t: Main OEKF structure containing state, covariance matrix, and timestamps
Octonion: Octonion structure with real part r and 7 imaginary parts i[0]-i[6]
OEKF_State: State wrapper with octonion attitude q, velocity v[3], and position p[3]
Key Functions
Function	Purpose
oekf_initialize()	Initialize OEKF with initial covariance
oekf_predict()	Predict state using motion model
oekf_update()	Update state with sensor observations
octonion_mult()	Octonion multiplication (non-associative)
octonion_rotate()	Rotate vector using octonion attitude
7. License
TinyOEKF is released under the MIT License. See LICENSE.md for details.
Copyright (C) 2024 Simon D. Levy
Modifications Copyright (C) 2025 ZuoCen Liu
8. Contributing
Contributions are welcome! Please submit issues or pull requests to the GitHub repository. For questions, contact maintainer ZuoCen Liu at liouzuocen@qq.com.
