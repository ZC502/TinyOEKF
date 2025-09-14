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

TinyOEKF is a lightweight, header-only C/C++ implementation of the Extended Kalman Filter (EKF) optimized for high-speed UAV maneuvering. It uses octonions to model attitude-motion coupling and supports asynchronous multi-sensor fusion, making it suitable for microcontrollers like Arduino Uno/Teensy.

## Key Features  
- **Octonion-Based Modeling**: 7 imaginary parts encode rotation-translation sequence, disturbance coupling, and other physical relationships, improving dynamic accuracy.  
- **Asynchronous Fusion**: Handles multi-rate sensors (e.g., 200Hz IMU + 10Hz GPS) via timestamp alignment.  
- **Lightweight**: Pure C implementation with static memory allocation (no `malloc`/`new`), ideal for embedded systems.  


## Basic Examples  
### 1. Sensor Fusion (BMP180 + LM35)  
Located in `examples/SensorFusion`, this example fuses data from a BMP180 barometer and LM35 temperature sensor to estimate altitude.  
- **Hardware**: Connect BMP180 to I2C (SDA/SCL) and LM35 to analog pin A0.  
- **Usage**: Upload `SensorFusion.ino` to an Arduino and monitor output via Serial Monitor (115200 baud).  
### 2. Python Prototyping  
The `python` folder includes:  
- `altitude_fuser.py`: Fuses barometer and sonar data using a Python EKF prototype (requires OpenCV).  
- `mouse_tracker.py`: Demonstrates EKF-based mouse position tracking with OpenCV.  


## Drone 4D Spatiotemporal Navigation & Multi-Sensor Fusion  
TinyOEKF excels in high-speed UAV scenarios, modeling attitude-motion coupling via octonion non-commutativity and enabling high-precision 3D space + time navigation.  


### Core Advantages  
- **Lightweight**: Pure C with no dynamic memory allocation, compatible with Arduino Uno/Teensy.  
- **Coupling Modeling**: 7 imaginary parts encode rotation-translation sequence, disturbance correlation, and other physical couplings, enhancing dynamic accuracy.  
- **Multi-Sensor Support**: Fuses IMU, GPS, barometer, and optical flow (asynchronous, auto-aligned by timestamp).  


### Hardware Preparation  
| Sensor               | Interface   | Recommended Model       |  
|----------------------|-------------|-------------------------|  
| IMU (Accel + Gyro)   | I2C         | MPU9250                 |  
| GPS                  | UART        | NEO-M8N (10Hz)          |  
| Barometer            | I2C         | BMP280                  |  
| (Optional) Optical Flow | I2C/UART | PX4Flow                 |  


### Quick Start  
1. Copy the `TinyOEKF` folder to your Arduino libraries directory.  
2. Open the example: `File > Examples > TinyOEKF > DroneFusion`.  
3. Modify pin definitions in `DroneFusion.ino` according to your hardware (see code comments).  
4. Upload to your flight controller and check output via Serial Monitor (115200 baud).  

### Validation Guide  
1. **Static Test**: When the drone is stationary, position `p[0]-p[2]` should be stable (fluctuation < 0.1m), and `i[5]-i[6]` (coupling terms) should be near 0.  
2. **Dynamic Test**:  
   - Rotation: `i[5]` (sequence bias) should increase with rotation speed.  
   - Translation: `i[3]-i[4]` (rotation-translation coupling) should change significantly during rapid movement.  
3. **Spatiotemporal Consistency**: Fusion trajectory should deviate from GPS trajectory by <1m at high speeds.  

### Parameter Tuning  
- If trajectory diverges: Increase noise in attitude dimensions of `Q` (allows larger corrections).  
- If sensitive to disturbances: Reduce `i[6]` coupling coefficient (e.g., 0.05 → 0.03 in examples).  
- If GPS jumps: Increase position noise in `R_gps` (e.g., 1e1 → 5e1).  

## Installation  
- **Arduino**: Copy the `TinyOEKF` folder to `Documents/Arduino/libraries`.  
- **C/C++**: Include `tinyoekf.h` and `tinyoekf_custom.h` in your project.  

## License  
MIT License (see `LICENSE.md` for details).  

## Acknowledgments  
Based on TinyEKF by Simon D. Levy, modified for octonion-based UAV navigation by ZuoCen Liu.
