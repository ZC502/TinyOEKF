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

## 1. Installation
### 1.1 Arduino/Embedded Platforms
- Download or clone the repository
- Drag the entire TinyOEKF folder into your Arduino libraries directory (typically Documents/Arduino/libraries/)
- Verify installation by opening File > Examples > TinyOEKF > SensorFusion in the Arduino IDE

### 1.2 Python Prototyping
For rapid algorithm validation, the Python implementation requires:

pip install numpy opencv-python  # Dependencies
cd python
python setup.py install          # Install TinyOEKF Python package


## 基础示例  
- **examples/SensorFusion**: 包含BMP180气压计与LM35温度传感器的融合示例...  
- **python文件夹**: 包含OpenCV鼠标跟踪、传感器融合等原型示例...  

## 无人机4维时空导航与多传感器融合  
TinyOEKF特别适用于无人机高速机动场景，通过八元数的非交换性建模姿态-运动耦合关系，支持3D空间+时间的高精度导航。  

### 核心优势  
- **轻量化**：纯C实现，无动态内存分配...  
- **耦合建模**：7个虚部编码旋转-平移顺序...  
- **多传感器融合**：支持IMU、GPS、气压计等异步传感器...  

<!-- 硬件准备、快速开始等补充内容 -->

## 安装与使用  
Arduino用户可将TinyOEKF文件夹复制到Arduino libraries目录...  
```{insert\_element\_0\_}### 1. 补充位置建议  
推荐将“无人机4维时空导航与多传感器融合”章节添加到原`README.md`中**现有示例介绍之后**，作为独立的核心应用章节。具体位置参考：  
- 原`README.md`已包含`examples/SensorFusion`和Python示例（鼠标跟踪等）的介绍，补充内容可紧随其后，形成“基础示例→无人机专项应用”的递进结构，既保留原有内容的完整性，又突出TinyOEKF在无人机场景的核心价值。  


### 2. 关于“四维时空符合相对论协变性”的说明建议  
**不建议补充此内容**，原因如下：  
- 从代码实现来看，“4维时空”实际指“3D空间+时间戳同步”（如`oekf_predict_async`处理传感器异步性），核心是工程上的时空数据对齐，而非严格的相对论协变性（涉及洛伦兹变换、时空弯曲等物理理论）。  
- 若无明确的理论推导或代码实现支撑（如未在状态方程中引入相对论修正项），添加此描述可能导致误解，建议保持“3D空间+时间同步”的工程化表述，更贴合实际功能。  


### 3. 原TinyEKF资料的处理建议  
**不建议删除，但需调整优先级**：  
- 原TinyEKF的资料（如Python鼠标跟踪示例）是项目的历史基础和功能补充，保留可体现兼容性和扩展性。  
- 可通过排版调整突出TinyOEKF的新特性：在介绍示例时，先重点描述`DroneFusion`等无人机相关示例，再简要提及原有`SensorFusion`和Python示例（如“除无人机场景外，项目还包含基础传感器融合、鼠标跟踪等示例，详见对应文件夹”）。  


### 调整后的README结构示例（简化）  
```markdown
# TinyOEKF: Lightweight C/C++ Extended Kalman Filter with Python for prototyping

<!-- 原有简介、图片等 -->

TinyOEKF is a simple, header-only C/C++ implementation of the Extended Kalman Filter...  

## 基础示例  
- **examples/SensorFusion**: 包含BMP180气压计与LM35温度传感器的融合示例...  
- **python文件夹**: 包含OpenCV鼠标跟踪、传感器融合等原型示例...  

## 无人机4维时空导航与多传感器融合  
TinyOEKF特别适用于无人机高速机动场景，通过八元数的非交换性建模姿态-运动耦合关系，支持3D空间+时间的高精度导航。  

### 核心优势  
- **轻量化**：纯C实现，无动态内存分配...  
- **耦合建模**：7个虚部编码旋转-平移顺序...  
- **多传感器融合**：支持IMU、GPS、气压计等异步传感器...  

<!-- 硬件准备、快速开始等补充内容 -->

## 安装与使用  
Arduino用户可将TinyOEKF文件夹复制到Arduino libraries目录...  
