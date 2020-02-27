# USRP Radar
This is a C++ program for transmit/receive, processing, display and storing radar data using a USRP software defined radio. It has the following dependencies:
* [UHD library](https://github.com/EttusResearch/uhd)
* [QT](https://www.qt.io/developers/)
* [QWT](http://qwt.sourceforge.net/)
* [Boost](https://www.boost.org/)
* [GNU libgomp](https://gcc.gnu.org/onlinedocs/libgomp/)
* [OpenCV](https://opencv.org/)
* [CUDA](https://developer.nvidia.com/cuda-downloads)

It has been compiled on a Ubuntu Linux 16.04 LTS platform, but has earlier been compiled on both linux and windows platforms. It has been tested with the following USRP configuration:
* [Ettus X310](https://kb.ettus.com/X300/X310)
* [UBX-160](https://www.ettus.com/product/details/UBX160)/[SBX-120](https://www.ettus.com/product/details/SBX120) in port RF0
* [TwinRX](https://www.ettus.com/product/details/TwinRX) in port RF1
* [2x10GB ethernet connection](https://www.ettus.com/product/details/10GIGE-KIT)

One way to compile the program is to open the .pro file in QT Creator and compile through this. The CUDA files is compiled and linked to the program.

RF frontend consists of an amplifier on transmit, and LNA's on receive, combined with filters on transmit and receive. Antenna rig consists of one transmit antenna and two receive antennas in an array. Program uses phase monopulse technique to measure angle of arrival.

For more information on the USRP based radar, see conference article [USRP based cognitive radar testbed](https://publications.ffi.no/handle/123456789/837) from *2017 IEEE Radar Conference (RadarConf)* in Seattle and magazine article [Development and Calibration of a Low-Cost Radar Testbed Based on the Universal Software Radio Periphera](https://www.researchgate.net/publication/338650084_Development_and_Calibration_of_a_Low-Cost_Radar_Testbed_Based_on_the_Universal_Software_Radio_Peripheral) from *IEEE Aerospace and Electronic Systems Magazine*, and publication using the USRP radar for experimental testing [Fully adaptive radar for track update-interval control](https://doi.org/10.1109/RADAR.2018.8378592) from *2018 IEEE Radar Conference (RadarConf18)* in Oklahoma City.


Jonas Myhre Christiansen <jonas-myhre.christiansen@ffi.no>
