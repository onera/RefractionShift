## RefractionShift
RefractionShift is a python module for computing the lateral shift due to atmospheric refraction.

## Installation
First, to install AstroAtmosphere (Python package), run: 
```bash
pip install AstroAtmosphere
```
then run:
```bash
pip install RefractionShift 
```

## Usage
This module allow you to calculate the lateral shift with four different methods. First, using numerical integration along the optical path of the light ray. We have employed the two-layer model of the atmosphere, with a constant temperature gradient in the troposphere and nil beyond. Second, it is possible to use three approximations of different order depending on the refractive index and the Earth's roundness.
Details on the use of the four methods are given in the file example.py .


## License
[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
