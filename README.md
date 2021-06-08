## RefractionShift
RefractionShift is a python module for computing the lateral shift due to atmospheric refraction.

## Installation
Check first that you have the Python package AstroAtmosphere, otherwise just run: 
```bash
pip install AstroAtmosphere
```
then run:
```bash
pip install RefractionShift 
python "test_refraction.py"
```
The file `example.py` contains a working example.

## Usage
This module allows you to calculate the lateral shift with four different methods. First, using numerical integration along the optical path of the light ray. To do so, we employ the two-layer model of the atmosphere, with a constant temperature gradient in the troposphere and nil beyond. In addition, there are three approximations of different order depending on the refractive index and Earth's roundness.
Details on the use of the four methods are given in the file `example.py`.


## License
[GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/)
