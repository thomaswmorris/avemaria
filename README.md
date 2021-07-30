# Modeling Auto-Regressive Integrated Atmosphere (MARIA)

MARIA is a python-based package that simulates turbulent atmospheric emission using a auto-regressive gaussian process framework, for applications in observational astronomy. Below: a distribution of turbulent water vapor moves through the field of view of the observer. 

![Watch the video](https://user-images.githubusercontent.com/41275226/117068746-acbf8980-acf9-11eb-8016-64fa01e12a77.mp4)

## Background and Methodology

Atmospheric modeling is an important step in both experiment design and subsequent data analysis for ground-based cosmological telescopes observing the cosmic microwave background (CMB). The next generation of ground-based CMB experiments will be marked by a huge increase in data acquisition, with some telescopes to consist of hundreds of thousands of superconducting polarization-sensitive bolometers sampling the sky. This necessitates new methods of efficiently modeling and simulating atmospheric emission at small angular resolutions, with algorithms than can keep up with the high volume of modern telescopes. 

ARAM simulates layers of turbulent atmospheric emission according to a statistical model derived from observations of the atmosphere in the Atacama Desert, from the [Atacama Cosmology Telescope (ACT)](https://lambda.gsfc.nasa.gov/product/act/) and the [Atacama B-Mode Search (ABS)](https://lambda.gsfc.nasa.gov/product/abs/). It uses a sparse-precision auto-regressive Gaussian process algorithm that allows for both fast simulation of high-resolution atmosphere, as well as the ability to simulate arbitrarily long periods of atmospheric evolution. 


A comprehensive outline of the method and products used is outlined [here](https://github.com/tomachito/aram/blob/main/README.md). 


## Examples and Usage 

To install ARAM with PyPi, run

```console
pip install aram
```
The main tool of the ARAM module is the model object. The default model can be easily intitialized as 

```python
from aram import aram

default_model = aram.model()
```

Different models can be initialized by configurating different aspects of the model.

### Arrays

The array defines the set of detectors that observe the sky, along with the properties of their optics and noise. The array is specified as a dictionary

```python
array_config = {'array_shape' : 'hex',      # shape of array
                  'array_fov' : .8,         # maximum span of array (deg)
                          'n' : 500,        # maximum number of detectors (deg)
                   'aperture' : 5.5,        # primary reflector size (meters)
                     'optics' : 'diff_lim', # type of optical system 
                       'band' : 150e9,      # observing band (in Hz)
                       'pink' : 1e0,        # pink noise scaling (dB per octave)
                      'white' : 1e0}        # white noise scaling (dB per Hz)
```
This supplies an diffraction-limited array with parameters similar to the Atacama Cosmology Telescope. Alternatively, the array can be configured manually by supplying an array of values for each parameter. In this case, the first three parameters are replaced by

```python
array_config={'offset_x' : some_array_of_offsets_x,  # in degrees
              'offset_y' : some_array_of_offsets_y}) # in degrees
```

### Observations

The observation defines the sampling time-ordered pointing
```python
obs_config = {'duration'  : 600,  # duration of observation (sec)
              'samp_freq' : 20,   # sampling frequency, in Hz (deg)
              'center_az' : 0,    # center azimuth of scan (deg)
              'center_el' : 45,   # elevation of scan (deg)
              'az_throw'  : 15,   # half of the scan width (deg)
              'az_speed'  : 1.5}  # scanning speed (deg/s)
```
This defines a 10-minute-long observation at 20 Hz, scanning over 30 degrees at an elevation 45 degrees. The scan is centered due north, and scans at 1.5 degrees per second. 

### Sites

The site determines the positions of celestial sources as the earth rotates under the telescope. For supported telescopes, the site also determines the weather generation that drives atmospheric fluctuations. 
```python
site_config = {'site'  : 'ACT',  # site name
               'time'  : 1.6e9}  # timestamp (unix)
```
Longitude, latitude, and altitude can also be entered manually. 

### Models

The model defines the 

```python
model_config = {'n_layers'     : 8,        # number of layers of atmosphere to simulate
                'weather_type' : 'random', # weather generation scheme 
                'min_height'   : 1000,     # minumum atmosphere height (meters)
                'max_height'   : 5000,     # maximum atmosphere height (meters)
                'res'          : .5,       # atmosphere generation resolution (fraction of beam resolution)
                'atm_rms'      : 50,       # desired total RMS of atmospheric signal
```

## Examples

Passing these dictionaries as arguments produces a customized model

```python
from aram import aram

my_model = aram.model(array_config=array_config,
                      site_config=site_config,
                      obs_config=obs_config,
                      model_config=model_config)
```
Data can then be simulated from the model by running 

```python
data = my_model.simulate(do_atmosphere=True,
                         do_cmb=True,
                         do_noise=True)
```
This produces a dictionary called "data" which is indexed by the keys "atmosphere", "cmb" and "noise". In each entry is an array where the first dimension corresponds to detector index, and the second dimension to the time sample index. 

### Caution

Gaussian process regression has cubic complexity, which scales poorly (especially when coded in Python). Simulating large swaths of atmosphere at high resolutions can be extremely slow, so don't go crazy with the input parameters. 

This package also produces large arrays: 1000 detectors sampling at 50 Hz for an hour is well over a gigabyte of data. 





Below: a theoretical array of 30,000 detectors (each with a resolution of 2 arcminutes) observes a simulated map of the cosmic microwave background through 16 simulated layers of atmospheric emission, at an elevation of 45 degrees and employing a constant-elevation scan of 2 degrees of azimuth per second. 

![Watch the video](https://user-images.githubusercontent.com/41275226/115489537-539c2400-a22a-11eb-9f3f-013b4c5e8f6a.mp4)

