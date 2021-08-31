# Modeling Auto-Regressive Integrated Atmosphere (maria)

maria is a python-based package that simulates turbulent atmospheric emission using a auto-regressive gaussian process framework, for applications in observational astronomy. Below: a distribution of turbulent water vapor moves through the field of view of the observer. 

![Watch the video](https://user-images.githubusercontent.com/41275226/117068746-acbf8980-acf9-11eb-8016-64fa01e12a77.mp4)

## Background

Atmospheric modeling is an important step in both experiment design and subsequent data analysis for ground-based cosmological telescopes observing the cosmic microwave background (CMB). The next generation of ground-based CMB experiments will be marked by a huge increase in data acquisition: telescopes like [AtLAST](https://www.atlast.uio.no) and [CMB-S4](https://cmb-s4.org) will consist of hundreds of thousands of superconducting polarization-sensitive bolometers sampling the sky. This necessitates new methods of efficiently modeling and simulating atmospheric emission at small angular resolutions, with algorithms than can keep up with the high throughput of modern telescopes.

maria simulates layers of turbulent atmospheric emission according to a statistical model derived from observations of the atmosphere in the Atacama Desert, from the [Atacama Cosmology Telescope (ACT)](https://lambda.gsfc.nasa.gov/product/act/) and the [Atacama B-Mode Search (ABS)](https://lambda.gsfc.nasa.gov/product/abs/). It uses a sparse-precision auto-regressive Gaussian process algorithm that allows for both fast simulation of high-resolution atmosphere, as well as the ability to simulate arbitrarily long periods of atmospheric evolution. 

## Methodology

maria auto-regressively simulates an multi-layeed two-dimensional "integrated" atmospheric model that is much cheaper to compute than a three-dimensional model, which can effectively describe time-evolving atmospheric emission. 

## Examples and Usage 

To install MARIA with PyPi, run

```console
pip install maria
```
The main tool of the maria module is the model object. The default model can be easily intitialized as 

```python
import maria

default_model = maria.model()
```

Different models can be initialized by configurating different aspects of the model.

### Arrays

The array config defines the set of detectors that observe the sky, along with the properties of their optics and noise. The array is specified as a dictionary

```python
array_config = {'array_shape' : 'hex',      # shape of array
                  'array_fov' : .8,         # maximum span of array (deg)
                          'n' : 500,        # maximum number of detectors (deg)
                       'band' : 150e9,      # observing band (in Hz)
                       
array_config = {'shape' : 'hex',   # The shape of the distribution of detectors. Supported shapes are `hex', 'square', and 'flower'. 
                    'n' : 10000,   # The number of detectors in the array.  
                  'fov' : 2,       # The maximum width of the array's field-of-view on the sky, in degrees. 
                 'band' : 1.5e11}  # The observing band of the detector, in Hz. 

```
Alternatively, the array can be configured manually by supplying an array of values for each parameter. In this case, the first three parameters are replaced by

```python
array_config={'offset_x' : some_array_of_offsets_x,  # in degrees
              'offset_y' : some_array_of_offsets_y}) # in degrees
```

### Observations

The pointing config defines the time-ordered parameters of the simulation. Below is the config for a constant-elevation scan (CES) that observes at 90+/-45 degrees of azimuth and 60 degrees of elevation, sampling at 100Hz for ten minutes. 

```python
pointing_config = {'scan_type' : 'CES', # scan pattern
                    'duration' : 600,   # duration of the observation, in seconds 
                   'samp_freq' : 100,   # sampling rate, in Hz
                 'center_azim' : 90,    # azimuth of the center of the scan, in degrees
                     az_throw' : 45,    # half of the azimuth width of the scan, in degrees
                 'center_elev' : 60,    # observing elevation of the scan, in degrees
                    'az_speed' : 1.5}   # scanning speed of the array, in degrees per second
```
Alternatively, the pointing data may be given manually

```python
pointing_config = {'time' : some_array_of_timestamps, # in seconds
             'focal_azim' : some_array_of_azimuths,   # in degrees
             'focal_elev' : some_array_of_elevations} # in degrees
```
where focal_azim and focal_elev describe the angular pointing of the center of the array. 

### Sites

The site determines the motion of celestial sources as the earth rotates under the telescope, as well as the  

```python
site_config = {'site' : 'ACT',
               'time' : datetime.now(timezone.utc).timestamp(),
 'weather_gen_method' : 'random'} 
```

For supported sites maria will generate weather data, which are used to inform the atmospheric simulation. Weather data are quantitatively realistic for a given site, altitude, time of day, and time of year, and are generated using the [weathergen](https://github.com/tomachito/weathergen) package. Longitude, latitude, and altitude may also be entered manually:

```python
site_config = {'site' : 'ACT',
               'time' : datetime.now(timezone.utc).timestamp(),
 'weather_gen_method' : 'random',
             'region' : 'atacama' } 
```



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
from maria import maria

my_model = maria.model(array_config=array_config,
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

