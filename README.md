# Auto-Regressive Atmospheric Modeling

ARAM is a python-based package that simulates turbulent atmospheric emission using a auto-regressive gaussian process framework, for applications in observational astronomy. 

Below: a theoretical array of 30,000 detectors (each with a resolution of 2 arcminutes) observes a simulated map of the cosmic microwave background through 16 simulated layers of atmospheric emission, at an elevation of 45 degrees and employing a constant-elevation scan of 2 degrees of azimuth per second. 

![Watch the video](https://user-images.githubusercontent.com/41275226/115489537-539c2400-a22a-11eb-9f3f-013b4c5e8f6a.mp4)

## Background and Methodology

Atmospheric modeling is an important step in both experiment design and subsequent data analysis for ground-based cosmological telescopes observing the cosmic microwave background (CMB). The next generation of CMB experiments will be marked by a huge increase in data acquisition.

By default, atmospheric emission simulated by ARAM is designed to resemble that observed by the Atacama Cosmology Telescope (ACT) and the Atacama B-Mode Search (ABS), two CMB experiments in the Atacama Desert. An outline of the method and products used is outlined [here](https://github.com/tomachito/aram/blob/main/README.md). 


## Examples and Usage 

The main tool of the ARAM module is the model object, which can be easily intitialized as 

```python
from aram import aram

my_model = aram.model()
```

Different models can be initialized by configurating different aspects of the model. 

### Arrays

The array defines the set of detectors that observe the sky. Detector positions can be supplied to the manually:

```python
my_model = aram.model(array_config={'offset_x' : some_array_of_offsets_x,  
                                    'offset_y' : some_array_of_offsets_y})
```
or automatically by the package 
```python
my_model = aram.model(array_config={'array_shape' : 'hex',
                                    'n_detectors' : 600,
                                    'array_fov'   : 2})
```
The other detectors characteristics are the beam resolution, observing band, and noise characteristics, which are supplied as 
```python
my_model = aram.model(array_config={'fwhm'  : 5 / 60, # the far-field full-width at half-maximum of the beams, in degrees. Defaults to 5 arcminutes. 
                                    'band'  : 1.5e11, # the observing band of the detectors, in Hz. Defaults to 150 GHz. 
                                    'pink'  : 0,     # scale factor for the pink noise spectrum, S(f) = a/f. Defaults to zero. 
                                    'white' : 0})     # scale factor for the white noise spectrum, S(f) = a. Defaults to zero. 
```


### Observations

The observation defines the 

### Sites

The site determines the pointing angles of 


