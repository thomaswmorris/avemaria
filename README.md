# Auto-Regressive Atmospheric Modeling (ARAM)

ARAM is a python-based package that simulates turbulent atmospheric emission using a auto-regressive gaussian process framework, for applications in observational astronomy. Below: a distribution of turbulent water vapor moves through the field of view of the observer. 

![Watch the video](https://user-images.githubusercontent.com/41275226/117068746-acbf8980-acf9-11eb-8016-64fa01e12a77.mp4)

## Background and Methodology

Atmospheric modeling is an important step in both experiment design and subsequent data analysis for ground-based cosmological telescopes observing the cosmic microwave background (CMB). The next generation of CMB experiments will be marked by a huge increase in data acquisition.

By default, atmospheric emission simulated by ARAM is designed to resemble that observed by the Atacama Cosmology Telescope (ACT) and the Atacama B-Mode Search (ABS), two CMB experiments in the Atacama Desert. An outline of the method and products used is outlined [here](https://github.com/tomachito/aram/blob/main/README.md). 




## Examples and Usage 

To install ARAM with PyPi, put 

```console
pip install aram
```
in the command line. The main tool of the ARAM module is the model object, which can be easily intitialized as 

```python
from aram import aram

my_model = aram.model()
```

Different models can be initialized by configurating different aspects of the model. 

### Arrays

The array defines the set of detectors that observe the sky. The input parameters, defined as their default values along with units and a brief description:

```python
my_model = aram.model(array_config={'array_shape' : 'hex',  # (none)  - the shape of the array. 
                                    'n_detectors' : 600,    # (none)  - the number of detectors. 
                                    'array_fov'   : 2,      # (deg)   - maximum separation of detectors, in degrees. 
                                    'fwhm'        : 5 / 60, # (deg)   - the full-width at half-maximum of the beams.
                                    'band'        : 1.5e11, # (Hz)    - the observing band of the detectors, in Hz. 
                                    'pink'        : 0,      # (mK/Hz) - scale factor for the pink noise spectrum. 
                                    'white'       : 0})     # (mK/Hz) - scale factor for the white noise spectrum. 
```
Alternatively, the array can be configured manually by supplying an array of values for each parameter. In this case, the first three parameters are replaced by
```python
my_model = aram.model(array_config={'offset_x' : some_array_of_offsets_x,  # in degrees
                                    'offset_y' : some_array_of_offsets_y}) # in degrees
```



### Observations

The observation defines the 

### Sites

The site determines the pointing angles of 



Below: a theoretical array of 30,000 detectors (each with a resolution of 2 arcminutes) observes a simulated map of the cosmic microwave background through 16 simulated layers of atmospheric emission, at an elevation of 45 degrees and employing a constant-elevation scan of 2 degrees of azimuth per second. 

![Watch the video](https://user-images.githubusercontent.com/41275226/115489537-539c2400-a22a-11eb-9f3f-013b4c5e8f6a.mp4)

