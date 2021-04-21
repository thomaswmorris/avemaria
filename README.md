# ARAM (Auto-Regressive Atmospheric Modeling)

Uneven emission of the lower atmosphere is one of the largest sources of interference for ground-based telescopes observing in millimeter and sub-millimeter regimes. Accurately and efficiently simulating this interference is crucial for advances in data-analysis techniques. 

ARAM is a Python package that simulates and realizes time-ordered atmospheric fluctuations in ground-based telescopes. This package uses a sparse auto-regressive precision model of the atmosphere in order to simulate the atmosphere at high resolution. As opposed to instantaneous methods of generation, which are limited to several tens of thousands atmosphere elements, ARAM can procedurally simulate several tens of millions of atmospheric elements and can thus achieve a spatial resolution of an arcminute and a temporal resolution of a few milliseconds. The key feature of ARAM is its speed: it can simulate days of comprehensively modeled time-ordered data in a few minutes. 

The atmospheric model is comprised of turbulent fluctuations at different atmospheric layers, which have variable physical parameters. The parameters describe different aspects of the atmosphere, and include wind speed, wind bearing, average water vapor mass density, and temparature. ARAM can automatically generate parameters using a weather model based on historical MERRA-2 reanalysis for supported sites (the Atacama Desert, Tibet, and the South Pole). 

The package has addition modules for simulating non-atmospheric signals. It can generate and return a time-ordered observation of a realization of the cosmic microwave background (CMB) for a given spectrum, but defaults to the T-T, E-T, and E-E spectra derived from ACTPol and WMAP datasets. It also has tools for reproducing noise characteristics native to polarization-sensitive bolometers. 

Below is an example of the package in action: atmospheric emission, cosmic background radiation, and noise were simulated at 100 Hz for a 180 by 180 (32,400) pixel array spanning a 3-degree-wide square on the sky, scanning at a speed of 2 degrees per second. Each pixel has a full-width at half-maximum of 1 arcminute, and atmosphere and CMB were simulated at a resolution of 0.25 arcminutes.  

![Watch the video](https://user-images.githubusercontent.com/41275226/115489537-539c2400-a22a-11eb-9f3f-013b4c5e8f6a.mp4)

