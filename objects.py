import sys
import time
import ephem
import pickle
import scipy as sp
import numpy as np
import numpy.linalg as la
from scipy import signal, stats
from datetime import datetime
from datetime import timezone
#from IPython.display import display,clear_output

from aram import resources
import gc
import os
from . import tools

import pkg_resources

# an array object describes the array of detectors. all of the arguments are dimensionful arrays of the same length 


#def set_default_arg(arg):
    

class array():
    def __init__(self, array_config=None, fov_tol=2e-2, verbose=True):
        default_array_config = {'array_shape' : 'hex',
                                          'n' : 600,      # maximum number of detectors
                                  'array_fov' : 5,        # maximum span of array
                                     'optics' : 'gaussian',
                                       'fwhm' : 5 / 60,   # in degrees
                                       'band' : 1.5e11,   # observing band (in Hz)
                                       'pink' : 4e-1,     # pink noise scaling 
                                      'white' : 1e0}      # white noise scaling 
        
        
        '''
        The array object describes a set of detectors and an optical system that houses them. There are several different
        ways to describe an array, and so we need to be flexible. Here are the independent groups of parameters:
        
        1 - Specify the arrangement of detectors. This can be entered as:
        
            a - 'offset_x' and 'offset_y' parameters, in degrees
            b - 'array_shape', 'n_detectors', and 'array_fov' parameters, for supported array shapes.
            
        2 - Specifiy the optics. There are quite a few ways to do this:
        
            a - 
            
            
        3 - Specify the noise levels. This can be a constant float or an array of floats with the same shape as the detectors. 
        
        '''
        
        # ================ OFFSETS ================
        
        # use the default config if none supplied
        if array_config==None:
            print('No array config specified, using default array.')
            array_config = default_array_config
            
        # which offset method has been entered? it should be exactly one!
        do_manual_offsets = np.all(np.isin(['offset_x','offset_y'],list(array_config)))
        do_auto_offsets   = np.all(np.isin(['array_shape','array_fov','n'],list(array_config)))
        if do_manual_offsets:
            if len(array_config['offset_x']) != len(array_config['offset_y']):
                raise Exception(f'ShapeError: shapes of x and y offsets do not match!')
        if do_manual_offsets & do_auto_offsets:# or (~manual_offsets or ~auto_offsets):
            raise Exception(f'Bad format! Please enter exactly one set of array specifications.')
        if not (do_manual_offsets or do_auto_offsets):
            print('No offsets specified, using default offsets.')
            array_config['array_shape'], array_config['array_fov'], array_config['n'] = default_array_config['array_shape'], default_array_config['array_fov'], default_array_config['n']
        if not do_manual_offsets:
            array_config['offset_x'], array_config['offset_y'] = tools.make_array(array_config['array_shape'], 
                                                                                  array_config['array_fov'], 
                                                                                  array_config['n']) 
            
        self.z = np.radians(array_config['offset_x']) + 1j*np.radians(array_config['offset_y']); self.z -= self.z.mean()
        self.x = np.real(self.z)
        self.y = np.imag(self.z)
        self.n = len(self.z)
        
        if not 'band' in list(array_config):
            array_config['band'] = default_array_config['band']
        if isinstance(array_config['band'],np.ndarray):
            if len(array_config['band']) != self.n:
                raise Exception(f'ShapeError: shapes of detector bands does not match shape of offsets.')
            self.band = array_config['band']
        else:
            self.band = np.repeat(array_config['band'], self.n)   
                     
        self.wavelength = 2.998e8 / self.band
                     
        # ================ OPTICS ================
        
        if not 'optics' in list(array_config):
            array_config['optics'] = default_array_config['optics']
                     
        self.beam_type = array_config['optics']
        
        if self.beam_type == 'diff_lim':
                     
            if not 'aperture' in list(array_config):
                array_config['aperture'] = default_array_config['aperture']
            self.aperture = array_config['aperture']
            self.fwhm = 1.22 * self.wavelength / self.aperture
                     
        if self.beam_type == 'gaussian':
            
            if not 'fwhm' in list(array_config):
                array_config['fwhm'] = default_array_config['fwhm']
            if isinstance(array_config['fwhm'],np.ndarray):
                if len(array_config['fwhm']) != self.n:
                    raise Exception(f'ShapeError: shapes of resolutions does not match shape of offsets.')
                self.fwhm = np.radians(array_config['fwhm'])
            else:
                self.fwhm = np.radians(np.repeat(array_config['fwhm'], self.n))
                     
        for arg in ['pink','white']:
            if not arg in list(array_config):
                array_config[arg] = default_array_config[arg]
            if isinstance(array_config[arg],np.ndarray):
                if len(array_config[arg]) != self.n:
                    raise Exception(f'ShapeError: shapes of resolutions does not match shape of offsets.')
                val = array_config[arg]; exec(f'self.{arg} = val')
            else:
                val = array_config[arg]; exec(f'self.{arg} = np.repeat(val, self.n) ')
                     
        self.fov_r = (1 + fov_tol) * (np.abs(self.z).max()) + self.fwhm.max() 
        
            
#arr = array(array_config)


class site():
    def __init__(self, site_config=None):
        
        with open(pkg_resources.resource_filename('aram', 'resources/sites/site_dict'), 'rb') as f:
            self.site_dict = pickle.load(f)

        site_list = '\nsite' + 5*' ' + 'location' + 4*' ' + 'longitude' + 2*' ' + 'latitude' + 2*' ' + 'height'
        site_list += '\n' + (len(site_list)-1)*'-'
        for sitename in list(self.site_dict):
            name,loc,lon,lat,hgt = [self.site_dict[sitename][key] for key in ['longname','location','longitude','latitude','altitude']]
            lon_name = f'{np.round(np.abs(lon),3):>8}°' + ['W','E'][int(lon>0)]
            lat_name = f'{np.round(np.abs(lat),3):>7}°' + ['S','N'][int(lon>0)]
            site_list += f'\n{sitename:<8} {loc:<10} {lon_name} {lat_name} {hgt:>6.0f}m'

        self.site = ephem.Observer()
        default_site_config = {'site' : 'ACT',
                               'time' : datetime.now(timezone.utc).timestamp()}
                
        if site_config==None:
            print('No site config specified, using the ACT site.')
            site_config = default_site_config
            
        manual_site = np.all(np.isin(['location','latitude','longitude','altitude'],list(site_config)))
        
        for arg in ['site','time']:
            if not arg in list(site_config):
                site_config[arg] = default_site_config[arg]
        if manual_site: 
            location, latitude, longitude, altitude = site_info['location'], site_config['latitude'], site_config['longitude'], site_config['altitude']
        else:
            if not site_config['site'] in list(self.site_dict):
                raise Exception(f'Not a supported site. Available sites are:\n' + site_list)
            else:
                site_info = self.site_dict[site_config['site']]
                self.is_auto = True
                location, latitude, longitude, altitude = site_info['location'], site_info['latitude'], site_info['longitude'], site_info['altitude']

        #print(latitude, longitude, altitude)
        self.location = location
        self.site.lat, self.site.lon, self.site.elevation = str(latitude), str(longitude), altitude
        self.timestamp = site_config['time']
        
        self.site.date = datetime.fromtimestamp(self.timestamp)
        
        
        
class observation():
    def __init__(self, site_config=None, obs_config=None):
        
        default_obs_config = {'duration'  : 600,'samp_freq' : 20,
                              'center_az' : 0, 'center_el'  : 90, 'az_throw' : 0, 'az_speed' : 0}

        if obs_config==None:
            print('No obs config specified, defaulting to a 10-minute zenith stare at 20 Hz.')
            obs_config = default_obs_config
            
        for arg in list(default_obs_config):
            if not arg in list(obs_config):
                obs_config[arg] = default_obs_config[arg]
            
        self.duration = obs_config['duration']
        self.dt    = 1/obs_config['samp_freq']
        self.t_    = np.arange(0,self.duration,self.dt); self.T = len(self.t_); self.f_ = np.fft.fftfreq(self.T,self.dt)
        self.c_az, self.c_el = np.radians(obs_config['center_az']), np.radians(obs_config['center_el'])
        self.r_az, self.v_az = np.radians(obs_config['az_throw']),  np.radians(obs_config['az_speed'])
        self.c_az_ = (self.c_az + self.r_az*sp.signal.sawtooth(np.pi/2 + 2*np.pi*self.v_az*(self.t_-self.t_.min())/(4*self.r_az+1e-16),width=.5)) % (2*np.pi)
        self.c_az_[self.c_az_>np.pi] -= 2*np.pi
        self.c_el_ = self.c_el + np.zeros(self.T)


