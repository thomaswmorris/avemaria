
import ephem
import pickle
import scipy as sp
import numpy as np
import pandas as pd
    
from datetime import datetime
from datetime import timezone
from maria import tools

from importlib import resources

import weathergen

# an array object describes the array of detectors. all of the arguments are dimensionful arrays of the same length 



default_atmosphere_config = {'n_layers'                       : 10,         # how many layers to simulate, based on the integrated atmospheric model 
                                  'weather_generation_method' : 'mean',  # 
                                  'min_depth'                 : 500,      # the height of the first layer 
                                  'max_depth'                 : 5000,      # 
                                  'atmosphere_rms'            : 50,  
                                  'turbulence_model'          : 'scale_invariant',
                                  'outer_scale'               : 500}

default_site_config = {'site' : 'ACT',
                       'time' : datetime.now(timezone.utc).timestamp(),
             'weather_gen_method' : 'random',
                     'region' : 'atacama' }


default_array_config = {'shape' : 'hex',
                            'n' : 600,       # maximum number of detectors
                          'fov' : 10,
                         'band' : 1.5e11 }         # maximum span of array

default_beams_config = {'optical_type' : 'diff_lim',
                           'beam_type' : 'tophat',
                        'primary_size' : 5.5,
                            'beam_res' : .5 }     
     
default_pointing_config = {'scan_type' : 'CES',
                           'duration'  : 600,'samp_freq' : 20,
                         'center_azim' : 0, 'center_elev'  : 90, 
                            'az_throw' : 0, 'az_speed' : 1.5,
                            'el_throw' : 0, 'el_speed' : 1.5}

class atmosphere():
    
    def __init__(self, config=None):
        
        if config==None:
            print('No atm config specified, using default layers.')
            self.config = default_atmosphere_config.copy()
        else:
            self.config = config.copy()
        
        if self.config==None:
            print('No site config specified, using the ACT site.')
            self.config = default_site_config
        
        use_auto_depths   = np.all(np.isin(['min_depth','max_depth','n_layers'],list(self.config)))
        use_manual_depths = np.all(np.isin(['depths'],list(self.config)))
        
        if use_manual_depths: 
            
            if isinstance(self.depths, np.ndarray):
                self.config['min_depth'] = self.depths.min()
                self.config['max_depth'] = self.depths.max()
                self.config['n_layers']  = len(self.depths)
            else:
                raise Exception('\'depths\' parameter must be a numpy array.') 
                
        if not (use_auto_depths or use_manual_depths):
            for arg in ['min_depth','max_depth','n_layers']:
                self.config[arg] = default_atmo_config[arg]
            use_auto_depths = True
        if use_auto_depths:  
            self.depths = np.linspace(self.config['min_depth'], self.config['max_depth'], self.config['n_layers'])
            
            #raise Exception('Could not build atmospheric layers. Please specify the \'min_depth\', \'max_depth\', and \'n_layers\' parameters, or else enter an array of heights.')
        
        necessary_args = ['turbulence_model','outer_scale','atmosphere_rms']
        for arg in necessary_args:
            if not arg in list(self.config):
                self.config[arg] = default_atmosphere_config[arg]
                
        if self.config['turbulence_model'] == 'scale_invariant':
            
            self.matern = lambda r,r0,nu : 2**(1-nu)/sp.special.gamma(nu)*sp.special.kv(nu,r/r0+1e-10)*(r/r0+1e-10)**nu
            
            
class site():
    
    def __init__(self, config=None):
        
        if config==None:
            print('No site specified, using the ACT site.')
            self.config = default_site_config.copy()
        else:
            self.config = config.copy()
        
        if 'site' in list(self.config):
            
            with resources.path("maria", "site_info.csv") as f:
                self.site_df = pd.read_csv(f, index_col=0)
    
            site_list = '\n\nsite' + 5*' ' + 'region' + 7*' ' + 'weather' +  3*' ' + 'longitude' + 3*' ' + 'latitude' + 2*' ' + 'height'
            site_list += '\n' + (len(site_list)-2)*'#'
            
            for sitename in list(self.site_df.index):
                name,sup,loc,lon,lat,hgt = [self.site_df.loc[sitename,key] for key in ['longname','supported','region','longitude','latitude','altitude']]
                lon_name = f'{np.round(np.abs(lon),3):>8}°' + ['W','E'][int(lon>0)]
                lat_name = f'{np.round(np.abs(lat),3):>8}°' + ['S','N'][int(lat>0)]
                site_list += f'\n{sitename:<8} {loc:<12} {sup:<8} {lon_name} {lat_name} {hgt:>6.0f}m'
                
            if not self.config['site'] in self.site_df.index:
                raise Exception('\'' + self.config['site'] + '\' is not a supported site! Supported sites are:' + site_list)
                
            site_info = self.site_df.loc[self.config['site']]
            region    = site_info['region']
            latitude  = site_info['latitude']
            longitude = site_info['longitude']
            altitude  = site_info['altitude']
                
        else:      
            parameters = ['time','location','latitude','longitude','altitude']
            if not np.all(np.isin(parameters),list(self.config)):
                par_error = 'Please supply '
                for par in parameters:
                    par_error += f'\'{par}\''
                raise Exception(par_error)
            else:
                region    = self.config['region']
                latitude  = self.config['latitude']
                longitude = self.config['longitude']
                altitude  = self.config['altitude']

        self.observer = ephem.Observer()
        self.observer.lat, self.observer.lon, self.observer.elevation = str(latitude), str(longitude), altitude
        self.region    = region
        self.timestamp = self.config['time']
        self.observer.date = datetime.fromtimestamp(self.timestamp)
            
        if 'weather_gen_method' in list(self.config):
            
            self.weather = weathergen.generate(region=self.region,
                                               time=self.timestamp,
                                               method=self.config['weather_gen_method'])
            
            print(self.weather['water_density'][0])
            
class array():
    
    def __init__(self, config=None):
       

        if config==None:
            print('No array specified, using the ACT array.')
            self.config = default_array_config.copy()
        else:
            self.config = config.copy()
            
        #print(self.config)
            
        # which offset method has been entered? it should be exactly one of these
        
        
        if 'shape' in list(self.config):
            
            self.config['offset_x'], self.config['offset_y'] = tools.make_array(self.config['shape'], 
                                                                                self.config['fov'], 
                                                                                self.config['n']) 
            
        else: 
            pass
                
            
        self.z = np.radians(self.config['offset_x']) + 1j*np.radians(self.config['offset_y']); self.z -= self.z.mean()
        self.x = np.real(self.z)
        self.y = np.imag(self.z)
        self.n = len(self.z)
        
        if not 'band' in list(self.config):
            self.config['band'] = default_array_config['band']
        if isinstance(self.config['band'],np.ndarray):
            if len(self.config['band']) != self.n:
                raise Exception(f'ShapeError: shapes of detector bands does not match shape of offsets.')
            self.band = self.config['band']
        else:
            self.band = np.repeat(self.config['band'], self.n)   
                     
        self.wavelength = 2.998e8 / self.band
                     
             
class pointing():
    
    def __init__(self, config=None):
        
        if config==None:
            print('No pointing specified, defaulting to a 10-minute zenith stare at 20 Hz.')
            self.config = default_pointing_config
        else:
            self.config = config.copy()
            
        if 'scan_type' in list(self.config):
            
            self.duration = self.config['duration']
            self.dt   = 1 / self.config['samp_freq']
            self.time = np.arange(0, self.duration, self.dt)
            self.nt   = len(self.time)
            self.f_   = np.fft.fftfreq(self.nt,self.dt)
            
            self.center_azim, self.center_elev = np.radians(self.config['center_azim']), np.radians(self.config['center_elev'])
            
            if self.config['scan_type']=='CES':
                
                self.scan_freq  = self.config['az_speed'] / (4*self.config['az_throw']+1e-16)
                self.focal_azim = (self.center_azim + np.radians(self.config['az_throw'])*sp.signal.sawtooth(np.pi/2 + 2*np.pi*self.scan_freq*self.time,width=.5)) % (2*np.pi)
                self.focal_elev = self.center_elev + np.zeros(self.nt)
                
            if self.config['scan_type']=='lissajous_box':

                focal_x = np.radians(self.config['x_throw']) * np.sin(2*np.pi*self.time/self.config['x_period'])
                focal_y = np.radians(self.config['y_throw']) * np.sin(2*np.pi*self.time/self.config['y_period'])
                
                self.focal_azim, self.focal_elev = tools.from_xy(focal_x,focal_y,self.center_azim,self.center_elev)
                
            if self.config['scan_type']=='lissajous_daisy':

                focal_r = np.radians(self.config['throw']) * np.sin(2*np.pi*self.time/self.config['r_period'])
                focal_p = 2*np.pi*self.time/self.config['p_period']
    
                focal_x, focal_y = focal_r * np.cos(focal_p), focal_r * np.sin(focal_p)
                
                self.focal_azim, self.focal_elev = tools.from_xy(focal_x,focal_y,self.center_azim,self.center_elev)
            
        else:
            
            self.focal_azim = self.config['focal_azim']
            self.focal_elev = self.config['focal_elev']
            self.time       = self.config['time']
            self.duration   = self.time.max() - self.time.min()
            self.dt = np.gradient(self.time).mean()
            self.nt = len(self.time)
            self.f_ = np.fft.fftfreq(self.nt,self.dt)
        
        
        #for arg in list(default_pointing_config):
        #    if not arg in list(self.config):
        #        self.config[arg] = default_pointing_config[arg]
        
        #self.obs.mlst0 = float(self.site.site.sidereal_time()) 
        #self.obs.mlst_ = self.obs.mlst0 % (2*np.pi) + ((2*np.pi / 86400) * self.obs.t_) % (2*np.pi) # mean local sidereal time 


class beams():
    
     def __init__(self, config=None):
        
        if config==None:
            print('No beams specified, defaulting to ACT beams.')
            self.config = default_beams_config
        else:
            self.config = config.copy()

        for arg in list(default_beams_config):
            if not arg in list(self.config):
                self.config[arg] = default_beams_config[arg]
                
                
        if self.config['optical_type'] == 'diff_lim':
            
            self.aperture = self.config['primary_size']
            self.beam_res = self.config['beam_res']
            self.waist = lambda z, w_0, f : w_0 * np.sqrt(1 + np.square(2.998e8 * z) / np.square(f * np.pi * np.square(w_0)))

    
