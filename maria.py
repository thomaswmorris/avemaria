import sys
import time
import ephem
import pickle
import scipy as sp
import numpy as np
import healpy as hp
import numpy.linalg as la
from scipy import signal, stats
from datetime import datetime
from datetime import timezone
from tqdm import tqdm
import gc
import os

from . import resources
from . import tools
from . import objects

import pkg_resources
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#exec(open('mpl_defs').read())

class model():
    
    def __init__(self,array_config=None, site_config=None, obs_config=None, model_config=None, other_config=None, full_compute=False, verbose=False):
        
        default_model_config = {'n_layers'     : 4,
                                'weather_type' : 'random',
                                'min_height'   : 1000,
                                'max_height'   : 3000,
                                'res'          : 1,
                                'half_height'  : 1500,
                                'atm_rms'      : 50,  
                                'wind_speed'   : 50,
                                'wind_bearing' : 0,
                                'cov_type'     : 'matern',
                                'outer_scale'  : 200,
                                'power_law'    : -8/3}  
        
        default_other_config = {'fov_tol' : 5e-2, 'n_dz_samp' : 16, 'n_ca_samp' : 16, 'n_st_samp' : 16, 'n_sam_max' : 10000}
        
        # MODEL CONFIGURATION
        if model_config==None:
            self.model_is_default = True  
            print('No model config specified, using default atmospheric model.')
            model_config = default_model_config
            
        if other_config==None:
            self.other_is_default = True  
            other_config = default_other_config
            
        self.array = objects.array(array_config=array_config)
        self.site  = objects.site(site_config=site_config)
        self.obs   = objects.observation(obs_config=obs_config)
        
        self.default_interp_method = 'rectilinear' if np.sum(self.array.fwhm[0]==self.array.fwhm)==self.array.n else 'brute_gaussian'

        self.obs.mlst0 = float(self.site.site.sidereal_time()) 
        self.obs.mlst_ = self.obs.mlst0 % (2*np.pi) + ((2*np.pi / 86400) * self.obs.t_) % (2*np.pi) # mean local sidereal time        
        self.ca_samp = self.obs.c_az + np.linspace(-self.obs.r_az,self.obs.r_az,other_config['n_ca_samp'])#+ np.pi/2
        self.dz_samp = self.array.fov_r * np.exp(1j*np.linspace(0,2*np.pi,other_config['n_dz_samp']+1)[:-1])
        self.st_samp = np.linspace(self.obs.mlst_[0],self.obs.mlst_[-1],other_config['n_st_samp'])
        self.dx_samp, self.dy_samp = np.real(self.dz_samp), np.imag(self.dz_samp)
        self.fov_simp = sp.spatial.Delaunay(np.vstack([self.dx_samp, self.dy_samp]).T, incremental=False, qhull_options=None)

        DZ_SAMP, ST_SAMP, CA_SAMP = np.meshgrid(self.dz_samp, self.st_samp, self.ca_samp)
        DX_SAMP, DY_SAMP = np.real(DZ_SAMP), np.imag(DZ_SAMP)
        AZ_SAMP, EL_SAMP = tools.from_xy(DX_SAMP, DY_SAMP, CA_SAMP, self.obs.c_el)
        RA_SAMP, DE_SAMP = tools.ae_to_rd(AZ_SAMP, EL_SAMP, ST_SAMP, self.site.site.lat)

        self.sp_dx,self.sp_dy  = DX_SAMP[0], DY_SAMP[0]
        self.sp_a,self.sp_e    = AZ_SAMP[0], EL_SAMP[0]
        self.rel_z             = np.exp(1j*(np.pi/2-AZ_SAMP[0])) / np.tan(EL_SAMP[0])
        self.rel_x,self.rel_y  = np.real(self.rel_z),np.imag(self.rel_z)
        self.sp_ra,self.sp_dec = RA_SAMP, DE_SAMP
        
        #self.az_, self.el_  = aram_tools.xy_to_ae(self.array.x[:,None],self.array.y[:,None],self.c_az_[None,:],self.c_el_[None,:])
        #self.ra_, self.dec_ = aram_tools.ae_to_rd(self.az_, self.el_,self.mlst_[None,:],self.site.lat)
        
        auto_heights   = np.all(np.isin(['min_height','max_height','n_layers'],list(model_config)))
        manual_heights = np.all(np.isin(['heights'],list(model_config)))
        if auto_heights:  
            self.heights = np.linspace(model_config['min_height'], model_config['max_height'], model_config['n_layers'])
        if manual_heights:  
            if isinstance(model_config['heights'], np.ndarray):
                self.heights = model_config['heights']
            else:
                raise Exception('\'heights\' parameter must be a numpy array.')
        if not (auto_heights or manual_heights):
            raise Exception('Could not build atmospheric layers. Please specify the \'min_height\', \'max_height\', and \'n_layers\' parameters, or else enter an array of heights.')
        
        for arg in list(default_model_config):
            if not arg in list(model_config):
                model_config[arg] = default_model_config[arg]
            thing = model_config[arg]
            exec(f'self.{arg} = thing')
                
        self.n_layers = len(self.heights)
        
        if self.site.is_auto:
            
            self.week = np.minimum(51,int((datetime.fromtimestamp(self.site.timestamp).timetuple().tm_yday - 1)/7))
            self.hour = datetime.fromtimestamp(self.site.timestamp).timetuple().tm_hour - 1
            site, week, hour = self.site.location, self.week, self.hour

            with open(pkg_resources.resource_filename('aram', 'resources/weather/gen'), 'rb') as handle:
                norm_gen  = pickle.load(handle)
            with open(pkg_resources.resource_filename('aram', 'resources/weather/avg'), 'rb') as handle:
                avg_merra = pickle.load(handle)
            with open(pkg_resources.resource_filename('aram', 'resources/weather/rms'), 'rb') as handle:
                rms_merra = pickle.load(handle)
            with open(pkg_resources.resource_filename('aram', 'resources/weather/min'), 'rb') as handle:
                min_merra = pickle.load(handle)
            with open(pkg_resources.resource_filename('aram', 'resources/weather/max'), 'rb') as handle:
                max_merra = pickle.load(handle)
            with open(pkg_resources.resource_filename('aram', 'resources/weather/hgt'), 'rb') as handle:
                hgt_merra = pickle.load(handle)
            n_gen = norm_gen[site][week,hour].shape[1]; gw = np.ones(n_gen)
                
            # generate weather, within the 40-year range 
            if model_config['weather_type'] == 'random':
                while np.any(gw > max_merra[site][week][hour]) or np.any(gw < min_merra[site][week][hour]):    
                    gw = np.matmul(norm_gen[site][week][hour], np.random.standard_normal(norm_gen[site][week][hour].shape[0]))
                    gw *= rms_merra[site][week][hour]
                    gw += avg_merra[site][week][hour]
                gw_spl = gw.reshape((4,int(n_gen/4)))
                
            if model_config['weather_type'] == 'mean':
                print('using mean')
                gw_spl = avg_merra[site][week][hour].reshape((4,int(n_gen/4)))
            
            # transfer to a generated weather dictionary 
            self.gw_dict = {}
            for icol,col in enumerate(['LWVMD','T','U','V']):
                self.gw_dict[col] = gw_spl[icol]
            self.gw_dict['WVMD']  = np.exp(self.gw_dict['LWVMD'])

            # resample to model convention
            self.weather = {}
            for col in ['WVMD','T','U','V']:
                self.weather[col] = np.interp(self.heights,hgt_merra,self.gw_dict[col]) 
            z = self.weather['V'] + 1j*self.weather['U']
            self.weather['ws'], self.weather['wa'] = np.abs(z), np.angle(z)
            self.wind_bearing = self.weather['wa']
            self.wind_speed   = self.weather['ws']
            self.temperature  = self.weather['T']
            self.wvmd         = self.weather['WVMD']

        else:
                
            if isinstance(model_config['wind_speed'], np.ndarray):
                self.wind_speed = model_config['wind_speed']
            else:
                self.wind_speed = model_config['wind_speed'] * np.ones(self.n_layers)

            if isinstance(model_config['wind_bearing'], np.ndarray):
                self.wind_bearing = np.radians(model_config['wind_bearing'])
            else:
                self.wind_bearing = np.radians(model_config['wind_bearing']) * np.ones(self.n_layers)

            if not (len(self.wind_speed) == self.n_layers and len(self.wind_bearing) == self.n_layers):
                raise Exception('Shape mismatch in the heights and the wind velocities.')

        if model_config['cov_type'] == 'matern':
            self.r0, self.nu = model_config['outer_scale'], (-1-model_config['power_law'])/2
            self.C = lambda r: 2**(1-self.nu)/sp.special.gamma(self.nu)*sp.special.kv(self.nu,r/self.r0+1e-10)*(r/self.r0+1e-10)**self.nu
        
        self.array.wavelength = 2.998e8 / self.array.band
                
        self.depths   = self.heights / np.sin(self.obs.c_el)
        
        if self.array.beam_type == 'diff_lim':
            self.ang_sig_ = np.maximum(self.array.aperture, self.array.fwhm[:,None]*self.depths) / self.depths / (2*np.sqrt(2*np.log(2)))
            
        if self.array.beam_type == 'gaussian':
            self.ff_sig   = self.array.fwhm / (2*np.sqrt(2*np.log(2)))
            self.w0       = self.array.wavelength / (self.ff_sig * np.pi * 1.003)
            self.ang_sig_   = tools.gaussian_beam(self.depths[None,:],self.w0[:,None],l=self.array.wavelength[:,None],n=1.003)
            
        self.phy_sig_ = self.ang_sig_.min(axis=0) * self.depths
        self.atm_res_ = self.phy_sig_ * model_config['res'] * (2*np.sqrt(2*np.log(2)))
        self.d_orth_  = self.atm_res_ 
        self.d_para_  = self.atm_res_ 
        self.sim_dt_  = self.atm_res_ / self.wind_speed
        
        # time and angle fields for the simulation
        self.sim_t_    = [self.obs.t_[0] + np.arange(0,self.obs.duration+dt, dt) for dt in self.sim_dt_]
        self.sim_T_    = [len(sim_t_) for sim_t_ in self.sim_t_]; self.sim_f_ = [np.fft.fftfreq(T,dt) for T,dt in zip(self.sim_T_,self.sim_dt_)]
        self.tot_sim_T = np.sum(self.sim_T_)
        self.sim_c_az_ = [np.interp(t_,self.obs.t_,self.obs.c_az_) for t_ in self.sim_t_]
        self.sim_c_el_ = [np.interp(t_,self.obs.t_,self.obs.c_el_) for t_ in self.sim_t_]
        
        self.az_0, self.el_0 = tools.from_xy(self.array.x,self.array.y,0,self.obs.c_el)
        self.fov_az_, self.fov_el_ = tools.from_xy(self.dx_samp[:,None],self.dy_samp[:,None],self.sim_c_az_[0][None,:],self.sim_c_el_[0][None,:])    
        self.sim_min_az_, self.sim_max_az_ = self.fov_az_.min(axis=0), self.fov_az_.max(axis=0)
        self.sim_min_el_, self.sim_max_el_ = self.fov_el_.min(axis=0), self.fov_el_.max(axis=0)
        
        # CMB time-ordered stuff
        self.cmb_res       = self.array.fwhm.min() / 2
        self.sim_cmb_dt    = np.maximum(self.cmb_res / (self.obs.v_az + 1e-16), self.obs.dt)
        self.sim_cmb_t_    = self.obs.t_[0] + np.arange(-self.sim_cmb_dt,self.obs.duration+self.sim_cmb_dt,self.sim_cmb_dt)
        self.sim_cmb_c_az_ = np.interp(self.sim_cmb_t_,self.obs.t_,self.obs.c_az_)
        self.sim_cmb_c_el_ = np.interp(self.sim_cmb_t_,self.obs.t_,self.obs.c_el_)
        self.sim_cmb_mlst_ = np.interp(self.sim_cmb_t_,self.obs.t_,self.obs.mlst_)
       
        # CMB spatial stuff
        self.ra_0, self.dec_0 = self.sp_ra.mean(), self.sp_dec.mean()
        self.sp_rdx, self.sp_rdy = tools.to_xy(self.sp_ra, self.sp_dec, self.ra_0, self.dec_0)

        x_bins = np.arange(self.sp_rdx.min(), self.sp_rdx.max() + self.cmb_res, self.cmb_res)
        y_bins = np.arange(self.sp_rdy.min(), self.sp_rdy.max() + self.cmb_res, self.cmb_res)
        self.rdx = x_bins[1:]/2 + x_bins[:-1]/2
        self.rdy = y_bins[1:]/2 + y_bins[:-1]/2

        self.RDX, self.RDY = np.meshgrid(self.rdx, self.rdy)
        self.RA,  self.DEC = tools.from_xy(self.RDX, self.RDY, self.ra_0, self.dec_0)
        
        del self.fov_az_, self.fov_el_
        gc.collect()
        self.lay_z_, self.lay_h_ = [], []
        self.gen_z_, self.gen_h_ = [], []
        self.n_para_, self.n_orth_ = [], []
        self.lay_para_, self.lay_orth_ = [], []
        self.sam_i_ = []
        
        # this is messy, but necessary
        if verbose:
            self.DATE_STRING, self.TIME_STRING = datetime.fromtimestamp(self.site.timestamp).isoformat().split('T')
            print(f'observing time : {self.DATE_STRING} {self.TIME_STRING[:8]} UTC\n')
            print('  # | height (m) | res. (m) | shape (p, o) | n_sam | W (g/m3) | T (C) | V (m/s) | A (deg) | dt (s) |')
            print(100*'-')

        self.good_layer_ = np.ones(self.n_layers).astype(bool)
        
        tot_n_sam, tot_n_lay = 0,0
        for i_h,h in enumerate(self.heights):
            adj_bear = 3*np.pi/2 + self.wind_bearing[i_h]
            self.rot_ch_z    = h * self.rel_z * np.exp(1j*adj_bear)
            ch_para, ch_orth = np.real(self.rot_ch_z),np.imag(self.rot_ch_z)
            ch_para_max = np.maximum(ch_para.max(),ch_para.min()+self.r0)
            self.lay_para_.append(np.arange(ch_para.min(),ch_para_max+self.d_para_[i_h], self.d_para_[i_h]))
            self.lay_orth_.append(np.arange(ch_orth.min(),ch_orth.max()+self.d_orth_[i_h], self.d_orth_[i_h]))
            ORTH_G, PARA_G = np.meshgrid(self.lay_orth_[-1],self.lay_para_[-1])
            n_para, n_orth = len(self.lay_para_[-1]), len(self.lay_orth_[-1])
            self.n_para_.append(n_para)
            self.n_orth_.append(n_orth)
            z_ = PARA_G + 1j*ORTH_G
            self.lay_z_.append(z_*np.exp(-1j*adj_bear))
            self.lay_h_.append(h*np.ones(z_.shape))
            self.gen_z_.append((z_[0]-self.d_para_[i_h])*np.exp(-1j*adj_bear))
            self.gen_h_.append(h*np.ones(z_[0].shape))
            i_para_, i_orth_ = [],[]
            for i_para in np.r_[0,2**np.arange(np.ceil(np.log(n_para)/np.log(2))),self.n_para_[-1]-1]:
                i_orth_.append(np.unique(np.linspace(0,n_orth-1,int(np.maximum(n_orth/(i_para+1),16))).astype(int)))
                i_para_.append(np.repeat(i_para,len(i_orth_[-1])).astype(int))
            i_max  = np.maximum(0,np.max(i_para_[-1])) 
            self.sam_i_.append((np.concatenate(i_para_),np.concatenate(i_orth_)))
            n_    = len(z_.ravel())
            n_sam = len(self.sam_i_[-1][0].ravel())
            n_gen = len(self.gen_z_[-1].ravel())
            tot_n_lay += n_
            tot_n_sam += n_sam
            
            if n_sam > other_config['n_sam_max']:
                self.good_layer_[i_h] = False
            
            if verbose:
                if self.good_layer_[i_h]:
                    print(f' {i_h:>2} | {h:>10.02f} | {self.atm_res_[i_h]:>8.02f} | {repr(z_.shape):12} | {n_sam:>5.0f} | {1e3*self.wvmd[i_h]:>8.3f} | {self.temperature[i_h]:>5.1f}'
                          + f' | {self.wind_speed[i_h]:>7.1f} | {np.degrees(self.wind_bearing[i_h]):>7.1f} | {self.sim_dt_[i_h]:>6.03f}')
                else:
                    print(f' {i_h:>2} | {h:>10.02f} | WARNING: LAYER EXCLUDED')
                    
        if verbose:
            print(100*'-')
            #print('  # | height (m) | res. (m) | shape (p, o) | n_sam | W (g/m3) | T (C) | V (m/s) | A (deg) | dt (s) |')
            print(f'tot_n   = {tot_n_lay}')
            print(f'tot_pwv = {np.trapz(self.wvmd,self.heights):.02f}')
            
        if (not self.default_interp_method == 'rectilinear') or full_compute:
            
            print('building trees...')
            
            #self.sim_dx_, self.sim_dy_ = list(zip(*[tools.to_xy(az,el,caz,cel) for az,el,caz,cel in zip(self.sim_az_, self.sim_el_, self.obs.c_az, self.obs.c_el)]))
            self.lay_x_, self.lay_y_ = [np.real(z)  for z in self.lay_z_], [np.imag(z) for z in self.lay_z_]
            self.lay_a_, self.lay_e_ = [np.pi/2 - np.angle(z) for z in self.lay_z_], [np.pi/2 - np.arctan(np.abs(z/h)) for z,h in zip(self.lay_z_,self.heights)]
            self.lay_dx_, self.lay_dy_, self.lay_dz_ = [], [], []
            for laz,lel in zip(self.lay_a_, self.lay_e_):
                x,y = tools.to_xy(laz,lel,self.obs.c_az,self.obs.c_el)
                self.lay_dx_.append(x), self.lay_dy_.append(y), self.lay_dz_.append(x + 1j*y)

            self.static_ckdt = [sp.spatial.cKDTree(np.vstack([x.ravel(),y.ravel()]).T) for x,y in zip(self.lay_dx_, self.lay_dy_)]

        # this is too 
        if verbose:
            print('computing transition matrices...')
        sam_args = [np.c_[np.real(lz[si]),np.imag(lz[si]),lh[si]] for lz,lh,si in zip(self.lay_z_,self.lay_h_,self.sam_i_)]
        gen_args = [np.c_[np.real(gz),np.imag(gz),gh]             for gz,gh    in zip(self.gen_z_,self.gen_h_)]
        lay_v_   = [np.zeros(lz.shape) for lz in self.lay_z_]
        self.prec_ = [la.inv(tools.make_covariance_matrix(self.C,np.real(lz[si]),np.imag(lz[si]),lh[si])) for lz,lh,si in zip(self.lay_z_,self.lay_h_,self.sam_i_)]
        self.csam_ = [tools.make_covariance_matrix(self.C,np.real(gz),np.imag(gz),gh,np.real(lz[si]),np.imag(lz[si]),lh[si],auto=False) 
                 for lz,lh,si,gz,gh in zip(self.lay_z_,self.lay_h_,self.sam_i_,self.gen_z_,self.gen_h_)]
        self.cgen_ = [tools.make_covariance_matrix(self.C,np.real(gz),np.imag(gz),gh) for gz,gh in zip(self.gen_z_,self.gen_h_)]
        self.A_ = [np.matmul(csam,prec) for prec,csam in zip(self.prec_,self.csam_)]
        self.B_ = [tools.msqrt(cgen-np.matmul(A,csam.T)) for A,cgen,csam in zip(self.A_,self.cgen_,self.csam_)]
        
        print(f'done! total atmospheric points : {tot_n_lay}')
        
    def atmosphere_timestep(self,i_layer):
        
        self.lay_v_[i_layer] = np.r_[(np.matmul(self.A_[i_layer],self.lay_v_[i_layer][self.sam_i_[i_layer]])
                                 + np.matmul(self.B_[i_layer],np.random.standard_normal(self.B_[i_layer].shape[0])))[None,:],self.lay_v_[i_layer][:-1]]

    def generate_atmosphere(self,blurred=True):

        self.lay_v_ = [np.zeros(lz.shape) for lz in self.lay_z_]
        n_init_ = [2*si[0].max() for si in self.sam_i_]
        n_ts_   = [T+n_para for T,n_para in zip(self.sim_T_,self.n_para_)]
        tot_n_init, tot_n_ts = np.sum(n_init_), np.sum(n_ts_)
        self.gen_data = [np.zeros((n_ts,v.shape[1])) for n_ts,v in zip(n_ts_,self.lay_v_)]

        with tqdm(total=tot_n_init,desc='Initializing atmosphere') as prog:
            for i_layer, n_init in enumerate(n_init_):
                for i_init in range(n_init):
                    prog.update(1)
                    self.atmosphere_timestep(i_layer)
            
        with tqdm(total=tot_n_ts,desc='Generating atmosphere') as prog:
            for i_layer,n_ts in enumerate(n_ts_):
                for i_ts in range(n_ts):
                    prog.update(1)
                    self.atmosphere_timestep(i_layer)
                    self.gen_data[i_layer][i_ts] = self.lay_v_[i_layer][0]

        if blurred:
            for i_h, h in enumerate(self.heights):
                self.gen_data[i_h] = sp.ndimage.gaussian_filter1d(self.gen_data[i_h],axis=0,sigma=self.phy_sig_[i_h] / self.d_orth_[i_h])
                self.gen_data[i_h] = sp.ndimage.gaussian_filter1d(self.gen_data[i_h],axis=1,sigma=self.phy_sig_[i_h] / self.d_para_[i_h])
                
    
        
            
    def generate_cmb(self,lmax=3200,nside=2048,kind='T',blurred=True):

        cmb_ps = np.load(pkg_resources.resource_filename('aram', 'resources/cmb_spectra/act+wmap.npy'))
        if kind=='T':
            ll_C_l = cmb_ps[1]
        if kind=='E':
            ll_C_l = cmb_ps[2]
            
        print('generating CMB...')
        alm  = hp.sphtfunc.synalm(ll_C_l / (cmb_ps[0] * (cmb_ps[0] + 1)), lmax=lmax, new=True)    
        self.CMB  = hp.sphtfunc.alm2map(alm, nside, lmax=lmax, pol=True)[hp.pixelfunc.ang2pix(nside=nside, theta=np.degrees(self.RA), phi=np.degrees(self.DEC), nest=False, lonlat=True)]
        
        if blurred:
            self.CMB = sp.ndimage.gaussian_filter(self.CMB,sigma=self.array.fwhm.min() / (2*np.sqrt(2*np.log(2))))

    def simulate(self, do_atmosphere=True, 
                       do_cmb=False, 
                       do_noise=False,
                       interp_method=None,
                       split_layers=True,
                       sky_per_det=16):
        
        if interp_method == None:
            interp_method = self.default_interp_method
        
        #print(f'using {interp_method} interpolation')
        gaussian = lambda offset, sig : np.exp(-.5*np.square(offset/sig))
          
        if interp_method=='rectilinear':
            if do_atmosphere:
                self.generate_atmosphere(blurred=True)
                with tqdm(total=self.tot_sim_T,desc='Observing atmosphere') as prog:
                    h_atm_sim_data = [np.zeros((self.array.n,sim_T)) for sim_T in self.sim_T_]
                    for i_h, (sim_t_, sim_c_az_, h, p_, o_) in enumerate(zip(self.sim_t_, self.sim_c_az_, self.heights, self.lay_para_, self.lay_orth_,)):
                        for i_t, t in enumerate(sim_t_):
                            prog.update(1)
                            rel_z_ = h*np.exp(1j*(np.pi/2-(self.az_0 + sim_c_az_[i_t]))) / np.tan(self.el_0)*np.exp(1j*(3*np.pi/2 + self.wind_bearing[i_h]))
                            rel_p_, rel_o_ = np.real(rel_z_), np.imag(rel_z_)
                            ip_min, ip_max = np.where(p_ < rel_p_.min())[0][-1], np.where(p_ > rel_p_.max())[0][0] + 1
                            io_min, io_max = np.where(o_ < rel_o_.min())[0][-1], np.where(o_ > rel_o_.max())[0][0] + 1
                            i_p = (p_ > rel_p_.min()) & (p_ < rel_p_.max())
                            i_o = (o_ > rel_o_.min()) & (o_ < rel_o_.max())
                            atm_interp = sp.interpolate.RegularGridInterpolator((p_[ip_min:ip_max], o_[io_min:io_max]), self.gen_data[i_h][i_t:i_t+self.n_para_[i_h]][ip_min:ip_max,io_min:io_max])
                            h_atm_sim_data[i_h][:,i_t] = atm_interp((rel_p_, rel_o_))
            if do_cmb:
                self.generate_cmb(blurred=True)
                with tqdm(total=len(self.sim_cmb_t_),desc='Observing CMB...') as prog:
                    cmb_sim_data = np.zeros((self.array.n,len(self.sim_cmb_t_)))
                    cmb_interp = sp.interpolate.RegularGridInterpolator((self.rdy,self.rdx),self.CMB)
                    for i_t, t in enumerate(self.sim_cmb_t_):
                        prog.update(1)
                        dx,dy = tools.to_xy(*tools.ae_to_rd(self.az_0 + self.sim_cmb_c_az_[i_t], self.el_0, self.sim_cmb_mlst_[i_t], self.site.site.lat), self.ra_0, self.dec_0)
                        cmb_sim_data[:,i_t] = cmb_interp((dy,dx))
    
        if do_atmosphere:
        
            h_atm_data = np.concatenate([np.concatenate([sp.interpolate.interp1d(sim_t_,atm_sim_data[i_det],kind='slinear')(self.obs.t_)[None,:] for i_det in range(self.array.n)])[None,:,:] 
                                         for sim_t_, atm_sim_data in zip(self.sim_t_,h_atm_sim_data)])
            self.var_prof  = np.exp(-2*self.heights/self.half_height)
            self.var_prof *= np.square(self.atm_rms) / self.var_prof.sum() 
            
            if split_layers:
                atm_data = h_atm_data * np.sqrt(self.var_prof[:,None,None]) / self.el_0[None,:,None]
            else:
                atm_data = (h_atm_data * np.sqrt(self.var_prof[:,None,None])).sum(axis=0) / self.el_0[:,None]
            
        if do_cmb:
            cmb_data  = 1e-3 * np.concatenate([sp.interpolate.interp1d(self.sim_cmb_t_,cmb_sim_data[i_det],kind='slinear')(self.obs.t_)[None,:] for i_det in range(self.array.n)])

        if do_noise:
            
            # simulate white noise
            white_noise = self.array.white[:,None] * np.random.standard_normal((self.array.n,self.obs.T))
            
            # simulate pink noise
            pn_ps       = 1 / (self.obs.f_ + 1e-32); pn_ps[0] = 0 # pink spectrum
            pink_noise  = np.concatenate([np.imag(np.fft.fft(self.array.pink[i_det] * pn_ps * np.fft.ifft(np.random.standard_normal(self.obs.T),norm='ortho'),norm='ortho'))[None,:] 
                                          for i_det in range(self.array.n)])
            noi_data    = pink_noise + white_noise
            
        data = {}
        if do_atmosphere:
            data['atmosphere'] = atm_data
        if do_cmb:
            data['cmb'] = cmb_data
        if do_noise:
            data['noise'] = noi_data

        return data
    


