# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg as la
import scipy as sp
from datetime import datetime
from datetime import timezone
import os

print(os.getcwd())

from maria import tools


# ================ ARRAY ================

       #self.build_focal_plane_layers()
        
atmosphere_config = {'n_layers'        : 10,         # how many layers to simulate, based on the integrated atmospheric model 
                    'min_depth'        : 500,      # the height of the first layer 
                    'max_depth'        : 10000,      # 
                    'atmosphere_rms'   : 50,  
                    'outer_scale'      : 500}


time_ = np.linspace(0,600,6000)
 
phase = 60

focal_x = np.radians(5)*np.cos(1.1*2*np.pi*time_/phase)
focal_y = np.radians(5)*np.cos(2*np.pi*time_/phase)

focal_azim, focal_elev = tools.from_xy(focal_x,focal_y,0,np.pi/4)


pointing_config = {'scan_type' : 'lissajous_daisy',
                    'duration' : 60,'samp_freq' : 20,
                 'center_azim' : -45, 'center_elev' : 30, 
                       'throw' : 5, 'r_period' : 21, 'p_period' : 29}


pointing_config = {'scan_type' : 'CES',
            'duration' : 600,'samp_freq' : 20,
         'center_azim' : 90, 'center_elev' : 60, 
            'az_throw' : 45, 'az_speed' : 1.5}

pointing_config = {'scan_type' : 'lissajous_box',
                    'duration' : 60,'samp_freq' : 10,
                 'center_azim' : -45, 'center_elev' : 30, 
                     'x_throw' : 15, 'x_period' : 21,
                     'y_throw' : 15, 'y_period' : 29}

array_config = {'shape' : 'flower',
                    'n' : 10000,      
                  'fov' : 2,
                 'band' : 1.5e11}    

beams_config = {'optical_type' : 'diff_lim',
                'primary_size' : 50,
                  'beam_model' : 'top_hat',
                    'beam_res' : .5 }    

site_config = {'site' : 'ACT',
               'time' : datetime.now(timezone.utc).timestamp(),
 'weather_gen_method' : 'random',
             'region' : 'atacama' } 

heights = np.linspace(0,10000,100)


import matplotlib as mpl
import matplotlib.pyplot as plt


mpl.rcParams['figure.dpi'] = 256
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})

def equalize(ax):
    xls,yls = ax.get_xlim(),ax.get_ylim()
    x0,y0 = np.mean(xls), np.mean(yls)
    r = np.maximum(-np.subtract(*xls),-np.subtract(*yls))/2
    ax.set_xlim(x0-r,x0+r); ax.set_ylim(y0-r,y0+r)
    return ax

new = True
sim = True

if new:
    
    #tm = model(verbose=True)

    tm = model(atmosphere_config=atmosphere_config,
               pointing_config=pointing_config,
               beams_config=beams_config,
               array_config=array_config,
               site_config=site_config,
               verbose=True)
    
    if sim:
        data = tm.sim()

if sim:

    
    fig,ax = plt.subplots(1,1,figsize=(8,8))
    
    ax.pcolormesh(np.degrees(tm.X[-1]),
                  np.degrees(tm.Y[-1]),
                  tm.vals[-1],shading='none',cmap='RdBu_r')
    
    ax.scatter(np.degrees(np.real(tm.pointing.theta_edge_z[-1]).T),
               np.degrees(np.imag(tm.pointing.theta_edge_z[-1]).T),s=1e-1,c='k')
    
    ax.plot(np.degrees(np.real(tm.pointing.focal_theta_z[-1])),
            np.degrees(np.imag(tm.pointing.focal_theta_z[-1])),c='k')
    
    equalize(ax)
    
    fig,axes = plt.subplots(2,1,figsize=(12,8))
    
    for idet,det in enumerate(np.random.choice(data.shape[0],16,replace=False)):
        axes[0].plot(tm.pointing.time,data[idet])
        axes[1].plot(tm.pointing.time,data[idet]-data.mean(axis=0))
        
    fig,axes = plt.subplots(1,2,figsize=(12,8))
        
    nf = 256
    
    fmids = np.geomspace(1e-3,1e1,nf)
    rfreq = np.exp(np.gradient(np.log(fmids))).mean()
    fbins = np.append(fmids/np.sqrt(rfreq),fmids[-1]*np.sqrt(rfreq)) 
    
    freq = np.fft.fftfreq(tm.pointing.nt,tm.pointing.dt)
    
    ps   = np.square(np.abs(np.fft.fft(data * np.hanning(data.shape[-1])[None,:],axis=-1)))
    mps  = ps.mean(axis=0)
    bmps = sp.stats.binned_statistic(freq,mps,bins=fbins,statistic='mean')[0]
    nn   = ~np.isnan(bmps)
    
    axes[0].plot(fmids[nn],bmps[nn])
    axes[0].plot(fmids[nn],1e2*fmids[nn]**(-8/3))
    axes[0].loglog()
    
    axes[1].scatter(tm.pointing.time,np.degrees(np.abs(tm.atmosphere.aam)))
    
fig,ax = plt.subplots(1,1,figsize=(6,6))

ax.scatter(np.degrees(tm.array.x),
           np.degrees(tm.array.y),
           s=1e0,c='k')

ax.plot(np.degrees(tm.array.edge_x).T,
        np.degrees(tm.array.edge_y).T)

pt = np.linspace(0,2*np.pi,64)

beam_plot_height = int(np.sqrt(len(tm.atmosphere.depths)))
beam_plot_length = int(np.ceil(len(tm.atmosphere.depths)/beam_plot_height))

if tm.array.n < 200:
    
    fig,axes = plt.subplots(beam_plot_height,beam_plot_length,
                            figsize=(2*beam_plot_length,2*beam_plot_height),constrained_layout=True)#,sharex=True,sharey=True)
    
    fig.suptitle(f'D = {tm.beams.aperture:.02f}m')
    for ilay,depth in enumerate(tm.atmosphere.depths):
        
        iy = ilay % beam_plot_length
        ix = int(ilay / beam_plot_length)
        
        axes[ix,iy].set_title(f'z = {depth:.02f}m')
        axes[ix,iy].plot(np.degrees(tm.array.x+tm.beams_waists[ilay]/depth*np.cos(pt)[:,None]),
                         np.degrees(tm.array.y+tm.beams_waists[ilay]/depth*np.sin(pt)[:,None]),lw=.5)
    




