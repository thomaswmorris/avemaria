# -*- coding: utf-8 -*-
from . import helpers
import numpy as np

def make_array(array_shape,max_fov,max_n_det):

    valid_array_types = 'flower', 'hex', 'square' 
    
    if array_shape=='phyllotaxis':
        phi = np.pi*(3.-np.sqrt(5.))  # golden angle in radians
        dzs = np.zeros(max_n_det).astype(complex)
        for i in range(max_n_det):
            dzs[i] = np.sqrt((i / (max_n_det - 1)) * 2 ) *np.exp(1j*phi*i)
        od = np.abs(np.subtract.outer(dzs,dzs))
        dzs *= max_fov / od.max()
        return np.real(dzs), np.imag(dzs)
    if array_shape=='hex':
        HEX  = lambda n : 3*n*(n+1) + 1 ; hex_layers = 0
        while HEX(hex_layers+1) < max_n_det: hex_layers += 1 
        dr   = max_fov/ (2 * hex_layers)
        dzs  = np.zeros(1).astype(complex)
        angs = np.arange(0,2*np.pi,np.pi/3)+ np.pi/6
        for ilay in range(hex_layers):
            for iz,z in enumerate(dzs):
                for iang,ang in enumerate(angs):
                    z_ = np.round(z+dr*np.exp(1j*ang),6)
                    if np.amin(np.abs(z_-dzs) > dr/4):
                        dzs = np.append(dzs,z_)
        return np.real(dzs), np.imag(dzs)
    if array_shape=='square':
        dxy_ = np.linspace(-max_fov,max_fov,int(np.floor(np.sqrt(max_n_det))))/(2*np.sqrt(2))
        DX, DY = np.meshgrid(dxy_,dxy_)
        return DX.ravel(), DY.ravel()
    
    except ValueError:
        print('Please specify a valid array type. Valid array types are:\n'+
              '\n'.join(valid_array_types))
        
    
