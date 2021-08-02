# -*- coding: utf-8 -*-
import numpy as np
import numpy.linalg as la
import scipy as sp

# ================ ARRAY ================

def make_array(array_shape,max_fov,max_n_det):

    valid_array_types = ['flower', 'hex', 'square']
    
    if array_shape=='flower':
        phi = np.pi*(3.-np.sqrt(5.))  # golden angle in radians
        dzs = np.zeros(max_n_det).astype(complex)
        for i in range(max_n_det):
            dzs[i] = np.sqrt((i / (max_n_det - 1)) * 2 ) *np.exp(1j*phi*i)
        od = np.abs(np.subtract.outer(dzs,dzs))
        dzs *= max_fov / od.max()
        return np.real(dzs), np.imag(dzs)
    if array_shape=='hex':
        nsides   = 6
        vertices = np.exp(1j*(np.pi/2-np.linspace(0,2*np.pi,nsides+1)))
        i_layer  = 0
        zs = np.array([0])
        while len(zs) < max_n_det:
            for i_edge,(i_vert,j_vert) in enumerate(zip(vertices,vertices[1:])):
                t_ = np.linspace(0,1,i_layer+1)[:-1]
                points = i_layer*(t_*i_vert + (1-t_)*j_vert)
                zs = np.append(zs,points)
            i_layer += 1
        zs = zs[:max_n_det]
        zs *= max_fov/np.abs(np.subtract.outer(zs,zs)).max()
        return np.real(zs), np.imag(zs)
    if array_shape=='square':
        dxy_ = np.linspace(-max_fov,max_fov,int(np.floor(np.sqrt(max_n_det))))/(2*np.sqrt(2))
        DX, DY = np.meshgrid(dxy_,dxy_)
        return DX.ravel(), DY.ravel()
    
    raise ValueError('Please specify a valid array type. Valid array types are:\n'+
              '\n'.join(valid_array_types))

# ================ STATS ================

msqrt  = lambda M : [np.matmul(u,np.diag(np.sqrt(s))) for u,s,vh in [la.svd(M)]][0]
matern = lambda r,r0,nu : 2**(1-nu)/sp.special.gamma(nu)*sp.special.kv(nu,r/r0+1e-10)*(r/r0+1e-10)**nu
    
gaussian_beam = lambda z, w0, l, n : np.sqrt(1/np.square(z) + np.square(l) / np.square(w0 * np.pi * n))

def make_covariance_matrix(C,x0,y0,z0,x1=None,y1=None,z1=None,auto=True):
    if auto:
        n = len(x0); i,j = np.triu_indices(n,1)
        o = C(np.sqrt(np.square(x0[i] - x0[j]) + np.square(y0[i] - y0[j]) + np.square(z0[i] - z0[j])))
        c = np.empty((n,n)) 
        c[i,j],c[j,i] = o,o
        c[np.eye(n).astype(bool)] = C(0)
    if not auto:
        n = len(x0); i,j = np.triu_indices(n,1)
        c = C(np.sqrt(np.square(np.subtract.outer(x0,x1))
                    + np.square(np.subtract.outer(y0,y1))
                    + np.square(np.subtract.outer(z0,z1))))
    return c

# ================ POINTING ================

def to_xy(az, el, c_az, c_el):
    ground_X, ground_Y, ground_Z = np.sin(az-c_az)*np.cos(el), np.cos(az-c_az)*np.cos(el), np.sin(el)
    return np.arcsin(ground_X), np.arcsin(-np.real((ground_Y+1j*ground_Z)*np.exp(1j*(np.pi/2-c_el))))

def from_xy(dx, dy, c_az, c_el):
    ground_X, Y, Z = np.sin(dx+1e-16), -np.sin(dy+1e-16), np.cos(np.sqrt(dx**2+dy**2))
    gyz = (Y+1j*Z)*np.exp(-1j*(np.pi/2-c_el))
    ground_Y, ground_Z = np.real(gyz), np.imag(gyz)
    return (np.angle(ground_Y+1j*ground_X) + c_az) % (2*np.pi), np.arcsin(ground_Z)

def ae_to_rd(az, el, lst, lat):
    NP_X, NP_Y, NP_Z = np.sin(az)*np.cos(el), np.cos(az)*np.cos(el), np.sin(el)
    lat_rot_YZ = (NP_Y + 1j*NP_Z)*np.exp(1j*(np.pi/2-lat))
    lat_rot_Y, globe_Z = np.real(lat_rot_YZ), np.imag(lat_rot_YZ)
    return np.arctan2(NP_X,-lat_rot_Y) + lst, np.arcsin(globe_Z)
    
def rd_to_ae(ra, de, lst, lat):
    NP_X, globe_Y, globe_Z = np.sin(ra-lst)*np.cos(de), -np.cos(ra-lst)*np.cos(de), np.sin(de)
    lat_rot_YZ = (globe_Y + 1j*globe_Z)*np.exp(-1j*(np.pi/2-lat))
    NP_Y, NP_Z = np.real(lat_rot_YZ), np.imag(lat_rot_YZ)
    return np.arctan2(NP_X,NP_Y), np.arcsin(NP_Z)
