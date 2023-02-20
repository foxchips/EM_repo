# -*- coding: utf-8 -*-
"""
====================Green cavity backed antenna Nano LED @530nm===========================
Trapezoidal cavity for enhanced LEE
-------------------------------------
Update: ?? include near to far DFT
Update: ?? include incoherent source that requires several simulations to be averaged
Update: 15/10/2020 tidy script - noob practice (sfox)

@author: sfox1
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import meep as mp
import numpy as np
import argparse
from mpi4py import MPI
import h5py
import json
from materials import *
import matplotlib.pyplot as plt

rank = MPI.COMM_WORLD.rank

def main(args):
    experiment = args.experiment_name
    arm = args.arm_name
    trial = args.trial_index
    param_dict = json.loads(args.param_dict)

    l = param_dict.get("l") #dimension of along x
    d = param_dict.get("d") #dimension along z
    w = param_dict.get("w") #dimension along y 
    slant = param_dict.get("slant") #sidewall slope angle (degrees) to create trapezoid

#=======================Optical Properties===================================
#from Anurag 06/24
    n_nGaN =2.42 #RI of n-type GaN
    n_pGaN = 2.42    #RI of p-type GaN
    n_QW=2.55    #RI of active region (InGaN)
    n_glass = 1.5    #RI of glass (SiO2?)
    n_Al2O3 = 1.77   #RI of n-contact insulator (Al2O3)

    nGaN = mp.Medium(index=n_nGaN)
    pGaN = mp.Medium(index=n_pGaN)
    Sapphire = mp.Medium(index =n_Al2O3)
    Glass = mp.Medium(index= n_glass)
    qw_GaN = mp.Medium(index = n_QW)

#used for checking geometry
    mat1 = mp.Medium(index = 3) #mp.Medium(epsilon=9)
    mat2 = mp.Medium(index = 6)
    mat3 = mp.Medium(index = 3)
    mat4 = mp.Medium(index = 7)

#====================Simulation parameters===================================
    um_scale = 1.0  
    eV_um_scale = um_scale//1.23984193  #for material libaray to get 1/um
    lambda_min = 0.51   #min source wavelength
    lambda_max = 0.54   #max source wavelength
    fmin = 1/lambda_max     #min source frequency
    fmax = 1/lambda_min     #max source frequency
    fcen = 0.5*(fmin+fmax)  #source center frequency
    df = (fmax-fmin)    #source frequency width
    nfreq = 20  #no of frequency points for DFT
    tpml = lambda_max/4.0   #thickness of absorber/PML
    res = 80    #resolution of simulation
    
#====================Geometry parameters=====================================
    t_Al2O3 = 0.01  #thickness of sidewall insulator
    h_Al2O3 = 0.9*d     #height of insulator (z)
    w_glass = 0.3   #width of glass
    t_glass = 2*t_Al2O3     #thickness of glass
    metal_pad = 0.1 #thickness of metal contact, 100nm should be enough for 100% reflectivity with Ag
    z_qw = d/2.0    #center position of active region (along z)
    t_qw = 0.05      #thickness of QW region
    h_ntype = 0     #thickness of n-type semiconductor
    h_ptype = d     #thickness of p-type semiconductor
        
    sx = l + 2*np.tan(np.deg2rad(slant))*d + 2*t_Al2O3 #width of semiconductor LED
    sy = w + 2*np.tan(np.deg2rad(slant))*d + 2*t_Al2O3
    if d>t_glass:   
        sz = d + metal_pad
    else:   #include any SiO2 planarisation on top
        sz = t_glass
    
    pad_x = lambda_max/2.0  #padding between structure and PML
    pad_y = lambda_max/2.0
    pad_z = lambda_max/2.0
    
    Lx = 2 * tpml + sx + pad_x  #cell size dimensions
    Ly = 2 * tpml + sy + pad_y
    Lz = 2 * tpml + sz + pad_z
    cell = mp.Vector3(Lx, Ly, Lz)
    boundary_layers = [mp.Absorber(tpml)]   #use absorber because more stable with evanescent fields
    dx=1/res
    
#===========================Set up geometry===================================

#create trapezoid by via many layers
#TBC: use mp.Prism
    p_type = []     #p-GaN
    t_block =  h_ptype
    N_layers = np.round(t_block*res)+1
    z0_p = Lz/2.0 - tpml - metal_pad - t_block/2.0
    dh = t_block/N_layers
    dx = dh*np.tan(np.deg2rad(slant))
    wx_p = l  + np.tan(np.deg2rad(slant))*t_block
    wy_p = w + np.tan(np.deg2rad(slant))*t_block
    for ii in np.arange(N_layers):
        p_type.append(mp.Block(center=mp.Vector3(0,0, z0_p-t_block/2.0 + dh/2.0+ii*dh ), size = mp.Vector3(wx_p - dx*ii, wy_p - dx*ii, dh), material = pGaN))

    n_type = []     #n-GaN
    t_block =  h_ntype
    N_layers = np.round(t_block*res)+1
    z0 =  Lz/2.0 - tpml - metal_pad - h_ptype - h_ntype/2.0
    dh = t_block/N_layers
    dx = dh*np.tan(np.deg2rad(slant))
    wx_n = wx_p + np.tan(np.deg2rad(slant))*t_block
    wy_n = wy_p + np.tan(np.deg2rad(slant))*t_block
    for ii in np.arange(N_layers):
        n_type.append(mp.Block(center=mp.Vector3(0,0, z0-t_block/2.0 + dh/2.0+ii*dh ), size = mp.Vector3(wx_n - dx*ii, wy_n - dx*ii, dh), material = nGaN))


    n_contact_insulator = []    #insulator
    t_block =  h_Al2O3
    N_layers = np.round(t_block*res/2.0)+1
    z0_ni =  z0_p - (h_ptype - h_Al2O3)/2.0
    dh = t_block/N_layers
    dx = dh*np.tan(np.deg2rad(slant))
    wx_i = wx_n + 2*t_Al2O3 #+ np.tan(np.deg2rad(slant))*t_block
    wy_i = wy_n + 2*t_Al2O3 #+ np.tan(np.deg2rad(slant))*t_block
    for ii in np.arange(N_layers):
        n_contact_insulator.append(mp.Block(center=mp.Vector3(0,0, z0_ni-t_block/2.0 + dh/2.0+ii*dh ), size = mp.Vector3(wx_i - dx*ii, wy_i - dx*ii, dh), material = Sapphire))

    qw = []     #active region
    t_block =  t_qw
    N_layers = np.round(t_block*res)+1
    z0_qw = z0_p - h_ptype/2.0 + t_qw/2.0
    dh = t_block/N_layers
    dx = dh*np.tan(np.deg2rad(slant))
    wx_qw = wx_n - np.tan(np.deg2rad(slant))*(h_ntype - t_qw/2.0)
    wy_qw = wy_n - np.tan(np.deg2rad(slant))*(h_ntype - t_qw/2.0)
    for ii in np.arange(N_layers):
        qw.append(mp.Block(center=mp.Vector3(0,0, z0_qw-t_block/2.0 + dh/2.0+ii*dh ), size = mp.Vector3(wx_qw - dx*ii, wy_qw - dx*ii, dh), material = qw_GaN))


    t_nGaN = Lz - 2*tpml - metal_pad - d
    z0_nGaN = Lz/2 - tpml -d - metal_pad - t_nGaN/2.0
    substrate = [mp.Block(center = mp.Vector3(0,0,0), size = mp.Vector3(Lx,Ly,Lz-2*tpml-2/res),material = Ag ),
            mp.Block(center = mp.Vector3(0,0, Lz/2 - tpml - d - metal_pad - t_nGaN/2.0 ), size = mp.Vector3(Lx,Ly,t_nGaN),material = nGaN), # substrate
                # metal p-contact
#                   mp.Block(center = mp.Vector3(0,0,d/2.0), size = mp.Vector3(l,w,d),material = mp.air),
               ]

   n_type_spacer = []
   t_block =  t_glass+2/res
   N_layers = np.round(t_block*res)+1
   z0_ns = z0_qw - t_qw/2.0 - t_glass/2.0 - 1/res
   dh = t_block/N_layers
   dx = dh*np.tan(np.deg2rad(slant))
   wx_spacer = wx_n + np.tan(np.deg2rad(slant))*t_block
   wy_spacer = wy_n + np.tan(np.deg2rad(slant))*t_block
   for ii in np.arange(N_layers):
   n_type_spacer.append(mp.Block(center=mp.Vector3(0,0, z0_ns-t_block/2.0 + dh/2.0+ii*dh ), size = mp.Vector3(wx_n - dx*ii, wy_n - dx*ii, dh), material =nGaN))



   spacer = [mp.Block(center = mp.Vector3(0, 0, z0_qw- t_qw/2.0 - t_glass/2.0), size = mp.Vector3(Lx, Ly, t_glass),material = Sapphire), # glass spacer
#                   mp.Block(center = mp.Vector3(0,0,-t_glass/2.0),size = mp.Vector3(wx_n_type,wy_n_type,t_glass),material = nGaN), # replace nGaN lost from glass spacer
        ]

   green_pixel = substrate + n_contact_insulator + n_type + p_type + spacer + qw + n_type_spacer





   geo_bg = [mp.Block(center = mp.Vector3(0,0,-Lz/4), size = mp.Vector3(1e20,1e20,Lz),material = nGaN),
          mp.Block(center = mp.Vector3(0,0,Lz/4.0), size = mp.Vector3(1e20,1e20,Lz),material =pGaN ),
          mp.Block(center = mp.Vector3(0,0,0), size = mp.Vector3(1e20,1e20,t_qw),material =qw_GaN ),]
