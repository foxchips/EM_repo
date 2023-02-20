# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 23:47:48 2020
#===================2D binary grating===========================
Blaze grating+coating to test Ax libraries on AWS and noob practice for sfox
Incident planewave from Glass substrate
- I discover and steal a lot of this from https://meep.readthedocs.io/en/latest/Python_Tutorials/Near_to_Far_Field_Spectra/#diffraction-spectrum-of-a-finite-binary-grating

@author: sfox1
"""

from __future__ import absolute_import, division, print_function, unicode_literals
import meep as mp
import numpy as np
import argparse
import math
import cmath
from mpi4py import MPI
import json
#from materials import *
#import matplotlib.pyplot as plt
from scipy.signal import find_peaks


rank = MPI.COMM_WORLD.rank

def main(args):
    experiment = args.experiment_name
    arm = args.arm_name
    trial = args.trial_index
    param_dict = json.loads(args.param_dict)

    gp = param_dict.get("gp")   #pitch (um)
    gh = param_dict.get("gh")   #height of element (um)
    ghc = param_dict.get("ghc")   #thickness of coating (um)
    gdc = param_dict.get("gdc")  #blaze ratio
    theta_i = param_dict.get("theta_i")  #incident angle (degrees)

#=======================Optical Properties===================================
    n_Glass = 1.7    #RI of high RI Glass
    Glass = mp.Medium(index= n_Glass)

    n_TiO2 = 2.36    #RI of TiO2 coating
    TiO2 = mp.Medium(index = n_TiO2)

#====================Simulation parameters===================================    
    lambda_cen = 0.62   #center source wavelength
    fcen = 1/lambda_cen  #source center frequency
    df = 0.05*fcen    #source frequency width
    nfreq = 1  #no of frequency points for DFT
    dpml = 1.0   #thickness of absorber/PML at +/-x boundaries
    resolution = 60    #resolution of simulation
    theta_inc = np.radians(theta_i)   # theta_inc is relative to normal of slab
    src_cmpt = mp.Ez   #polarisation of source

#====================Geometry parameters=====================================
    dsub = 3.0  #thickness of substrate
    dpad = 3.0  #distance from top of grating element to +x PML (superstrate thickness)
    
    Lx = 2 * dpml + dsub + gh + ghc + dpad  #cell size dimensions
    Ly = gp

    cell = mp.Vector3(Lx, Ly, 0)

#vertices to define blaze element
    vtx_coating = [mp.Vector3(-0.5*Lx+dpml+dsub,0.5*Ly),
             mp.Vector3(-0.5*Lx+dpml+dsub+ghc, 0.5*Ly),
             mp.Vector3(-0.5*Lx+dpml+dsub+gh+ghc,0.5*Ly-(gp*(1-gdc))),
             mp.Vector3(-0.5*Lx+dpml+dsub+ghc,-0.5*Ly),
             mp.Vector3(-0.5*Lx+dpml+dsub,-0.5*Ly)]

    vtx_Glass = [mp.Vector3(-0.5*Lx+dpml+dsub,-0.5*Ly),
             mp.Vector3(-0.5*Lx+dpml+dsub+gh,-0.5*Ly+gdc*Ly),
             mp.Vector3(-0.5*Lx+dpml+dsub,0.5*Ly)]

#===========================Set up geometry===================================
    geometry_bg = []  
   
    geometry_grat = [mp.Prism(vertices = vtx_coating, material=TiO2, height = mp.inf),
                 mp.Prism(vertices = vtx_Glass, material=Glass, height = mp.inf), 
                 mp.Block(material=Glass, size=mp.Vector3(dpml+dsub,mp.inf,mp.inf), center=mp.Vector3(-0.5*Lx+0.5*(dpml+dsub),0,0))]   

#===========================Boundary conditions===============================
    boundary_layers = [mp.PML(thickness=dpml,direction=mp.X)]       #boundary layers
    k = mp.Vector3(fcen*n_Glass).rotate(mp.Vector3(z=1), theta_inc)
    if theta_i == 0:
        k = mp.Vector3(0,0,0)

#===========================Set source===============================
    def pw_amp(k,x0):
        def _pw_amp(x):
            return cmath.exp(1j*2*math.pi*k.dot(x+x0))
        return _pw_amp

    src_x = -0.5*Lx+dpml+0.3*dsub
    src_r0 = mp.Vector3(src_x,0,0)

    sources = [mp.Source(mp.GaussianSource(fcen, fwidth=df), 
                     component=src_cmpt, 
                     center=src_r0, 
                     size=mp.Vector3(0,Ly,0),
                     amp_func=pw_amp(k,src_r0))]

#=========================Define computation cell=============================
    sim_bg = mp.Simulation(resolution = resolution,
                       cell_size=cell,
                       boundary_layers=boundary_layers,
                       geometry=geometry_bg,
                       sources=sources,
                       k_point=k,
                       default_material=Glass)    

    sim = mp.Simulation(resolution=resolution,
                    cell_size=cell,
                boundary_layers=boundary_layers,
                geometry=geometry_grat,
                sources=sources,
                k_point=k)


#===================Farfield parameters======================================
    ff_distance = 1e3      # far-field distance from near-field monitor
    ff_angle = 88          # far-field cone angle
    ff_npts = 1500         # number of far-field points
    nperiods = 10          # number of unit cells in finite periodic grating
    ff_length = ff_distance*math.tan(math.radians(ff_angle))
    ff_res = ff_npts/ff_length
    n2f_pt = mp.Vector3(-0.5*Lx+dpml+0.6*dsub,0,0) 

    n2f_obj_bg = sim_bg.add_near2far(fcen, 0, nfreq, mp.Near2FarRegion(center=n2f_pt, size=mp.Vector3(y=Ly), weight = -1), nperiods=nperiods)   #get near fields in norm run
    n2f_obj = sim.add_near2far(fcen, 0, nfreq, mp.Near2FarRegion(center=n2f_pt, size=mp.Vector3(y=Ly), weight=-1), nperiods=nperiods)  #near fields in grating run


#===================NORMALISATION RUN=========================================
    sim_bg.run(until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt, n2f_pt, 1e-9))
    input_flux = sim_bg.get_near2far_data(n2f_obj_bg)   #save data from normalization run
    ff_source = sim_bg.get_farfields(n2f_obj_bg, ff_res, center=mp.Vector3(-ff_distance,0.5*ff_length), size=mp.Vector3(y=ff_length))    #to calculate the enhancement

#============================RUN WITH GRATING=================================
    sim.load_minus_near2far_data(n2f_obj, input_flux)   #subtract incident fields
    sim.run(until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt, n2f_pt, 1e-9))
    ff_unitcell = sim.get_farfields(n2f_obj, ff_res, center=mp.Vector3(-ff_distance,0.5*ff_length), size=mp.Vector3(y=ff_length))
    
#=======================Calculate angular far field================================
    def calc_angle(m,lambda_cen,n,gp,theta_i):
        "Calculate angle of diffracted m'th order"
        lambda_n = lambda_cen/n_Glass
        theta_m=math.degrees(math.asin(math.sin(math.radians(theta_i))+m*lambda_n/gp))   #angle of mth order
        return theta_m
    
    def find_nearest(a, a0):
        "Element in nd array `a` closest to the scalar value `a0`"
        idx = np.abs(a - a0).argmin()
        return idx

    order = 1   #mth order of interest
    theta_m=calc_angle(order,lambda_cen,n_Glass,gp,theta_i)   #angle of mth order
    #N.B. angular space is not evenly spaced
    freq = mp.get_near2far_freqs(n2f_obj)
    wvl = np.divide(1,freq)
    ff_lengths = np.linspace(0,ff_length,ff_npts)
    angles = [math.degrees(math.atan(f)) for f in ff_lengths/ff_distance]
    F_rel = np.absolute(ff_unitcell['Ez'])**2/np.absolute(ff_source['Ez'])**2   #relative enhancement
    F_rel_arr = np.array(F_rel)
    angles_arr = np.array(angles)
    peaks, _ = find_peaks(F_rel, height=1e6)    #find peaks
    angles_pks = angles_arr[peaks]
    indx=find_nearest(angles_pks, theta_m)
    ind_pk=peaks[indx]
    F_rel_order = F_rel_arr[ind_pk]     #output relative enhancement of order
    print("wavelength:, {:.3f}".format(wvl[0]))
    print("max angle range:, {:.3f}".format(angles[-1]))
    print("relative enhancement at {:.3f} degrees:, {:.1f}".format(angles[ind_pk],F_rel_order))
    #print(F_rel_arr[ind_pk])
    # np.savetxt('F_rel.txt',F_rel_arr,fmt='%.3f')
    # np.savetxt('angles.txt',angles_arr,fmt='%.3f')
    # np.savetxt('peak.txt',F_rel_arr[ind_pk],fmt='%.3f')
    
    #Save to database
    metric = F_rel_order
    data = {                                #data to be serialized
        "experiment_name": experiment,  
        "arm_name": arm,
        "trial_index": trial,
        "data": {"metric": metric,},
        "raw_data": {
            'relative enhancement factor' : F_rel_arr.tolist(),
            'thetas' : angles_arr.tolist(),
            'center wavelength': lambda_cen.tolist(),
        }
    }

    with open(f'data_{trial}.txt', 'w') as outfile: #write information to outfile whcih will store as bytes
        json.dump(data, outfile)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--experiment_name",
        type=str,
    default="experiment",
        help="experiment_name for AEMySQLDataWriter",
    )
    parser.add_argument(
        "--arm_name", type=str, default="arm", help="arm_name for AEMySQLDataWriter"
    )
    parser.add_argument(
        "--trial_index", type=int, default=0, help="trial_index for AEMySQLDataWriter"
    )
    parser.add_argument(
        "--param_dict", type=str, default="{}", help="parameters to optimize over"
    )
    args = parser.parse_args()
    main(args)