import meep as mp
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 20})
import argparse








def pw_amp(k, x0):
    def _pw_amp(x):
        return np.exp(1j*k.dot(x+x0))
    return _pw_amp


res = 20


lamda_min = .4
lamda_max = 1.5
fmin = 1/lamda_max
fmax = 1/lamda_min
fcen = (fmin+fmax)/2.0
df = fmax-fmin
Nfreq = 1000

pad = 1
tpml = lamda_max/4.0
pml_layers = [mp.Absorber(tpml, direction=mp.X)]


theta_inc = 0 * np.pi/180  ## theta_inc is relative to normal of slab
kdir = mp.Vector3(np.cos(theta_inc), -1*np.sin(theta_inc))  ## direction of k 
k0 = kdir.unit().scale(2*np.pi*fcen)  ## give k correct length
k = kdir.unit().scale(fcen)

src_cmpt = mp.Ey


n1 = 1.74
n2 = 1.56

mat1 = mp.Medium(index = n1)
mat2 = mp.Medium(index = n2)
glass = mp.Medium(index = n2)
#slab_w = 1.5
#slab_h = 2*tpml+Ly

##### Etalon Paramaters #####

t_grat = 5
 
pitch = .31
theta = 18*np.pi/180
Np = np.round(t_grat/pitch/np.cos(theta)) 

Px = pitch / np.cos(theta)
Py = pitch / np.sin(theta)



def appodize_func(p):
     phi = 2*np.pi*p.x/Px + 2*np.pi*p.y/Py
     sig=t_grat/4.0
     n = n2 + (n1-n2)*np.cos(phi)*np.exp(-(p.x**2)/(2*sig**2))
     return mp.Medium(index=n)
def no_appodize_func(p):
     phi = 2*np.pi*p.x/Px + 2*np.pi*p.y/Py
     sig=t_grat/4.0
     n = n2 + (n1-n2)*np.cos(phi)
     return mp.Medium(index=n)

#Lx = 2*tpml+2*pad+Np*pitch
Lx = 2*tpml + 2*pad + t_grat
Ly = Py



a1 = mp.Vector3(1,0,0)
a2 = mp.Vector3(0,1,0)
a3 = mp.Vector3(0,0,1)


dx = t_grat/Np
dy = Ly/Np

appodized = [mp.Block(center=mp.Vector3(0,0), size = mp.Vector3(t_grat,1e20,1e20),material=appodize_func)]

no_appodize = [mp.Block(center=mp.Vector3(0,0), size = mp.Vector3(t_grat,1e20,1e20),material=no_appodize_func)]

#for jj in np.arange(Np):
#     etalon.append(mp.Block(center=mp.Vector3(-t_grat/2+jj*dx,-Ly/2.0+jj*dy), size=mp.Vector3(pitch/2.0, 1e20, 1e20),e1=a1.rotate(a3,theta), e2=a2.rotate(a3,theta), material=mat1))

spheres = []
for jj in np.arange(Np):
     spheres.append(mp.Sphere(center=mp.Vector3(-t_grat/2+jj*dx,-Ly/2.+jj*dy), radius=0.1, material=glass))

#etalon1 = [mp.Block(center=mp.Vector3(0,0), size=mp.Vector3(t_grat,L))]

w = (Lx - t_grat)/2
appodized.append(mp.Block(center = mp.Vector3(-t_grat/2.0 - w/2), size = mp.Vector3(w,Ly),material=glass))
appodized.append(mp.Block(center = mp.Vector3(t_grat/2 + w/2), size = mp.Vector3(w,Ly),material=glass))

no_appodize.append(mp.Block(center = mp.Vector3(-t_grat/2.0 - w/2), size = mp.Vector3(w,Ly),material=glass))
no_appodize.append(mp.Block(center = mp.Vector3(t_grat/2 + w/2), size = mp.Vector3(w,Ly),material=glass))

cell = mp.Vector3(Lx, Ly, 0)
#geo = [mp.Block(mp.Vector3(slab_w, slab_h, 1e20), center=mp.Vector3(0,0), material=mat1)]
geo_app = appodized #+ spheres
geo_no_app = no_appodize


geo_background = []

src_r0 = mp.Vector3(-Lx/2.0+1.001*tpml)
sources = [mp.Source(mp.GaussianSource(fcen,fwidth=df), component=src_cmpt, \
        center=src_r0, size=mp.Vector3(0,Ly), amp_func=pw_amp(k0,mp.Vector3(x=-Lx/2+1.001*tpml)))]


sim_background = mp.Simulation(resolution = res, cell_size = cell,
        boundary_layers=pml_layers, geometry = geo_background, sources = sources, k_point=k, default_material=glass, ensure_periodicity=False, Courant=0.2)
    
sim_app = mp.Simulation(resolution = res, cell_size = cell,
        boundary_layers=pml_layers, geometry = geo_app, sources = sources, k_point=k, default_material=glass, ensure_periodicity=False, Courant=0.2)

sim_no_app = mp.Simulation(resolution = res, cell_size = cell,
        boundary_layers=pml_layers, geometry = geo_no_app, sources = sources, k_point=k, default_material=glass, ensure_periodicity=False, Courant=0.2)



# 
## add flux monitors
refl_bg = sim_background.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(-Lx/2+1.1*tpml,0), size=mp.Vector3(0,Ly), direction=mp.X))
trans_bg = sim_background.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(Lx/2-1.1*tpml,0), size = mp.Vector3(0,Ly), direction=mp.X))

sim_background.run( #mp.at_time(10, mp.output_efield_z),
        until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt,src_r0, 1e-7))
sim_background.save_flux('refl_bg', refl_bg)
sim_background.save_flux('trans_bg',trans_bg)

refl_app = sim_app.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(-Lx/2+1.1*tpml,0), size=mp.Vector3(0,Ly), direction=mp.X))
trans_app = sim_app.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(Lx/2-1.1*tpml,0), size = mp.Vector3(0,Ly), direction=mp.X))

refl_noapp = sim_no_app.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(-Lx/2+1.1*tpml,0), size=mp.Vector3(0,Ly), direction=mp.X))
trans_noapp = sim_no_app.add_flux(fcen, df, Nfreq, mp.FluxRegion(center=mp.Vector3(Lx/2-1.1*tpml,0), size = mp.Vector3(0,Ly), direction=mp.X))

sim_app.load_minus_flux('refl_bg',refl_app)
sim_no_app.load_minus_flux('refl_bg',refl_noapp)
#sim.load_minus_flux('trans_bg', trans)

sim_app.run(mp.at_beginning(mp.output_epsilon), #mp.at_every(1, mp.output_efield_z),
        until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt,src_r0, 1e-7))
sim_no_app.run(mp.at_beginning(mp.output_epsilon), #mp.at_every(1, mp.output_efield_z),
        until_after_sources=mp.stop_when_fields_decayed(50, src_cmpt,src_r0, 1e-7))


freqs = np.array(mp.get_flux_freqs(trans_bg))
r_bg_flux = np.array(mp.get_fluxes(refl_bg))
t_bg_flux = np.array(mp.get_fluxes(trans_bg))
r_flux = np.array(mp.get_fluxes(refl_app))
t_flux = np.array(mp.get_fluxes(trans_app))
r_flux_noapp = np.array(mp.get_fluxes(refl_noapp))
t_flux_noapp = np.array(mp.get_fluxes(trans_noapp))


T = t_flux/t_bg_flux
R = -1*r_flux/t_bg_flux
T_noapp = t_flux_noapp/t_bg_flux
R_noapp = -1*r_flux_noapp/t_bg_flux

if mp.am_master():
    
     plt.figure(figsize=(8, 8));
     plt.plot(1/freqs,R,label='R appodized')
     # plt.plot(1/freqs,R,label='R')
     plt.plot(1/freqs,R_noapp, label='R no appodization')
#     plt.plot(1/freqs,1-(R+T),label='Abs')
     plt.xlabel('Wavelength (um)',fontsize=20)
     plt.ylabel('R/T',fontsize=20)
     plt.legend()
     plt.show()




#sim.run(until=10)
