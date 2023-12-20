from numpy import *
from scipy.special import erf

########################################################################

# number of nerons in the network:

orient_map=loadtxt('/orient_map.txt',delimiter=',') # import the orientation map as an array                            

num_el=orient_map.size  # number of elements in the array 

# Assume matrix corresponds to 16deg x 16deg cortical area and use cortical magnification of 0.6mm/deg

ncol=orient_map.shape[1]  # number of columns

delta_theta=16./(ncol-1)  #distance between two points in the matrix in degrees

delta_x= delta_theta*0.6  # distance between two points in the matrix in mm


n_loc=2.                     # number of neurons at each location of the array or map
n_loc_e=1.                   # number of E neurons at each location 
N_pop= n_loc*num_el         # total number of neurons
Ne= n_loc_e*num_el          # total number of E (excitatory)  neurons
Ni= (n_loc-n_loc_e)*num_el  # total number of I (inhibitory) cells


# neuron parameters:
mtaue= 10. 
mtaui=6.67 
mne=2.2 
mni=2.2
maK=0.01 
srf=0.25/0.6  
        


# connection parameters:
Rc=3. # the critical distance determining the separation between "Near" and "Far"

sig_spatial_Far_ie =6.  # standard deviation of the spatial connectivity of E to I neurons   
sig_spatial_Far_ee =3.  # standard deviation of the spatial connectivity of E to E neurons
sig_spatial_i=2.   # standard deviation of the spatial connectivity for I neurons 
sig_orient_Near_e=55.     
sig_orient_Far_e=25.     # standard deviation of the orientation dependent part of the connectivity of E cells
sig_orient_i=55.         # standard deviation of the orientation dependent part of the connectivity of the I cells



# Connection weights 

psi=1.2

#excitatory
Jee_Near=psi*0.06        
Jie_Near=psi*0.05        

Jee_Far=psi*0.03 
Jie_Far=psi*0.03 

#inhibitory
Jei=psi*0.044   
Jii=psi*0.024   


# input:                                                                                                                                                                                                    
def NR(C): 
 Rmax=50.
 C50=11.
 response=(Rmax* C**3.5)/( C50**3.5 + C**3.5)
 return response

sig_fori=20. # standard deviation in orientation space of the input                                                                                                                                       
                                                                                                                                                                                                            

# stimulus parameters
C1=16.4 # center stim  contrast 



# Simulation parameters:


run_time=500. 
dt=0.5
T=int(float(run_time/dt)) 


