from numpy import *
from scipy.special import erf

########################################################################


orient_map=loadtxt('/orient_map.txt',delimiter=',') # import the orientation map as an array                            

num_el=orient_map.size 

ncol=orient_map.shape[1]  

delta_theta=16./(ncol-1) 
delta_x= delta_theta*0.6  


n_loc=2.                    
n_loc_e=1.                  
N_pop= n_loc*num_el        
Ne= n_loc_e*num_el         
Ni= (n_loc-n_loc_e)*num_el  

# neuron parameters:
mtaue= 10. 
mtaui=6.67 
mne=2.2 
mni=2.2
maK=0.01 
srf=0.25/0.6  
        


# connection parameters:
Rc=3.
sig_spatial_Far_ie =6. 
sig_spatial_Far_ee =3.  
sig_spatial_i=2.   
sig_orient_Near_e=55.     
sig_orient_Far_e=25.     
sig_orient_i=55.         



# connection weights 
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

sig_fori=20.                                                                                                                                     
                                                                                                                                                                                                            

# stimulus parameters
C1=16.4 



# simulation parameters:

run_time=500. 
dt=0.5
T=int(float(run_time/dt)) 


