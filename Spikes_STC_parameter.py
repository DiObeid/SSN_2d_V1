from brian import *
from numpy import *
from scipy.special import erf

########################################################################

# number of nerons in the network:

orient_map=loadtxt('orient_map8c.txt',delimiter=',')

num_el=orient_map.size 


ncol=orient_map.shape[1] 

delta_theta=16./(ncol-1) 

delta_x= delta_theta*0.6 


n_loc=2                    
n_loc_e=1                  
N_pop= n_loc*num_el        
Ne= n_loc_e*num_el         
Ni= (n_loc-n_loc_e)*num_el 

# connections
Rc=3.
sig_spatial_Far_ie=6.
sig_spatial_Far_ee=3.
sig_spatial_i=2.
sig_orient_Near_e=55.
sig_orient_Far_e=25.
sig_orient_i=55.

# neuron parameters:
srf=0.25/0.6
taue=3*ms
taui=3*ms 
taumE=15*ms
taumI=15*ms
tau_ref=3*ms
Vt=-50*mV
Vr=-56*mV
El=-70*mV
Ee=0*mV
Ei=-80*mV
gl=10.*nS

psi=1.
# conductaces 
gee_near=psi*1.8*nS     
gie_near=psi*1.76*nS 

gee_far=psi*0.7*nS 
gie_far=psi*0.65*nS   

gei=psi*3.3*nS   
gii=psi*2.*nS 


#  input :

def NR(C):
 Rmax=50.*Hz
 C50=11.
 response=(Rmax* C**3.5)/( C50**3.5 + C**3.5)
 return response
 
gext=0.1*nS 
sig_fori=20.      
Ninput=200.   
sigma_bg=2.5*mV

# stimuli parameters:                                                                                                                                                                                                               
C1=100.
C2=100.
lc=6.
a=20. 
l=100.

# stimulus at which we monitor network behavior                                                                                                                                                                                     
l_monitor_network=25.

# run
run_time=8.*second 

