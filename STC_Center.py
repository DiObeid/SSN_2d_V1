from numpy import *
import sys
sys.path.append('path') # path to saved weights 
from STC_parameter import *
from scipy.special import erf     
from scipy import sparse
from functools import reduce


#*****************
# Extracting Map Information:

putmask(orient_map, orient_map > 179.9, 0.0) 
# Place neurons on the grid; each neuron has a location and a preferred orientation
x = arange(num_el) // ncol
y = arange(num_el) % ncol
theta = orient_map[x,y]

#*****************
# load the weight matrices:

Wee=load('path/Wee.npy')
Wie=load('path/Wie.npy')
Wei=load('path/Wei.npy')
Wii=load('path/Wii.npy')



#******************
# neuronal parameters                                                                                                                                                

taue=mtaue*ones(num_el)  
taui=mtaui*ones(num_el)  
aK=maK*ones(num_el)  
ne=mne*ones(num_el)  
ni=mni*ones(num_el)  


# SIMULATION LOOP
#****************************************************************************************************************************
 
def find_rate(alphacen): 
    
    
    # center
    def fc(x,y):
     profcenter= (erf(((lc/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((lc/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((lc/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((lc/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.) * profcenter
     return profile
  
    # annulus:                                                                                                                                                                                          
    # def fs(x,y):
    #  profout= (erf(((l/2.) + (x-xo))/(sqrt(2)*srf))+erf(((l/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((l/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((l/2.) - (y-yo))/(sqrt(2)*srf)))
    #  profin= (erf(((a/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((a/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((a/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((a/2.) - (y-yo))/(sqrt(2)*srf)))
    #  profile= (1./4.)*(profout-profin)
    #  return profile


    # Set the external inputs:
    diffangcen=absolute(alphacen-theta)
    diffangcen[diffangcen>90.]= 180.-diffangcen[diffangcen>90.]

    # diffangs=absolute(alphas-theta)
    # diffangs[diffangs>90.]= 180.-diffangs[diffangs>90.]

    he=NR(C1)*fc(x,y)*exp(-(diffangcen**2.)/(2*sig_fori**2.)) 
    hi=NR(C1)*fc(x,y)*exp(-(diffangcen**2.)/(2*sig_fori**2.)) 


    # Integrate the rate equations:
    re=zeros(num_el) 
    ri=zeros(num_el)
    zer=zeros(num_el)
    mre=zeros(num_el)
    mri=zeros(num_el)

    for j in range(T): 
       
       rece=dot(Wee,re)-dot(Wei,ri)
       reci=dot(Wie,re)-dot(Wii,ri)

       re= re + (dt/taue)* (-re +  aK* ( maximum(zer,he+rece) )**ne)
       ri= ri + (dt/taui)* (-ri +  aK* ( maximum(zer,hi+reci) )**ni)
       # if j>=Tss:      
       #  mre=mre+re
       #  mri=mri+ri

       unst=where(re>200.)
       size_unst=size(unst)
       if size_unst>0.:
        print "number of unstable cells", size_unst
        print "unstable cells", unst
        exit(1)       

 
    pe=re[recorded_neuron]
    pi=ri[recorded_neuron] 



     

    return(pe,pi)
   

#*************************************************************************************                                                                                                                                                                                   
if __name__ == "__main__":    

 # recorded neuron                                                                                                                                         
 S=loadtxt('sample_cells.txt')
 coordinate=int(sys.argv[1])                                                                                                                                       
 xo,yo= S[coordinate-1]
 recorded_neuron= 75*xo + yo
 recorded_neuron=int(recorded_neuron)

 # stimulus orientations                                                                                                                             
 pref=orient_map[int(xo),int(yo)]
 delta=50.
 rotate=union1d(arange(pref-delta,pref,10),arange(pref,pref+delta+10,10))
 rotate[rotate>=180.]=rotate[rotate>=180.]-180.
 rotate[rotate<0.]=180.+rotate[rotate<0.]
 
 
 

 # run experiment:
 alphacen=rotate[::-1]
 Center_tune=zeros((len(alphacen),3))   

 position=0
 for val in alphacen:
  response_e, response_i= find_rate(val)       
  Center_tune[position,0]=val
  Center_tune[position,1]=response_e
  Center_tune[position,2]=response_i
  position+=1   

 # save data
 savefitpath='/Center/'
 savefit='neuron_%.1f_pref-orientation_%.3f' % (recorded_neuron,theta[recorded_neuron])
 savefit.replace(' ','')

 fitname=savefitpath+ savefit +'.txt'
 savetxt(fitname,Center_tune)
    
 
