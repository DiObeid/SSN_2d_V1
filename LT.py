from numpy import *
import sys
sys.path.append('path') # path to saved weights
from LT_parameter import *
from scipy.special import erf     
from scipy import sparse
from functools import reduce



#**********************
# Extracting Map Information:

putmask(orient_map, orient_map > 179.9, 0.0) 
# Place neurons on the grid; each neuron has a location and a preferred orientation
x = arange(num_el) // ncol
y = arange(num_el) % ncol
theta = orient_map[x,y]

#***********************
# Load  weight matrices:

Wee=load('path/Wee.npy')
Wie=load('path/Wie.npy')
Wei=load('path/Wei.npy')
Wii=load('path/Wii.npy')


#*********************
# neuronal parameters                                                                                                                                                

taue=mtaue*ones(num_el)  
taui=mtaui*ones(num_el)  
aK=maK*ones(num_el)  
ne=mne*ones(num_el)  
ni=mni*ones(num_el)  


# SIMULATION LOOP
#****************************************************************************************************************************
 
def find_rate(lc): 
    
    
    # center
    def f1(x,y):
     profcenter= (erf(((lc/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((lc/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((lc/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((lc/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.) * profcenter
     return profile

    # annulus:                                                                                                                                                                                          
    # def f2(x,y):
    #  profout= (erf(((l/2.) + (x-xo))/(sqrt(2)*srf))+erf(((l/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((l/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((l/2.) - (y-yo))/(sqrt(2)*srf)))
    #  profin= (erf(((a/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((a/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((a/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((a/2.) - (y-yo))/(sqrt(2)*srf)))
    #  profile= (1./4.)*(profout-profin)
    #  return profile


    # Set the external inputs:
    diffangs1=absolute(alpha1-theta)
    diffangs1[diffangs1>90.]= 180.-diffangs1[diffangs1>90.]

    # diffangs2=absolute(alpha2-theta)
    # diffangs2[diffangs2>90.]= 180.-diffangs2[diffangs2>90.]

    he= NR(C1)*f1(x,y)*exp(-(diffangs1**2.)/(2*sig_fori**2.)) 
    hi= NR(C1)*f1(x,y)*exp(-(diffangs1**2.)/(2*sig_fori**2.)) 


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


    # Analysis:
    # mre=(1./(T-Tss))*mre 
    # mri=(1./(T-Tss))*mri


    #inputs currents to the two populations
    exc_input_e= he + dot(Wee,re)
    inh_input_e= dot(Wei,ri)
 
    exc_input_i= hi + dot(Wie,re)
    inh_input_i= dot(Wii,ri)
  
    length.append(lc)
    save_e.append(exc_input_e[recorded_neuron].item())
    save_i.append(inh_input_e[recorded_neuron].item())

    savi_e.append(exc_input_i[recorded_neuron].item())
    savi_i.append(inh_input_i[recorded_neuron].item())
 

    return(pe,pi)
   

 
#*************************************************************************************                                                                                                                                                                                   
if __name__ == "__main__":    

 # recorded neuron                                                                                                                                         
 S=loadtxt('sample_cells.txt')
 coordinate=int(sys.argv[1])                                                                                                                                       
 xo,yo= S[coordinate-1]
 recorded_neuron= 75*xo + yo
 recorded_neuron=int(recorded_neuron)

 # stimulus orientation                                                                                                                             
 alpha1=orient_map[int(xo),int(yo)]

 # lists to save currents:
 length=[]
 save_e=[]
 save_i=[]
 savi_e=[]
 savi_i=[]
 
 # run the exp. for various  stimulus sizes:
 lc=reduce(union1d,(arange(1.,12.,1.),arange(12.,20.,2.),arange(20.,80.,4.)))
 L_tune=zeros((len(lc),3))    

 pos=0
 for val in lc:
  response_e, response_i= find_rate(val)       
  L_tune[pos,0]=val
  L_tune[pos,1]=response_e
  L_tune[pos,2]=response_i
  pos+=1   

 
 savefit='C1_%.1f_neuron_%.1f_pref-orientation_%.3f' % (C1,recorded_neuron,theta[recorded_neuron])
 savefit.replace(' ','')

 fitname='/Length_Tuning/length_tuning_'+ savefit +'.txt'
 savetxt(fitname,L_tune)
    
 save=[length,save_e,save_i]
 asave=transpose(array(save))

 savi=[length,savi_e,savi_i]
 asavi=transpose(array(savi))


 fitname1='/LT_Inputs/E_cells/inputs_to_Ecell_'+ savefit +'.txt'
 savetxt(fitname1,asave)

 fitname2='/LT_Inputs/I_cells/inputs_to_Icell_'+ savefit +'.txt'
 savetxt(fitname2,asavi)
