from numpy import *
import sys
sys.path.append('path') # path to saved weights
from FSS_parameter import *
from scipy.special import erf     
from scipy import sparse
from scipy.optimize import curve_fit
from functools import reduce
import itertools as it



#*****************
# Extracting Map Information:

putmask(orient_map, orient_map > 179.9, 0.0)  
# Place neurons on the grid; each neuron has a location and a preferred orientation
x = arange(num_el) // ncol
y = arange(num_el) % ncol
theta = orient_map[x,y]

#*****************
# load weight matrices:

Wee=load('/path/Wee.npy')
Wie=load('/path/Wie.npy')
Wei=load('path/Wei.npy')
Wii=load('path/Wii.npy')



#*****************
# neuronal parameters                                                                                                                                                

taue=mtaue*ones(num_el) 
taui=mtaui*ones(num_el) 
aK=maK*ones(num_el)  
ne=mne*ones(num_el)  
ni=mni*ones(num_el)  


# SIMULATION LOOP
#****************************************************************************************************************************
 
def find_rate(alphacen1,alphacen2,alphas,s1,s2,s3): 
    
    
    # center
    def fc(x,y):
     profcenter= (erf(((lc/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((lc/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((lc/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((lc/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.) * profcenter
     return profile
  

    # annulus:                                                                                                                                                                                          
    def fs(x,y):
     profout= (erf(((l/2.) + (x-xo))/(sqrt(2)*srf))+erf(((l/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((l/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((l/2.) - (y-yo))/(sqrt(2)*srf)))
     profin= (erf(((a/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((a/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((a/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((a/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.)*(profout-profin)
     return profile



    # Set the external inputs:
    diffangcen1=absolute(alphacen1-theta)
    diffangcen1[diffangcen1>90.]= 180.-diffangcen1[diffangcen1>90.]
 
    diffangcen2=absolute(alphacen2-theta)
    diffangcen2[diffangcen2>90.]= 180.-diffangcen2[diffangcen2>90.]
    
    diffangs=absolute(alphas-theta)
    diffangs[diffangs>90.]= 180.-diffangs[diffangs>90.]

    he= s1*NR(C1)*fc(x,y)*exp(-(diffangcen1**2.)/(2*sig_fori**2.)) + s2*NR(C1)*fc(x,y)*exp(-(diffangcen2**2.)/(2*sig_fori**2.)) + s3*NR(C2)*fs(x,y)* exp(-(diffangs**2.)/(2*sig_fori**2.))
    hi= s1*NR(C1)*fc(x,y)*exp(-(diffangcen1**2.)/(2*sig_fori**2.)) + s2*NR(C1)*fc(x,y)*exp(-(diffangcen2**2.)/(2*sig_fori**2.)) + s3*NR(C2)*fs(x,y)* exp(-(diffangs**2.)/(2*sig_fori**2.))

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
       if size_unst>0:
        print "number of unstable cells", size_unst
        print "unstable cells", unst
        print "re_unstable", re[unst], "ri_unstable", ri[unst]
        exit(1)       


    # Analysis
    # mre=(1./(T-Tss))*mre 
    # mri=(1./(T-Tss))*mri

    # center population
    re_pop=re[pop]
    ri_pop=ri[pop]
    theta_pop=theta[pop]

    # # if you want the row population activity vector
    # orient_av_rate=zeros((len(re[pop]),3))
    # orient_av_rate[:,0]=theta_pop
    # orient_av_rate[:,1]=re_pop
    # orient_av_rate[:,2]=ri_pop
    
    # binned population activity by orientation: average rate of neurons with similar orientations
    numbins=36
    bin=linspace(0,180,numbins)
    B=digitize(theta_pop,bin)
    orient_av_rate=zeros((numbins,3))
    for b in range(numbins):
     tmp=((B-1)==b)
     if any(tmp):
      orient_av_rate[b,0]=mean(theta_pop[tmp])
      orient_av_rate[b,1]=mean(re_pop[tmp])
      orient_av_rate[b,2]=mean(ri_pop[tmp])


      
    return(orient_av_rate)
  
 
#*************************************************************************************                                                                                                                                              
## non-zero angle plaid

def Run(ang,alphacen2,alphas):
 try:
  # conditions:
  stim1=[1,0,1,1]
  stim2=[0,1,1,1]
  stim3=[0,0,0,1]

  for s1,s2,s3 in zip(stim1,stim2,stim3):
     av_pop_rate=find_rate(ang,alphacen2,alphas,s1,s2,s3)
     xdata=av_pop_rate[:,0] 
     ydata=av_pop_rate[:,1]

     if sum(ydata)<30.:
      break

     else:
     # Statistical fitting: 
      def VonMises(x,a,b,mu):
         f=a*exp(b*cos(2*0.0175*x-mu))
         return f
      def double_VonMises(x,w1,w2,a1,b1,mu1,a2,b2,mu2):
        f= w1*a1*exp(b1*cos(2*0.0175*x-mu1)) + w2*a2*exp(b2*cos(2*0.0175*x-mu2))
        return f
                      
      if [s1,s2,s3]==[1,0,0]:
       po=[0.5,0.5,0.5]
       inv_weights=sqrt(mean(ydata)/(ydata+1e-10))  # sigma in python = sqrt(1/w) where w is the weight vector in matlab fit
       popt,popcov=curve_fit(VonMises,xdata,ydata,po,sigma=inv_weights)
       a1=popt[0]
       b1=popt[1]
       mu1=popt[2]

      if [s1,s2,s3]==[0,1,0]:
       po=[0.5,0.5,0.5]
       inv_weights=sqrt(mean(ydata)/(ydata+1e-10))
       popt,popcov=curve_fit(VonMises,xdata,ydata,po,sigma=inv_weights)
       a2=popt[0]
       b2=popt[1]
       mu2=popt[2]
                             
      if [s1,s2,s3]==[1,1,0]:
       von2= lambda x,w1,w2 : double_VonMises(x,w1,w2,a1,b1,mu1,a2,b2,mu2)
       pod=[0.5,0.5]
       inv_weights=sqrt(mean(ydata)/(ydata+1e-10))
       poptd,popcovd=curve_fit(von2,xdata,ydata,pod,sigma=inv_weights,bounds=(0,2))
       w1_p=poptd[0]
       w2_p=poptd[1]
       W_P[0].append(w1_p)
       W_P[1].append(w2_p)
     
      if [s1,s2,s3]==[1,1,1]:
       von2= lambda x,w1,w2 : double_VonMises(x,w1,w2,a1,b1,mu1,a2,b2,mu2)
       pod=[0.5,0.5]
       inv_weights=sqrt(mean(ydata)/(ydata+1e-10))
       poptd,popcovd=curve_fit(von2,xdata,ydata,pod,sigma=inv_weights,bounds=(0,2))
       w1_ps=poptd[0]
       w2_ps=poptd[1]
       W_PS[0].append(w1_ps)
       W_PS[1].append(w2_ps) 
 except RuntimeError:
  print ' Run_time_Error','alphacen1=',ang,'alphacen2',alphacen2
 return
 #***********************************************


if __name__ == "__main__":    

 # recorded neuron                                                                                                                                         
 S=loadtxt('sample_cells.txt')
 coordinate=int(sys.argv[1])                                                                                                                                       
 xo,yo= S[coordinate-1]
 recorded_neuron= 75*xo + yo
 recorded_neuron=int(recorded_neuron)


 # define local population                                                                                                                                                                                  
 dpop=4.
 pop=[]
 for posx,posy in it.product(arange(xo-rint(dpop/2),xo+1+rint(dpop/2)),arange(yo-rint(dpop/2),yo+1+rint(dpop/2))):
  pop.append(75*posx+posy)
 pop=array(pop).astype(int)

 # center stimuli
 plaid_angle=float(sys.argv[2])
 # alphacen1=arange(0.,180.-plaid_angle,10.) # positive plaid angles
 alphacen1=arange(0.-plaid_angle,180.,10.) # negative plaid angles 

 # create lists for saving data
 W_P=[[] for i in range(2)]
 W_PS=[[] for i in range(2)]
 
 for ang in alphacen1:
  alphacen2= ang + plaid_angle
  alphas=alphacen2
  Run(ang,alphacen2,alphas)

 
 savepar='plaid_angle_%.1f_neuron_%.1f' % (plaid_angle,recorded_neuron)
 savepar.replace(' ','')
 

 filename1='/Plaid/Plaid_Weights_'+ savepar +'.txt'
 savetxt(filename1,W_P)

 filename2='/Plaid_Surround/Plaid+Surround_Weights_'+ savepar +'.txt'             
 savetxt(filename2,W_PS)     



