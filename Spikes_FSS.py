from numpy import *
import sys
from brian import *
sys.path.append('path')   # path to saved weights
from Spikes_FSS_parameter import *
from functools import reduce
import itertools as it
from scipy.special import erf     
from scipy.optimize import curve_fit
import gc



#**********
# Extracting Map Information:                                                          
putmask(orient_map, orient_map > 179.9, 0.0)
x = arange(num_el) // ncol
y = arange(num_el) % ncol
theta = orient_map[x,y]

#**********
# load weight matrices                                                                                                                                                                                      
Wee=load('path/Wee.npy')
Wie=load('path/Wie.npy')
Wei=load('path/Wei.npy')
Wii=load('path/Wii.npy')


# SIMULATION LOOP
#****************************************************************************************************************************

def find_rate(alphacen1,alphacen2,alphas,s1,s2,s3):

    clear(True,True)
    gc.collect()
    reinit_default_clock(t=0.0*second)
    defaultclock.dt=0.1*ms
    #mon_clock=Clock(dt=2.0*ms)
    random.seed()


    def fc(x,y):
     profcenter= (erf(((lc/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((lc/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((lc/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((lc/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.)*profcenter
     return profile
 
    def fs(x,y):
      profout= (erf(((l/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((l/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((l/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((l/2.) - (y-yo))/(sqrt(2)*srf)))
      profin= (erf(((a/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((a/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((a/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((a/2.) - (y-yo))/(sqrt(2)*srf)))
      profile= (1./4.)*(profout-profin)
      return profile   


    eqs= Equations('''
               dV/dt= (-(V-El) + (ge/gl)*(Ee-V) + (gi/gl)*(Ei-V) + (gin/gl)*(Ee-V))/taum + sigma_bg*xi/(ms**0.5)  : volt   
               taum: second
               dge/dt=-ge/taue            : siemens
               dq/dt=-q/taue + (1/sigmag)*(mean_gin/taue)+ xi/taue    : second**(-0.5)
               gin= sigmag*q               : siemens
               mean_gin                    : siemens
               sigmag                      : siemens*second**(0.5)
               dgi/dt=-gi/taui            : siemens
               Ie_N=ge*(Ee-V)             : amp
               Ii_N=gi*(Ei-V)             : amp
               Ie_in=gin*(Ee-V)           : amp
               ''')

    G= NeuronGroup(N=N_pop, model=eqs, threshold=Vt , reset=Vr, refractory=tau_ref)

    # initialize the value of V for each neuron in the group:
    G.V= El+ (Vt-Vr)*rand(len(G))
    G.ge=zeros(len(G))
    G.gi=zeros(len(G))
    G.q=zeros(len(G))
    
    # Divide the population into Ne excitatory neurons and Ni inhibitory neurons:
    Ge=G.subgroup(Ne)
    Gi=G.subgroup(Ni)
    Ge.taum=taumE*ones(len(Ge))
    Gi.taum=taumI*ones(len(Gi))
   
    # Set the Connections and the Connectivity:
    Cee=Connection(Ge,Ge,'ge')
    Cie=Connection(Ge,Gi,'ge')
    Cii=Connection(Gi,Gi,'gi')
    Cei=Connection(Gi,Ge,'gi')

    Cee.connect(Ge,Ge,Wee)
    Cie.connect(Ge,Gi,Wie)
    Cii.connect(Gi,Gi,Wii)
    Cei.connect(Gi,Ge,Wei)

    # Input:

    diffangcen1=absolute(alphacen1-theta)
    diffangcen1[diffangcen1>90.]= 180. - diffangcen1[diffangcen1>90.]

    diffangcen2=absolute(alphacen2-theta)
    diffangcen2[diffangcen2>90.]= 180. - diffangcen2[diffangcen2>90.]

    diffangs=absolute(alphas-theta)
    diffangs[diffangs>90.]= 180. - diffangs[diffangs>90.]

    rate=s1*NR(C1)*exp(-(diffangcen1**2.)/(2*sig_fori**2.))*fc(x,y) + s2*NR(C1)*exp(-(diffangcen2**2.)/(2*sig_fori**2.))*fc(x,y) + s3*NR(C2)*exp(-(diffangs**2.)/(2*sig_fori**2.))*fs(x,y) + 1e-8*Hz
    
    Mean=Ninput*gext*taue*rate

    SD=(Ninput*rate*(taue**2)*(gext**2))**0.5

    Ge.mean_gin= Mean
    Ge.sigmag= SD
    Gi.mean_gin= Mean
    Gi.sigmag= SD
    
    #  Monitors: 
    Se = SpikeMonitor(Ge)
    Si = SpikeMonitor(Gi)
    # E_V=StateMonitor(Ge,'V', record=recorded_neuron)
    # E_Ie=StateMonitor(Ge,'Ie_N',record= samecells)
    # E_Ii=StateMonitor(Ge,'Ii_N',record= samecells)
    # I_Ie=StateMonitor(Gi,'Ie_N',record= samecells)
    # I_Ii=StateMonitor(Gi,'Ii_N',record=samecells)
    # E_Iin=StateMonitor(Ge,'Ie_in',record=samecells)
    # I_Iin=StateMonitor(Gi,'Ie_in',record=samecells)
    # E_ge=StateMonitor(Ge,'ge',record= samecells)
    # E_gi=StateMonitor(Ge,'gi',record= samecells)
    # I_ge=StateMonitor(Gi,'ge',record= samecells)
    # I_gi=StateMonitor(Gi,'gi',record= samecells)
    
    # Run the simulation :
    run(run_time)


    #Statistics :

    #1- rates: 
    rate_e=zeros(len(samecells))
    rate_i=zeros(len(samecells))
    index1=0
    for i in samecells:
      if len(Se[i]) < 1:
         rate_e[index1]=0.0
      else :
         rate_e[index1]= len(Se[i])/run_time
      index1+=1
      
    index11=0  
    for i in samecells:
      if len(Si[i])< 1:
         rate_i[index11]=0.0
      else:    
         rate_i[index11]=len(Si[i])/run_time
      index11+=1
      
    # binned rate data: average over neurons within an orientation bin   
    numbins=36
    bin=linspace(0,180,numbins)
    theta_pop=theta[samecells]
    B=digitize(theta_pop,bin)
    orient_avg_rate=zeros((numbins,3))
    for b in range(numbins):
     tmp=((B-1)==b)
     if any(tmp):
      orient_avg_rate[b,0]=mean(theta_pop[tmp])
      orient_avg_rate[b,1]=mean(rate_e[tmp])
      orient_avg_rate[b,2]=mean(rate_i[tmp])

 
    # #2- currents:

    # Mean_E_Ie=mean(E_Ie.values,axis=1)   
    # Mean_E_Ii=mean(E_Ii.values,axis=1)
    # Mean_I_Ie=mean(I_Ie.values,axis=1)
    # Mean_I_Ii=mean(I_Ii.values,axis=1)
    # Mean_E_Iin=mean(E_Iin.values,axis=1)
    # Mean_I_Iin=mean(I_Iin.values,axis=1)
    # Mean_E_ge=mean(E_ge.values,axis=1)
    # Mean_E_gi=mean(E_gi.values,axis=1)
    # Mean_I_ge=mean(I_ge.values,axis=1)
    # Mean_I_gi=mean(I_gi.values,axis=1)

    # Avg_E_Ie=mean(Mean_E_Ie)
    # Avg_E_Ii=mean(Mean_E_Ii)
    # Avg_I_Ie=mean(Mean_I_Ie)
    # Avg_I_Ii=mean(Mean_I_Ii)
    # Avg_E_Iin=mean(Mean_E_Iin)
    # Avg_I_Iin=mean(Mean_I_Iin)
    # Avg_E_ge=mean(Mean_E_ge)
    # Avg_E_gi=mean(Mean_E_gi)
    # Avg_I_ge=mean(Mean_I_ge)
    # Avg_I_gi=mean(Mean_I_gi)
    
    # #3- cv's 
    # if lc==l_monitor_network:
    #   index3=0
    #   cv_e=zeros(len(samecells))
    #   cv_i=zeros(len(samecells))
    #   for i in samecells:
    #     cv_e[index3]=(3. if len(Se[i])<=1 else CV(Se[i]))
    #     cv_i[index3]=(3. if len(Si[i])<=1 else CV(Si[i]))
    #     index3+=1
    #   cv_e=nan_to_num(cv_e)
    #   cv_i=nan_to_num(cv_i)



    # #append currents:
    # Sur_angle.append(alphas)
    # save_e.append(Avg_E_Ie)
    # save_i.append(Avg_E_Ii)

    # savi_e.append(Avg_I_Ie)
    # savi_i.append(Avg_I_Ii)    

    # save_in.append(Avg_E_Iin)
    # savi_in.append(Avg_I_Iin)

    # save_ge.append(Avg_E_ge)
    # save_gi.append(Avg_E_gi)

    # savi_ge.append(Avg_I_ge)
    # savi_gi.append(Avg_I_gi)
    
    return(orient_avg_rate) 

#*************************************************************************************                                                                                                                                              

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
  print 'Runtime_Error', 'alphacen1=',ang,'alphacen2=',alphacen2
 return
#************************************************************************************************                                     

if __name__ == "__main__":
 #recorded neuron
 S=loadtxt('sample_cells.txt')
 coordinate=int(sys.argv[1]) 
 xo,yo=S[coordinate-1]
 xo=int(xo)
 yo=int(yo)
 recorded_neuron=75*xo+yo
 pref=orient_map[xo,yo]

 # define local population:
 dpop=4.
 pop=[]
 for posx,posy in it.product(arange(xo-rint(dpop/2.),xo+1+rint(dpop/2.)),arange(yo-rint(dpop/2.),yo+1+rint(dpop/2.))):
  pop.append(75*posx+posy)
 pop=array(pop).astype(int)
 samecells=pop.tolist()

 # center stimulus  
 plaid_angle=int(sys.argv[2])
 #alphacen1=arange(0.,180.-plaid_angle,10.) 
 alphacen1=arange(0.-plaid_angle,180.,10.)   
 
 # create lists for saving data                                                                                                                     
 W_P=[[] for i in range(2)]
 W_PS=[[] for i in range(2)]

 for ang in alphacen1:
  alphacen2= ang + plaid_angle
  alphas=alphacen2
  Run(ang,alphacen2,alphas)

 # save weights:
 savepar='plaid_angle_%.1f_neuron_%.1f' % (plaid_angle,recorded_neuron)
 savepar.replace(' ','')

 filename1='/path/Plaid/Plaid_Weights_'+ savepar +'.txt'
 savetxt(filename1,W_P)

 filename2='/path/Plaid_Surround/Plaid+Surround_Weights_'+ savepar +'.txt'
 savetxt(filename2,W_PS)

  
