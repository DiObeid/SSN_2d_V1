from numpy import *
from brian import *
import sys
sys.path.append('path')  # path to saved weights
from Spikes_LT_parameter import *
from functools import reduce
import itertools as it
from scipy.special import erf     
from scipy import stats
import gc
import time
   

#*****************
#Extracting Map information:

putmask(orient_map, orient_map > 179.9, 0.0)  
x = arange(num_el) // ncol
y = arange(num_el) % ncol
theta = orient_map[x,y]

#*****************
# load weight matrices                                                                                                                                      
Wee=load('path/Wee.npy')
Wie=load('path/Wie.npy')
Wei=load('path/Wei.npy')
Wii=load('path/Wii.npy')



# SIMULATION LOOP
#****************************************************************************************************************************

def find_rate(lc):

    clear(True,True)
    gc.collect()
    reinit_default_clock(t=0.0*second)
    defaultclock.dt=0.1*ms
    #mon_clock=Clock(dt=2.0*ms)
    random.seed()    

    
    def f1(x,y):
     profcenter= (erf(((lc/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((lc/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((lc/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((lc/2.) - (y-yo))/(sqrt(2)*srf)))
     profile= (1./4.)*profcenter
     return profile
  

    # def f2(x,y):
    #   profout= (erf(((l/2.) + (x-xo))/(sqrt(2)*srf))+erf(((l/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((l/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((l/2.) - (y-yo))/(sqrt(2)*srf)))
    #   profin= (erf(((a/2.) + (x-xo))/(sqrt(2)*srf))+ erf(((a/2.) - (x-xo))/(sqrt(2)*srf))) * (erf(((a/2.) + (y-yo))/(sqrt(2)*srf))+ erf(((a/2.) - (y-yo))/(sqrt(2)*srf)))
    #   profile= (1./4.)*(profout-profin)
    #   return profile   



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

    # initialize the value of the parameters for each neuron in the group:
    G.V= El+(Vt-Vr)*rand(len(G))
    G.ge=zeros(len(G))
    G.gi=zeros(len(G))
    G.q=zeros(len(G))
    
    # Divide the population into Ne excitatory neurons and Ni inhibitory neurons:
    Ge=G.subgroup(Ne)
    Gi=G.subgroup(Ni)
    Ge.taum=taumE*ones(len(Ge))
    Gi.taum=taumI*ones(len(Gi))
    subgpe=Ge[recorded_neuron]
    subgpi=Gi[recorded_neuron]
   
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

    diffangs1=absolute(alpha1-theta)
    diffangs1[diffangs1>90.]= 180. - diffangs1[diffangs1>90.]

    #diffangs2=absolute(alpha2-theta)
    #diffangs2[diffangs2>90.]= 180. - diffangs2[diffangs2>90.]

    rate=NR(C1)*exp(-(diffangs1**2.)/(2*sig_fori**2.))*f1(x,y) + 1e-6*Hz
    
    Mean=Ninput*gext*taue*rate

    SD=(Ninput*rate*(taue**2)*(gext**2))**0.5

    Ge.mean_gin= Mean
    Ge.sigmag= SD
    Gi.mean_gin= Mean
    Gi.sigmag= SD
    
    #  Monitors: 
    Se = SpikeMonitor(Ge)
    Si = SpikeMonitor(Gi)
    E_Ie=StateMonitor(Ge,'Ie_N',record= samecells)
    E_Ii=StateMonitor(Ge,'Ii_N',record= samecells)
    I_Ie=StateMonitor(Gi,'Ie_N',record= samecells)
    I_Ii=StateMonitor(Gi,'Ii_N',record=samecells)
    E_Iin=StateMonitor(Ge,'Ie_in',record=samecells)
    I_Iin=StateMonitor(Gi,'Ie_in',record=samecells)
    #E_V=StateMonitor(Ge,'V',record=samecells)
    E_ge=StateMonitor(Ge,'ge',record= samecells)
    E_gi=StateMonitor(Ge,'gi',record= samecells)
    I_ge=StateMonitor(Gi,'ge',record= samecells)
    I_gi=StateMonitor(Gi,'gi',record= samecells)
    Rte=PopulationRateMonitor(subgpe,bin=500*ms)
    Rti=PopulationRateMonitor(subgpi,bin=500*ms)

    
    # Run the simulation :
    run(run_time)



    #Statistics :

    # 1- rates: 
    # rate_e=zeros(len(samecells))
    # rate_i=zeros(len(samecells))
    # index1=0
    # for i in samecells:
    #   if len(Se[i]) < 1:
    #      rate_e[index1]=0.0
    #   else :
    #      rate_e[index1]= len(Se[i])/run_time
    #   index1+=1
      
    # index11=0  
    # for i in samecells:
    #   if len(Si[i])< 1:
    #      rate_i[index11]=0.0
    #   else:    
    #      rate_i[index11]=len(Si[i])/run_time
    #   index11+=1
      
    # Avg_rate_e=mean(rate_e)  
    # Avg_rate_i=mean(rate_i)

    #method 2 for computing the rates with sem  
    mRte=mean(Rte.rate[:-1])
    mRti=mean(Rti.rate[:-1]) 
    #std_Rte=std(Rte.rate,ddof=1)
    #std_Rti=std(Rti.rate,ddof=1)
    sem_Rte=stats.sem(Rte.rate[:-1])
    sem_Rti=stats.sem(Rti.rate[:-1])

    #2- currents:
    
    Mean_E_Ie=mean(E_Ie.values,axis=1)   
    Mean_E_Ii=mean(E_Ii.values,axis=1)
    Mean_I_Ie=mean(I_Ie.values,axis=1)
    Mean_I_Ii=mean(I_Ii.values,axis=1)
    Mean_E_Iin=mean(E_Iin.values,axis=1)
    Mean_I_Iin=mean(I_Iin.values,axis=1)
    Mean_E_ge=mean(E_ge.values,axis=1)
    Mean_E_gi=mean(E_gi.values,axis=1)
    Mean_I_ge=mean(I_ge.values,axis=1)
    Mean_I_gi=mean(I_gi.values,axis=1)
    
    Avg_E_Ie=mean(Mean_E_Ie)
    Avg_E_Ii=mean(Mean_E_Ii)
    Avg_I_Ie=mean(Mean_I_Ie)
    Avg_I_Ii=mean(Mean_I_Ii)
    Avg_E_Iin=mean(Mean_E_Iin)
    Avg_I_Iin=mean(Mean_I_Iin)
    Avg_E_ge=mean(Mean_E_ge)
    Avg_E_gi=mean(Mean_E_gi)
    Avg_I_ge=mean(Mean_I_ge)
    Avg_I_gi=mean(Mean_I_gi)
    

    #3- cv's 
    if lc in l_monitor_network:           
      index3=0
      cv_e=zeros(len(samecells))
      cv_i=zeros(len(samecells))
      for i in samecells:
        cv_e[index3]=(3. if len(Se[i])<=1 else CV(Se[i]))
        cv_i[index3]=(3. if len(Si[i])<=1 else CV(Si[i]))
        index3+=1
      cv_e=nan_to_num(cv_e)
      cv_i=nan_to_num(cv_i)
      print 'lc=', lc
      print 'cv_e=', cv_e
      print 'cv_i=', cv_i
      
        
    #append currents:
    length.append(lc)
    save_e.append(Avg_E_Ie)
    save_i.append(Avg_E_Ii)

    savi_e.append(Avg_I_Ie)
    savi_i.append(Avg_I_Ii)    

    save_in.append(Avg_E_Iin)
    savi_in.append(Avg_I_Iin)

    save_ge.append(Avg_E_ge)
    save_gi.append(Avg_E_gi)

    savi_ge.append(Avg_I_ge)
    savi_gi.append(Avg_I_gi)
    
    sav_gin.append(Mean[recorded_neuron])
    
    return(mRte,mRti,sem_Rte,sem_Rti) 

#*************************************************************************************                                                                                                                                                                                   

if __name__ == "__main__":
 #recorded neuron
 S=loadtxt('sample_cells.txt')
 coordinate=int(sys.argv[1])
 xo,yo=S[coordinate-1]
 recorded_neuron=75*xo+yo
 recorded_neuron=int(recorded_neuron)
 alpha1=orient_map[int(xo),int(yo)]
 
 #*******
 # define local population:
 #dpop=2.
 #pop=[]
 #for posx,posy in it.product(arange(xo-rint(dpop/2.),xo+1+rint(dpop/2.)),arange(yo-rint(dpop/2.),yo+1+rint(dpop/2.))):
 # pop.append(75*posx+posy)
 #samecells=pop
 
 # # neurons of similar orientation
 # logic=where(logical_and(theta>=(theta[recorded_neuron]-0.5),theta<=(theta[recorded_neuron]+0.5))) 
 # samecells=(logic[0]).astype(int)
 # samecells=samecells.tolist()
 #*******

 samecells=[recorded_neuron]  

 # lists to save currents and conductances
 length=[]
 save_e=[]
 save_i=[]
 savi_e=[]
 savi_i=[]
 save_in=[]
 savi_in=[]
 save_ge=[]
 save_gi=[]
 savi_ge=[]
 savi_gi=[]
 sav_gin=[]

 # run the exp. for various stimulus sizes: 
 lc=reduce(union1d,(arange(1.,25.,1.),arange(25.,80.,5.)))
 L_tune= zeros((len(lc),5))    

 pos=0
 for val in lc:
  response_e, response_i, sem_re, sem_ri= find_rate(val)       
  L_tune[pos,0]=val
  L_tune[pos,1]=response_e
  L_tune[pos,2]=response_i
  L_tune[pos,3]=sem_re
  L_tune[pos,4]=sem_ri
  pos+=1   


 
 savefit='C1_%.1f_neuron_%.1f_pref-orientation_%.3f' % (C1,recorded_neuron,theta[recorded_neuron])
 savefit.replace(' ','')

 fitname='path'+ savefit +'.txt' 
 savetxt(fitname,L_tune)

 save_cur=[length,save_e,save_i]
 asave_cur=transpose(array(save_cur))

 savi_cur=[length,savi_e,savi_i]
 asavi_cur=transpose(array(savi_cur))

 savIn_cur=[length,save_in,savi_in]
 asavIn_cur=transpose(array(savIn_cur))

 save_g=[length,save_ge,save_gi]
 asave_g=transpose(array(save_g))

 savi_g=[length,savi_ge,savi_gi]
 asavi_g=transpose(array(savi_g))

 s_gin=[length,sav_gin]
 as_gin=transpose(array(s_gin))
 
 fitname1='/path/LT_Rec_Inputs/E_cells/recurrent_inputs_E_cell_'+ savefit+'.txt'
 savetxt(fitname1,asave_cur)

 fitname2='/path/LT_Rec_Inputs/I_cells/recurrent_inputs_I_cell_'+ savefit+'.txt'
 savetxt(fitname2,asavi_cur)

 fitname3='/path/LT_External_Inputs/external_inputs_'+ savefit+'.txt'
 savetxt(fitname3,asavIn_cur)

 fitname4='/path/LT_Conductances/E_cells/conductances_E_cell_'+ savefit+'.txt'
 savetxt(fitname4,asave_g)

 fitname5='/path/LT_Conductances/I_cells/conductances_I_cell_'+ savefit+'.txt'
 savetxt(fitname5,asavi_g)
 
 fitname6='/path/LT_Conductances/External/external_conductance_'+ savefit+'.txt'
 savetxt(fitname6,as_gin)
