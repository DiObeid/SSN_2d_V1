# SSN_2d_V1
Stabilized Supralinear Network Model of Responses to Surround Stimuli in Primary Visual Cortex. This repository conatins the following codes,
- LT.py : Computes the length tuning curves of selected cells in the network.
- STC_{}.py: These two codes compute the response to a center and a center+surround stimulus at various relative orientations.
- FSS.py: Computes the local population response to a plaid stimulus at the center and a surround stimulus of orientation matching the orientation of the second component of the plaid.
- LT_parameter.py, STC_parameter.py and FSS_parameter.py contain the paramters for each of the codes above respectively. The user can merge these files into a single one. We decided to create
  a spearte paramter file for each in-silico experiment since some users may only be inetrested in a specific experiment.
- Codes with names starting with "Spikes" are for the Spiking model. These codes are written using Brian simulator (https://briansimulator.org/) (1) Goodman DFM and Brette R (2013) "Brian simulator" Scholarpedia 8(1):10883 and 
  (2)Stimberg M, Goodman DFM, Benichoux V, Brette R (2014) "Equation-oriented specification of neural models for simulations" Frontiers in Neuroinformatics 8:6. doi: 10.3389/fninf.2014.00006.
  
