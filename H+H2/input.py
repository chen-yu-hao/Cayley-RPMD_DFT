#!/usr/bin/env python 
# encoding: utf-8 
  
from h3_bkmp2 import get_potential 
 
################################################################################
 
label = 'H + H2 -> HH + H' 
tranums=1
evotime=50
eqi=10
reactants( 
    atoms = ['H', 'H', 'H'],
    reactant1Atoms = [1,2], 
    reactant2Atoms = [3], 
    Rinf = (30 * 0.52918,"angstrom"), 
) 
 
transitionState( 
    geometry = ( 
        [[  0.000000,  0.000000, -1.7570],  
         [  0.000000,  0.000000,  0.000000],  
         [  0.000000,  0.000000,  1.7570]], 
        "bohr", 
    ), 
    formingBonds = [(2,3)],  
    breakingBonds = [(1,2)], 
) 
equivalentTransitionState( 
    formingBonds=[(1,3)],  
    breakingBonds=[(2,1)], 
) 
 
#thermostat('GLE', A=('gle_A.txt','s^-1')) 
 
################################################################################# 
 

thermostat('Andersen') 
#thermostat('Andersen',samplingTime=(150,'fs'))
################################################################################# 
xi_list = numpy.arange(-0.05, 0.90, 0.01)
xi_list = numpy.append(xi_list,numpy.arange(0.905, 0.989, 0.005))
xi_list = numpy.append(xi_list,numpy.arange(0.988, 1.012, 0.002))
xi_list = numpy.append(xi_list,numpy.arange(1.012, 1.05, 0.005))

#xi_list = [0.0,1.0] 
#xi_list = numpy.arange(-0.05, 1.05, 0.01)
#xi_list=[1.0] 
## Can be used as 1-bead energy conservation test.
generateUmbrellaConfigurations(  
    dt = (0.1,"fs"), 
    evolutionTime = (20,"ps"), 
    xi_list = xi_list, 
    kforce = 0.5 * T, 
    #kforce = 0.0 * T, 
)

################################################################################
#xi_list = [-0.05]
#xi_list = [1.0]
xi_list = numpy.arange(-0.05, 1.05, 0.01) 
windows = [] 

# for xi in numpy.arange(-0.05, 1.05, 0.01): 
#     window = Window(xi=xi, kforce=0.1*T, trajectories=tranums, equilibrationTime=(eqi,"ps"), evolutionTime=(evotime,"ps")) 
#     #window = Window(xi=xi, kforce=0.0*T, trajectories=1, equilibrationTime=(1,"ps"), evolutionTime=(1,"ps"))  ## Energy conversation test only!!

#     windows.append(window) 


for xi in numpy.arange(-0.02, 1.05, 0.01): 
    window = Window(xi=xi, kforce=0.10*T, trajectories=tranums, equilibrationTime=(eqi,"ps"), evolutionTime=(evotime,"ps")) 
    #window = Window(xi=xi, kforce=0.0*T, trajectories=1, equilibrationTime=(1,"ps"), evolutionTime=(1,"ps"))  ## Energy conversation test only!!

    windows.append(window) 

conductUmbrellaSampling( 
    dt = (0.0005,"ps"), 
    windows = windows, 
    saveTrajectories = False 
    ) 

# computePotentialOfMeanForce(windows=windows, xi_min=-0.02, xi_max=1.02, bins=50000) 
ABF()
# computeRecrossingFactor( 
#     dt = (0.5,"fs"), 
#     equilibrationTime = (50,"ps"), 
#     childTrajectories = 128000, 
#     childSamplingTime = (0.5,"ps"), 
#     childrenPerSampling =96*4, 
#     childEvolutionTime = (40,"ps"), 
#  ) 

# computeRateCoefficient() 
