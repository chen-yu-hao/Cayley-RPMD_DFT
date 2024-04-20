#!/usr/bin/env python 
# encoding: utf-8 
from NEW import *
get_potential=pes.get_potential
initialize_potential=pes.initialize_potential



label = 'NH4' 
tranums=200
evotime=30
eqi=10
reactants( 
    atoms = ['H', 'H', 'H', 'H', 'N'],
    reactant1Atoms = [1,2], 
    reactant2Atoms = [3,4,5], 
    Rinf = (30,"angstrom"), 
) 
 
transitionState( 
    geometry = ([[ 4.59141283e+00,  1.12784325e-04, -1.48385126e-04],
                [ 2.85382609e+00, -1.62554847e-04,  2.13923384e-04],
                [-3.87636146e-01,  1.90555496e+00,  6.12161992e-05],
                [-3.87635850e-01, -5.11826143e-01, -1.83553348e+00],
                [ 1.99457004e-04,  5.82583769e-05, -7.86309265e-05]],
        "bohr", 
    ), 
    formingBonds = [(2,5)] ,
    breakingBonds = [(1,2)],
) 

equivalentTransitionState( 
    formingBonds=[(1,5)],
    breakingBonds=[(2,1)], )


#thermostat('GLE', A=('gle_A.txt','s^-1')) 

thermostat('Andersen') 
#thermostat('Andersen',samplingTime=(150,'fs'))
################################################################################# 
xi_list = numpy.arange(-0.05, 1.05, 0.01)

generateUmbrellaConfigurations(  
    dt = (0.1,"fs"), 
    evolutionTime = (10,"ps"), 
    xi_list = xi_list, 
    kforce = 0.2 * T, 
    #kforce = 0.0 * T, 
)

windows = [] 

for xi in numpy.arange(-0.05, 1.05, 0.01): 
    window = Window(xi=xi, kforce=0.10*T, trajectories=tranums, equilibrationTime=(eqi,"ps"), evolutionTime=(evotime,"ps")) 
    #window = Window(xi=xi, kforce=0.0*T, trajectories=1, equilibrationTime=(1,"ps"), evolutionTime=(1,"ps"))  ## Energy conversation test only!!

    windows.append(window) 

conductUmbrellaSampling( 
    dt = (0.0005,"ps"), 
    windows = windows, 
    saveTrajectories = False 
    ) 

computePotentialOfMeanForce(windows=windows, xi_min=-0.05, xi_max=1.05, bins=20000) 

computeRecrossingFactor( 
    dt = (0.5,"fs"), 
    equilibrationTime = (50,"ps"), 
    childTrajectories = 128000, 
    childSamplingTime = (0.5,"ps"), 
    childrenPerSampling =96*4, 
    childEvolutionTime = (40,"ps"), 
 ) 

computeRateCoefficient() 
