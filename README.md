# Quantum Diffusion Monte Carlo 
Quantum Diffusion Monte Carlo (QDMC) generator used to simulate the evolution of quantum systems with a given potential function, such as the Simple Harmonic Oscillator (SHO). 

The initial state of the repo was used as a research project, but will be undergoing more refactoring
to improve readability and usability.

## Getting Started

Currently, the script has to be manually compiled depending on which example you wish to run. 
For example, DMCv3SimplePotentials2.cpp contains the toy examples with many known quantum potential systems for verification of the algorithm. The more interesting Holstein Dimer model needs to be compiled separately, DMCv3HolsteinDimerSource.cpp. 

to compile: 

'''
g++ DMCv3SimplePotentials2.cpp -o SimplePotentials
'''
Or

'''
g++ DMCv3HolsteinDimerSource.cpp -o HolsteinDimer
'''

Hopefully this project will be brought up to better standards soon.