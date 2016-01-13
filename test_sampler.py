# THIS SCRIPT:
# Tests repeated simulation for a single parameter set. Does Timme inference, saves result. Times whole process.
## Varies: number of samples (taken from two oscillation dynamic). 
## for pID=0 can go up to 50,000 samples. Do from 500 to 50,000 in steps of 500..

from datetime import datetime
from configure import NOISE_VALUE, NUMBER_OF_PARAMETERS, PARAMETER_FILE, DT, REPEATS_PER_PARAMETER_SET, NUMBER_OF_SAMPLES, NUMBER_OF_BINS
from linear_simulator import linear_simulator
from sampler import sampler
from timme_calculator import timme_calculator
import numpy as np
import os, sys


parameters = np.genfromtxt(PARAMETER_FILE, delimiter=',')
p = 0  ## pID

noise = 50
#samples = range(500,50500,500) 
#samples = range(50,5050,50)
#samples = np.logspace(2,15,num=14,base=2)
#samples = samples.astype(int)  
samples = range(3,20)

a = parameters[p,1]
b = parameters[p,2]
c = parameters[p,3]
       
x00 = parameters[p,13]
x10 = parameters[p,14]
T2P = parameters[p,12]
       
print "b = ", -b

for si in samples:
        

        ls = linear_simulator(a, b, c, x00, x10, DT, T2P, noise, False)
        ls.run()
        D = ls.get_dynamics()
        
        S = sampler(si, D)
        S.sample()
        calc = timme_calculator(S, NUMBER_OF_BINS)
	results, err = calc.calculate()
	results = results[0:6]

	print(results[2])
	#print(S.sampled_dynamics)


