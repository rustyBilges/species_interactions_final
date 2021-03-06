## produces plot of estimated parameter values; variance; and error functions versus noise intensity
## FOR A SINGLE PARAMETER SET:
## 	p0_versus_noise.txt

import numpy as np
import matplotlib.pyplot as plt

plot_type = "line"
xlog = True 

pID = 0
parameters = np.genfromtxt("parameter_log.txt", delimiter=',')
#results = np.genfromtxt("test_single_params_sample_size_1.txt", delimiter=',')
#results = np.genfromtxt("test_single_params_sample_size_2.txt", delimiter=',')
results = np.genfromtxt("test_single_params_sample_size_3.txt", delimiter=',')

fsa = 15

#fig, axes = plt.subplots(1,3, figsize=(15,6))
fig, axes = plt.subplots(1,2, figsize=(10,6))
AX = axes.flatten()

if plot_type=="line":
	AX[0].plot(results[:,0], results[:,1], 'r-', label='$\hat{r}_0$')
	AX[0].axhline(parameters[0,1], color='r', linestyle='--')

	AX[0].plot(results[:,0], results[:,2], 'b-', label='$\hat{J}_{00}$')
	AX[0].axhline(-parameters[0,1], color='b', linestyle='--')

	AX[0].plot(results[:,0], results[:,3], 'g-', label='$\hat{J}_{01}$')
	AX[0].axhline(-parameters[0,2], color='g', linestyle='--')

	AX[0].plot(results[:,0], results[:,4], 'k-', label='$\hat{r}_{1}$')
	AX[0].axhline(-1.0, color='k', linestyle='--')

	AX[0].plot(results[:,0], results[:,5], 'm-', label='$\hat{J}_{10}$')
	AX[0].axhline(parameters[0,3], color='m', linestyle='--')

	AX[0].plot(results[:,0], results[:,6], 'c-', label='$\hat{J}_{11}$')
	AX[0].axhline(0, color='c', linestyle='--')

	AX[0].grid()
	AX[0].set_xlabel('number of samples', fontsize = fsa)
	AX[0].set_ylabel('mean estimated parameter value', fontsize = fsa)

	if xlog:
		AX[0].set_xscale('log', basex=2)

	AX[1].plot(results[:,0], results[:,1+8], 'r-', label='$\hat{r}_0$')
	AX[1].plot(results[:,0], results[:,2+8], 'b-', label='$\hat{J}_{00}$')
	AX[1].plot(results[:,0], results[:,3+8], 'g-', label='$\hat{J}_{01}$')
	AX[1].plot(results[:,0], results[:,4+8], 'k-', label='$\hat{r}_{1}$')
	AX[1].plot(results[:,0], results[:,5+8], 'm-', label='$\hat{J}_{10}$')
	AX[1].plot(results[:,0], results[:,6+8], 'c-', label='$\hat{J}_{11}$')
	AX[1].legend(loc=2, fontsize=fsa)
	AX[1].grid()
	AX[1].set_xlabel('number of samples', fontsize = fsa)
	AX[1].set_ylabel('variance in estimates', fontsize = fsa)
	if xlog:
		AX[1].set_xscale('log', basex=2)


if plot_type=="scatter":
	AX[0].plot(results[:,0], results[:,1], 'ro', label='$\hat{r}_0$')
	AX[0].axhline(parameters[0,1], color='r', linestyle='--')

	AX[0].plot(results[:,0], results[:,2], 'bo', label='$\hat{J}_{00}$')
	AX[0].axhline(-parameters[0,1], color='b', linestyle='--')

	AX[0].plot(results[:,0], results[:,3], 'go', label='$\hat{J}_{01}$')
	AX[0].axhline(-parameters[0,2], color='g', linestyle='--')

	AX[0].plot(results[:,0], results[:,4], 'ko', label='$\hat{r}_{1}$')
	AX[0].axhline(-1.0, color='k', linestyle='--')

	AX[0].plot(results[:,0], results[:,5], 'mo', label='$\hat{J}_{10}$')
	AX[0].axhline(parameters[0,3], color='m', linestyle='--')

	AX[0].plot(results[:,0], results[:,6], 'co', label='$\hat{J}_{11}$')
	AX[0].axhline(0, color='c', linestyle='--')

	AX[0].grid()
	AX[0].set_xlabel('number of samples', fontsize = fsa)
	AX[0].set_ylabel('mean estimated parameter value', fontsize = fsa)
	if xlog:
		AX[0].set_xscale('log', basex=2)


	AX[1].plot(results[:,0], results[:,1+8], 'r-', label='$\hat{r}_0$')
	AX[1].plot(results[:,0], results[:,2+8], 'b-', label='$\hat{J}_{00}$')
	AX[1].plot(results[:,0], results[:,3+8], 'g-', label='$\hat{J}_{01}$')
	AX[1].plot(results[:,0], results[:,4+8], 'k-', label='$\hat{r}_{1}$')
	AX[1].plot(results[:,0], results[:,5+8], 'm-', label='$\hat{J}_{10}$')
	AX[1].plot(results[:,0], results[:,6+8], 'c-', label='$\hat{J}_{11}$')
	AX[1].legend(loc=2, fontsize=fsa)
	AX[1].grid()
	AX[1].set_xlabel('number of samples', fontsize = fsa)
	AX[1].set_ylabel('variance in estimates', fontsize = fsa)
	if xlog:
		AX[1].set_xscale('log', basex=2)

#nsamples = results[:,0]
#nsamples = nsamples.astype(float)
#AX[2].errorbar(results[:,0], results[:,7]/nsamples, yerr=np.sqrt(results[:,7+8])/nsamples , fmt='ro', label='error 1')
#AX[2].errorbar(results[:,0], results[:,8]/nsamples, yerr=np.sqrt(results[:,8+8])/nsamples , fmt='bo', label='error 2')
#AX[2].legend(loc=2, fontsize=fsa)
#AX[2].grid()
#AX[2].set_xlabel('number of samples', fontsize = fsa)
#AX[2].set_ylabel('error function value', fontsize = fsa)

plt.tight_layout()
plt.show()
