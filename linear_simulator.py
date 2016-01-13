# this class creates and runs a simualtion using linear FR
# the run method returns the population dynamics
# it also implements all the check son the dyanmics that we require for parameter selection

import numpy as np
import matplotlib.pyplot as plt
import math as math
import random as rnd

class linear_simulator_parameter_tester():
    
    def __init__(self, a, b, c, x00, x10, x0s, x1s, plot=False):
        
        self.tsteps = 100000
        self.dt = 0.001
        
        self.a = a
        self.b = b
        self.c = c
        
        self.x0s = x0s
        self.x1s = x1s
        
        self.dynamics = np.zeros((3,self.tsteps+1))
        self.dynamics[0,0] = 0.00
        self.dynamics[1,0] = x00
        self.dynamics[2,0] = x10
        
        self.relaxation_threshold = 0.05  # if within 5% of equilibrium
        self.population_ratio_tol = 10
        self.theta = 0
        self.params_ok = True
        
        self.plot = plot
        self.t = None
        
    def run(self):
        
        #for t in range(self.tsteps):
        finish=False
        t = 0
        while finish==False:
            
            x0 = self.dynamics[1,t]
            x1 = self.dynamics[2,t]
            
            self.dynamics[0,t+1] = self.dynamics[0,t] + self.dt
            self.dynamics[1,t+1] = x0 + self.dt*( self.a*x0*(1.0-x0) - self.b*x1*x0)
            self.dynamics[2,t+1] = x1 + self.dt*( -1.0*x1 + self.c*x0*x1)
            
            if self.check_ratio(t):
                self.params_ok = False
                return False            
            if self.check_relaxation(x0,x1):
                self.params_ok = False
                #print("relaxation too fast")
                #break
                return False
            
            if self.check_theta(t):
                finish=True
                self.t = t
                
            if t>=self.tsteps-1:
                self.params_ok = False
                #finish = True
                return False
            
            t += 1
        
        #print self.theta/np.pi
        if self.plot:
            plt.subplot(1,2,1)    
            plt.plot(self.dynamics[1,0:t],'g')
            plt.plot(self.dynamics[2,0:t],'r')
            plt.subplot(1,2,2)
            plt.plot(self.dynamics[1,0:t], self.dynamics[2,0:t])    
            plt.show()
     
    def get_dynamics(self):
        return self.dynamics[:,0:self.t]
           
    def check_relaxation(self, x0, x1):
        
        
        #print(abs(self.x1s - x1))
        
        if (abs(self.x0s - x0)/self.x0s)<self.relaxation_threshold and (abs(self.x1s - x1)/self.x1s)<self.relaxation_threshold:
            return True
        else:
            return False
    
##    def check_theta(self, t):
##            
##        delta = np.sqrt( (self.dynamics[1,t+1]-self.dynamics[1,t])**2 + (self.dynamics[2,t+1]-self.dynamics[2,t])**2  )
##        r = np.sqrt ( (self.dynamics[1,t]-self.x0s)**2 + (self.dynamics[2,t]-self.x1s)**2 )
##        self.theta += delta/r
##        
##        if self.theta < 4*np.pi:
##            return False
##        else:
##            return True
## using the cosine rule is more accurate than this apprroximation:
    def check_theta(self, t):
        
        x00 = self.dynamics[1,t]
        x01 = self.dynamics[1,t+1]
        x10 = self.dynamics[2,t]
        x11 = self.dynamics[2,t+1]
        
        a = np.sqrt( (x00 - self.x0s)**2 + (x10-self.x1s)**2  )
        b = np.sqrt( (x01 - self.x0s)**2 + (x11-self.x1s)**2  )
        c = np.sqrt( (x01 - x00)**2 + (x11 - x10)**2  )    
        
        self.theta += math.acos( (a**2 + b**2 - c**2) /(2.0*a*b) )
        
        if self.theta < 4*np.pi:
            return False
        else:
            return True

    def check_ratio(self,t):
        
        ratio = self.dynamics[1,t]/self.dynamics[2,t]
        if ratio>self.population_ratio_tol or ratio<(1.0/self.population_ratio_tol):
            return True
        else:
            return False


class linear_simulator():
    
    def __init__(self, a, b, c, x00, x10, dt, T2P, noise_intensity=0.0, plot=False):
        
        self.tsteps = 100000000  # overallocate array intially, trim after simiualtion is finished.
        self.dt = dt
        self.T2P = T2P
        self.noise_intensity = noise_intensity
        
        self.a = a
        self.b = b
        self.c = c
        
        self.dynamics = np.zeros((3,self.tsteps+1))
        self.dynamics[0,0] = 0.00
        self.dynamics[1,0] = x00
        self.dynamics[2,0] = x10
        
        self.plot = plot
        
        self.ext_prey = []  # log the iteration numbers for which prey extinctions occured
        self.ext_pred = []  # likewise for predators - so that these can be excluded from inference calculations later on
        
    def run(self):
        
        i = 0
        t = 0
        while t < self.T2P:
            
            self.check_extinction(i)
            
            x0 = self.dynamics[1,i]
            x1 = self.dynamics[2,i]
            
            self.dynamics[0,i+1] = self.dynamics[0,i] + self.dt
            self.dynamics[1,i+1] = x0 + self.dt*( self.a*x0*(1.0-x0) - self.b*x1*x0) + self.noise(x0)
            self.dynamics[2,i+1] = x1 + self.dt*( -1.0*x1 + self.c*x0*x1) + self.noise(x1)
            
            i += 1
            t += self.dt
        
        self.t = i-1   # to be used for trimming overallocated dynamics array
        
        if self.plot:
            plt.subplot(1,2,1)    
            plt.plot(self.dynamics[0,0:i], self.dynamics[1,0:i],'g')
            plt.plot(self.dynamics[0,0:i], self.dynamics[2,0:i],'r')
            plt.subplot(1,2,2)
            plt.plot(self.dynamics[1,0:i], self.dynamics[2,0:i])    
            plt.show()
    
    def noise(self,X):
        return rnd.gauss(0, self.dt*self.noise_intensity*X)
        
    def check_extinction(self, t):
        
        if self.dynamics[1,t]<=0:
            self.dynamics[1,t] = self.dynamics[1,0]
            self.ext_prey.append(t)
        if self.dynamics[2,t]<=0:
            self.dynamics[2,t] = self.dynamics[2,0]
            self.ext_pred.append(t)

    def get_dynamics(self):
        return self.dynamics[:,0:self.t]
           
if __name__=='__main__':
    
    a = 6.808062820097990908e+00
    b = 3.660581875207621749e+01
    c = 1.449003296756300685e+01
    
    x00 = 6.901295547350187742e-02 / 2.0
    x10 = 1.731478355056400020e-01 / 2.0
    T2P = 5.141000000000051529e+00
    
    hs = linear_simulator(a, b, c, x00, x10, 0.0001, T2P, 0.0, True)
    hs.run()
    print(len(hs.ext_prey))
    print(len(hs.ext_pred))
    
    D = hs.get_dynamics()
    E_prey = np.asarray(hs.ext_prey)
    E_pred = np.asarray(hs.ext_pred)
    np.savetxt("example_run.dynamics", D, delimiter=',')
    np.savetxt("example_prey.extinctions", E_prey, delimiter=',')
    np.savetxt("example_pred.extinctions", E_pred, delimiter=',')
