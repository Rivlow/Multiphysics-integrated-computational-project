import numpy as np
import matplotlib.pyplot as plt
import os
import sys

current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))



SMALL_SIZE = 8
MEDIUM_SIZE = 14
BIGGER_SIZE = 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
#plt.rc('text', usetex=True)
#plt.rc('font', family='lmodern')

def Cubic(W,dW):

    for idx, r in enumerate(r_span):
        
        if (r/h < 1.0):
            W[idx] = alpha_cubic * (2/3 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h))
            dW[idx] = alpha_cubic / h * (1.5 * r * r / (h * h) - 2.0 * r / h)
            
        elif (r/h >= 1.0 and r/h <= 2.0):
            W[idx] = alpha_cubic * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h))
            dW[idx] =  alpha_cubic / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h))
            
        else:
            W[idx] = 0
            dW[idx] = 0


def Quintic(W, dW):

    for idx, r in enumerate(r_span):

        factor = (1.0 - 0.5*(r/h))

        W[idx] = alpha_quintic * (factor*factor*factor*factor) * (2.0 * (r/h) + 1.0)
        dW[idx] = -5*alpha_quintic/h * (factor*factor*factor) * (r/h) 
    


s = 0.001
h = 1.2*s
kappa = 2
pi = np.pi
alpha_cubic = 15.0 / (7.0 * pi  * h * h)
alpha_quintic = 7.0 /( 4.0 * pi * h * h )

r_span = np.linspace(0, kappa*h, 100)

W_cubic = np.zeros(len(r_span))
dW_cubic = np.zeros(len(r_span))
Cubic(W_cubic, dW_cubic)

W_quintic = np.zeros(len(r_span))
dW_quintic = np.zeros(len(r_span))
Quintic(W_quintic, dW_quintic)


plt.figure(figsize=(6.9, 4))
#plt.plot(r_span, W_cubic/alpha_cubic, label=r'W(r,h)$ / \alpha$', c='b')
#plt.plot(r_span, dW_cubic/(alpha_cubic/h), ls = "--", label=r'dW(r,h)/dr/($\alpha$/h)', c='r')
#plt.plot(r_span, W_quintic/alpha_quintic, label=r'W(r,h)$ / \alpha$', c='b')
#plt.plot(r_span, dW_quintic/(alpha_quintic/h), ls = "--", label=r'dW(r,h)/dr/($\alpha$/h)', c='r')
plt.plot(r_span, W_cubic, label= "W_cubic", c='b')
plt.plot(r_span, dW_cubic, ls = "--", label= "dW_cubic", c='b')
plt.plot(r_span, W_quintic, label="W_quintic", c='r')
plt.plot(r_span, dW_quintic, ls = "--", label= "dW_quintic", c='r')
plt.xlabel('r/h')
plt.grid(True)
plt.legend(loc = 'best')

plt.tight_layout()
#plt.savefig(f"{current_directory}/Pictures/kernel_cubic_spline.PDF")
plt.show()
        
 
