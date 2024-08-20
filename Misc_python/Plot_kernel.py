import numpy as np
import matplotlib.pyplot as plt
import os
import sys

current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))


def isLatex(latex):
    if latex:
        
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
        plt.rc('text', usetex=True)
        plt.rc('font', family='lmodern')     



def plotKernel(plot, save, dimension):

    # Assign normalization coefficient
    if dimension == 2:
        alpha_cubic = 15.0 / (7.0 * pi  * h * h)
        alpha_quintic = 7.0 /( 4.0 * pi * h * h )
        alpha_coh = 40/(pi*h*h*h*h*h*h*h*h)
        alpha_adh = 16/(4* pi*np.power(h, 2.25))

    elif dimension == 3:
        alpha_cubic = 3.0 / (2.0 * pi * h * h * h)
        alpha_quintic = 7.0 /( 4.0 * pi * h * h * h)
        alpha_coh = 32/(pi*h*h*h*h*h*h*h*h*h)
        alpha_adh = 0.0007/np.power(h,3.25)


    # Fill value for the kernels and their derivatives
    for idx, r in enumerate(r_span):

        if (r/h < kappa/2):
            W_cubic[idx] = alpha_cubic * (2/3 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h))
            dW_cubic[idx] = alpha_cubic / h * (1.5 * r * r / (h * h) - 2.0 * r / h)
            W_coh[idx] = alpha_coh*(2.0 * (kappa*h-r) * (kappa*h-r) * (kappa*h-r) * (r*r*r)- np.power(kappa*h, 6)/64)
            dW_coh[idx] = alpha_coh * 6*(kappa*h-r)*(kappa*h-r)*(kappa*h-r)* (r*r) * (kappa*h-2*r)
                      
        elif (r/h >= kappa/2 and r/h < kappa):
            W_cubic[idx] = alpha_cubic * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h))
            dW_cubic[idx] =  alpha_cubic / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h))
            W_coh[idx] = alpha_coh*(kappa*h-r)*(kappa*h-r)*(kappa*h-r)*r*r*r
            dW_coh[idx] = alpha_coh * 3*(kappa*h-r)*(kappa*h-r)*(kappa*h-r) * (r*r) * (kappa*h-2*r)
            W_adh[idx] = alpha_adh* np.power(-4*(r*r)/(kappa*h) + 6*r - 2*(kappa*h), 1/4)
            dW_adh[idx] = alpha_adh * 0.25 * (-8*r/(kappa*h) + 6)/np.power(-4*(r*r)/(kappa*h) + 6*r - 2*(kappa*h), 3/4)
            
        factor = (1.0 - 0.5*(r/h))
        W_quintic[idx] = alpha_quintic * (factor*factor*factor*factor) * (2.0 * (r/h) + 1.0)
        dW_quintic[idx] = -5*alpha_quintic/h * (factor*factor*factor) * (r/h) 

    # Plot desired kernels
    if plot["Cubic"]:
        plt.figure()
        plt.plot(r_span, W_cubic/alpha_cubic, label=r'W(r,h)$ / \alpha$', c='b')
        plt.plot(r_span, dW_cubic/(alpha_cubic/h), ls = "--", label=r'dW(r,h)/dr/($\alpha$/h)', c='r')
        plt.xlabel('r/h')
        plt.grid(True)
        plt.legend(loc = 'best')
        if save:
            plt.savefig(rf"{current_directory}\Pictures\cubic_kernel.PDF")

    if plot["Quintic"]:
        plt.figure()
        plt.plot(r_span, W_quintic/alpha_quintic, label=r'W(r,h)$ / \alpha$', c='b')
        plt.plot(r_span, dW_quintic/(alpha_quintic/h), ls = "--", label=r'dW(r,h)/dr/($\alpha$/h)', c='r')
        plt.xlabel('r/h')
        plt.grid(True)
        plt.legend(loc = 'best')
        if save:
            plt.savefig(rf"{current_directory}\Pictures\quintic_kernel.PDF")

    if plot["Adhesion"]:
        plt.figure()
        plt.plot(r_span, W_adh/alpha_adh, label=r'$W(r,h)/0.007/h^{3.25}$', c='b')
        plt.xlabel('r/h')
        plt.grid(True)
        plt.legend(loc = 'best')
        if save:
            plt.savefig(rf"{current_directory}\Pictures\adh_kernel.PDF")

    if plot["Cohesion"]:
        plt.figure()
        plt.plot(r_span, W_coh/alpha_coh, label=r'$W(r,h)/32/h^9$', c='b')
        plt.xlabel('r/h')
        plt.grid(True)
        plt.legend(loc = 'best')
        if save:
            plt.savefig(rf"{current_directory}\Pictures\coh_kernel.PDF")



# Geometrical parameters
s = 0.025
h = 1.2*s
kappa = 2
pi = np.pi
dimension = 3

r_span = np.linspace(0, kappa*h, 100)
W_cubic = np.zeros(len(r_span))
dW_cubic = np.zeros(len(r_span))
W_quintic = np.zeros(len(r_span))
dW_quintic = np.zeros(len(r_span))
W_adh = np.zeros(len(r_span))
dW_adh = np.zeros(len(r_span))
W_coh = np.zeros(len(r_span))
dW_coh = np.zeros(len(r_span))

plot = {"Cubic":False, "Quintic":False, "Adhesion":True, "Cohesion":True}
save = False
isLatex(False)
plotKernel(plot, save, dimension)
plt.show()
        
