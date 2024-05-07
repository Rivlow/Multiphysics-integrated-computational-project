import numpy as np
import matplotlib.pyplot as plt

s = 0.25
h = 1.2*s
kappa = 2


pi = np.pi
alpha = 3.0 / (2.0 * pi * h * h * h)
r_span = np.linspace(0, kappa*h, 100)
W = np.zeros(len(r_span))
dW = np.zeros(len(r_span))

for idx, r in enumerate(r_span):
    
    if (r/h < 1.0):
        W[idx] = alpha * (2/3 - r * r / (h * h) + 0.5 * r * r * r / (h * h * h))
        dW[idx] = alpha / h * (1.5 * r * r / (h * h) - 2.0 * r / h)
        
    elif (r/h >= 1.0 and r/h <= 2.0):
        W[idx] = alpha * 1.0/6.0 * ((2.0 - r / h) * (2.0 - r / h) * (2.0 - r / h))
        dW[idx] =  alpha / h * (-0.5 * (2.0 - r / h) * (2.0 - r / h))
        
    else:
        W[idx] = 0
        dW[idx] = 0


plt.figure()
plt.plot(r_span, W/alpha, label=r'W(r,h)$ / \alpha$', c='b')
plt.plot(r_span, dW/(alpha/h), ls = "--", label=r'dW(r,h)/dr/($\alpha$/h)', c='r')
plt.xlabel('r/h')
plt.grid(True)
plt.legend(loc = 'best')

plt.tight_layout()
plt.show()
plt.savefig("Pictures/kernel_cubic_spline.PDF")
        
 
