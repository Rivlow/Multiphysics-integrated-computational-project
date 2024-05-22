import numpy as np
import matplotlib.pyplot as plt

nstepT = 20000


s =  np.array([0.1, 0.05, 0.025, 0.0125, 0.01])
nb_part = np.array([122, 357, 1187, 4192, 6319, 24151])


OMP_time =  np.array([33.10, 45.74, 110.60, 331.61, 460])
serial_time =  np.array([6.98, 24.14, 74.65, 298.96, 465.71])


plt.scatter(s, OMP_time, label = "OMP")
plt.scatter(s, serial_time, label = "serial")
plt.legend(loc = "best")
plt.show()