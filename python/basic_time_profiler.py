import numpy as np
import matplotlib.pyplot as plt


particle = np.array([125, 1331, 3375, 5832, 9261])
linked_time = np.array([0.000802, 0.054912, 0.372272, 1.11271, 2.88625])
naive_time = np.array([0.000332, 0.019394, 0.11395, 0.336295, 0.833985])

plt.scatter(particle, linked_time, color = 'red', label='linked list algorithm')
plt.scatter(particle, naive_time, color = 'blue', label='naive algorithm')
plt.xlabel('number of particles [-]')
plt.ylabel('elapsed time [s]')

plt.legend(loc = 'best')
plt.savefig("time_performance.PDF")
plt.show()
