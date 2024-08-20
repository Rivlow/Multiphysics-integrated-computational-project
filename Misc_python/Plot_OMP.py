import numpy as np
import matplotlib.pyplot as plt
import os
import sys

SMALL_SIZE = 8
MEDIUM_SIZE = 14
BIGGER_SIZE = 18
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)   # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)   # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rc('text', usetex=True)
plt.rc('font', family='lmodern')
    
current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))


threads = [1, 5, 10, 22]
particles = [2508, 9330, 27498, 377538]
time_1 = [194.38, 963.52, 3552.42, 83073]
time_5 = [75.99, 569.25, 1290, 26500]
time_10 = [73.80, 571.64, 1146.862, 18674]
time_22 = [78.35, 274.44, 1067.17, 13319]

tab_1 = [time_1[0], time_5[0], time_10[0], time_22[0]]
tab_5 = [time_1[1], time_5[1], time_10[1], time_22[1]]
tab_10 = [time_1[2], time_5[2], time_10[2], time_22[2]]
tab_22 = [time_1[3], time_5[3], time_10[3], time_22[3]]

plot_threads = True
plot_time = True

if plot_threads:
    plt.figure(figsize=(7, 5))
    plt.xlabel("Number of threads [-]")
    plt.ylabel("Total computation time [s]")
    plt.yscale("log")
    plt.scatter(threads, tab_1, label=f"{particles[0]} particles")
    plt.plot(threads, tab_1, ls="--")
    plt.scatter(threads, tab_5, label=f"{particles[1]} particles")
    plt.plot(threads, tab_5, ls="--")
    plt.scatter(threads, tab_10, label=f"{particles[2]} particles")
    plt.plot(threads, tab_10, ls="--")
    plt.scatter(threads, tab_22, label=f"{particles[3]} particles")
    plt.plot(threads, tab_22, ls="--")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    plt.savefig(f"{current_directory}/../Pictures/Scaling_OMP_threads.PDF")

if plot_time:
    plt.figure(figsize=(7, 5))
    plt.xlabel("Number of particles [-]")
    plt.ylabel("Total computation time [s]")
    plt.yscale("log")
    plt.xscale("log")
    plt.scatter(particles, time_1, label=f"{threads[0]} threads")
    plt.plot(particles, time_1, ls="--")
    plt.scatter(particles, time_5, label=f"{threads[1]} threads")
    plt.plot(particles, time_5, ls="--")
    plt.scatter(particles, time_10, label=f"{threads[2]} threads")
    plt.plot(particles, time_10, ls="--")
    plt.scatter(particles, time_22, label=f"{threads[3]} threads")
    plt.plot(particles, time_22, ls="--")
    plt.grid(True, alpha = 0.5)
    plt.legend(loc = "best")
    plt.tight_layout()
    plt.savefig(f"{current_directory}/../Pictures/Scaling_OMP_particles.PDF")

plt.show()



