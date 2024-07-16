import numpy as np
import os
import sys
import matplotlib.pyplot as plt

def isLatex(latex):
    if latex:
        SMALL_SIZE = 8
        MEDIUM_SIZE = 14
        BIGGER_SIZE = 18
        plt.rc('font', size=SMALL_SIZE)
        plt.rc('axes', titlesize=MEDIUM_SIZE)
        plt.rc('axes', labelsize=MEDIUM_SIZE)
        plt.rc('xtick', labelsize=MEDIUM_SIZE)
        plt.rc('ytick', labelsize=MEDIUM_SIZE)
        plt.rc('legend', fontsize=MEDIUM_SIZE)
        plt.rc('figure', titlesize=BIGGER_SIZE)
        plt.rc('text', usetex=True)
        plt.rc('font', family='lmodern')

# EDO system
def f(t, var):
    x, y = var
    dx_dt = x*y*t
    dy_dt =  - 2 * x + np.sin(y)
    return [dx_dt, dy_dt]

# Explicit Euler method
def euler_explicit(f, xy0, t):
    y = np.zeros((len(t), len(xy0)))
    y[0] = xy0
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        y[i] = y[i-1] + dt * np.array(f(t[i-1], y[i-1]))
    return y

# RK22 method
def runge_kutta_2(f, xy0, t):
    y = np.zeros((len(t), len(xy0)))
    y[0] = xy0
    for i in range(1, len(t)):
        dt = t[i] - t[i-1]
        k1 = dt * np.array(f(t[i-1], y[i-1]))
        k2 = dt * np.array(f(t[i-1] + dt, y[i-1] + k1))
        y[i] = y[i-1] + 0.5 * (k1 + k2)
    return y

def main():
    
    isLatex(True)
    current_directory = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    plot_solution = True
    plot_error = False
    
    
    
    xy0 = [1, 0.5]  # Initial conditions
    #steps = np.array([10, 50, 100, 200, 500, 1000, 5000, 10000, 50000, 100000]) # Different step sizes for time vector
    steps = np.array([10, 100, 1000, 10000]) # Different step sizes for time vector

    if plot_solution:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    mean_errors_x, mean_errors_y = [], []  


    for i, step in enumerate(steps):
        t = np.linspace(0, 10, step) 

        # Solver methods
        sol_euler = euler_explicit(f, xy0, t)
        sol_rk22 = runge_kutta_2(f, xy0, t)
        
        if plot_error:
            # Compute the mean error
            mean_error_x = np.mean(np.abs((sol_rk22[:, 0] - sol_euler[:, 0])/sol_rk22[:, 0]))
            mean_error_y = np.mean(np.abs((sol_rk22[:, 1] - sol_euler[:, 1])/sol_rk22[:, 1]))
            mean_errors_x.append(mean_error_x)
            mean_errors_y.append(mean_error_y)


        if plot_solution:

            row = i // 2
            col = i % 2
            ax = axes[row, col]
            ax.scatter(t, sol_euler[:, 0], s = 6, label='Euler explicit (x)', color="red")
            ax.scatter(t, sol_euler[:, 1], s = 6, label='Euler explicit (y)', color="red", linestyle=":")
            ax.plot(t, sol_rk22[:, 0], label='Runge-Kutta 22 (x)', color="blue")
            ax.plot(t, sol_rk22[:, 1], label='Runge-Kutta 22 (y)', color="blue", linestyle="--")
            ax.set_xlabel('time t')
            ax.set_title(f'Timestep = {10/step}')
            ax.legend(loc="best")
            ax.grid(True)
    
    if plot_error:

        plt.figure()
        plt.yscale('log')
        plt.xscale('log')

        plt.scatter(steps, mean_errors_x, label = "x variable", color = 'red')
        plt.plot(steps, mean_errors_x, color = 'red', linestyle = "--")

        plt.scatter(steps, mean_errors_y, label = "y variable", color = "blue")
        plt.plot(steps, mean_errors_y, color = 'blue', linestyle = "--")

        plt.xlabel("Number of steps")
        plt.ylabel("relative error")
        plt.legend(loc = "best")
        plt.savefig(f"{current_directory}/Pictures/Euler_RK22_relative_error.PDF")
        plt.show()
        

    
    
    
    plt.tight_layout()
    #plt.savefig(f"{current_directory}/Pictures/Euler_RK22_comparison.PDF")
    plt.show()

if __name__ == "__main__":
    main()