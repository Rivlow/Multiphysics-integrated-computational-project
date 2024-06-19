import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_random_points_3d(n):
    points = []
    for _ in range(n):
        u = np.random.rand()
        v = np.random.rand()
        theta = 2 * np.pi * u
        phi = np.arccos(2 * v - 1)
        r = np.cbrt(np.random.rand())  
        x = r * np.sin(phi) * np.cos(theta)
        y = r * np.sin(phi) * np.sin(theta)
        z = r * np.cos(phi)
        points.append((x, y, z))
    return points

def determine_sector(points):
    sectors = []
    for x, y, z in points:
        theta = np.arctan2(y, x)
        if theta < 0:
            theta += 2 * np.pi
        phi = np.arccos(z / np.sqrt(x**2 + y**2 + z**2))
        
        theta_sector = int(theta // (2 * np.pi / 8))
        phi_sector = int(phi // (np.pi / 4))
        
        sector_number = phi_sector * 8 + theta_sector + 1
        sectors.append(sector_number)
    return sectors

def plot_sphere_with_sectors(points, sector_numbers):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('auto')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones(np.size(u)), np.cos(v))

    ax.plot_surface(x, y, z, color='b', alpha=0.1, rstride=5, cstride=5, linewidth=0)

    for i in range(8):
        theta = i * (2 * np.pi / 8)
        x = np.cos(theta) * np.sin(v)
        y = np.sin(theta) * np.sin(v)
        z = np.cos(v)
        ax.plot(x, y, z, color='k', linewidth=0.5)

    for j in range(1, 4):
        phi = j * (np.pi / 4)
        x = np.outer(np.cos(u), np.sin(phi))
        y = np.outer(np.sin(u), np.sin(phi))
        z = np.outer(np.ones(np.size(u)), np.cos(phi))
        ax.plot_wireframe(x, y, z, color='k', linewidth=0.5, rstride=10, cstride=10)
    
    x_points, y_points, z_points = zip(*points)
    ax.scatter(x_points, y_points, z_points, s=20, c='red')

    for i, (x, y, z) in enumerate(points):
        ax.text(x, y, z, str(sector_numbers[i]), color='blue', fontsize=9, ha='center', va='center')

    max_range = 1.0
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range, max_range])

def main(n):
    points = generate_random_points_3d(n)
    sector_numbers = determine_sector(points)
    
    print(sector_numbers)
    
    for i, (part, point) in enumerate(zip(sector_numbers, points)):
        print(f"Point {i} : sector {part}")
    
    plot_sphere_with_sectors(points, sector_numbers)
    plt.show()
    
    
    free_surf = np.zeros(32)
    
    for val in sector_numbers:
        free_surf[val] += 1
        
    for i in range(len(free_surf)):
        if free_surf[i] != 0:
            print(f"Particules present in sector {i}")
            
    print(free_surf)

main(6)