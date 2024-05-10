import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import math

def generate_points_inside_sphere(num_points, radius):
    points = []
    np.random.seed(5)
    while len(points) < num_points:
        # Generate random points inside a cube
        x = np.random.uniform(-radius, radius)
        y = np.random.uniform(-radius, radius)
        z = np.random.uniform(-radius, radius)
        
        # Check if the point is inside the sphere
        if x**2 + y**2 + z**2 <= radius**2:
            # Calculate spherical coordinates
            r = np.sqrt(x**2 + y**2 + z**2)
            theta = np.arccos(z / r)  # Polar angle
            phi = np.arctan2(y, x)     # Azimuthal angle
            
            # Convert angles to degrees and ensure they are positive
            theta_deg = np.degrees(theta)
            phi_deg = np.degrees(phi)
            
            # Ensure angles are positive and within [0, 360)
            theta_deg_positif = (theta_deg + 360) % 360
            phi_deg_positif = (phi_deg + 360) % 360
            
            points.append((x, y, z, theta_deg_positif, phi_deg_positif))
    
    return points

def get_cone_index(point, num_cones):
    # Convert Cartesian coordinates to spherical coordinates
    r = np.linalg.norm(point)
    theta = np.arccos(point[2] / r)
    phi = np.arctan2(point[1], point[0])
    
    # Convert phi to [0, 2*pi)
    phi = (phi + 2 * np.pi) % (2 * np.pi)
    
    # Calculate cone index
    cone_index = int(np.floor(phi / (2 * np.pi / num_cones)))
    return cone_index


# Parameters
radius = 5
num_points_inside = 8
num_cones = 32

# Generate sphere surface points
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

# Generate points inside the sphere
inside_points = generate_points_inside_sphere(num_points_inside, radius)
inside_points = np.array(inside_points)


xyz = inside_points[:,:-2]
theta_arr = inside_points[:, 3]
phi_arr = inside_points[:, 4]

print(xyz)
print('\n')
print(theta_arr)
print('\n')
print(phi_arr)

hash = []
matrix = np.zeros((4, 8))
print(matrix[1])


for theta, phi in zip(theta_arr, phi_arr):
    
    a = int(theta//45)
    b = int(phi//45)
    print(a)
    print(b)
    
    matrix[a][b] += 1
    
print(matrix)




plot = True
# Plot
if (plot):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot sphere
    ax.plot_surface(x, y, z, color='b', alpha=0.3)

    # Plot points inside the sphere
    ax.scatter(inside_points[:,0], inside_points[:,1], inside_points[:,2], color='r')

    # Plot cones

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

    # Assign points to cones
    cone_assignments = [get_cone_index(point, num_cones) for point in inside_points]
    print("Cone assignments for each point:", cone_assignments)
