import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_points_inside_sphere(num_points, radius, num_segments_longitudinaux, num_segments_latitudinaux):
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
            
            points.append((x, y, z, theta_deg_positif, phi_deg_positif, False))  # Set initial color to False
            
            # Generate points on the surfaces of segments
            for i in range(num_segments_longitudinaux):
                for j in range(1, num_segments_latitudinaux):
                    phi_segment = np.radians(i * 360 / num_segments_longitudinaux)
                    theta_segment = np.radians(j * 180 / num_segments_latitudinaux)
                    x_seg = radius * np.cos(phi_segment) * np.sin(theta_segment)
                    y_seg = radius * np.sin(phi_segment) * np.sin(theta_segment)
                    z_seg = radius * np.cos(theta_segment)
                    points.append((x_seg, y_seg, z_seg, np.degrees(theta_segment), np.degrees(phi_segment), True))  # Set color to True for points at intersections
    
    return points

# Parameters
radius = 5
num_points_inside = 8
num_segments_longitudinaux = 8
num_segments_latitudinaux = 4

# Generate points inside the sphere
inside_points = generate_points_inside_sphere(num_points_inside, radius, num_segments_longitudinaux, num_segments_latitudinaux)
inside_points = np.array(inside_points)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='b', alpha=0.3)

# Plot points inside the sphere
for point in inside_points:
    ax.scatter(point[0], point[1], point[2], color='k')

# Plot lines connecting intersections to origin
for point in inside_points:
    ax.plot([0, point[0]], [0, point[1]], [0, point[2]], color='k')

# Plot cones
for i in range(num_segments_longitudinaux):
    phi = np.radians(i * 360 / num_segments_longitudinaux)
    x_vals = radius * np.cos(phi) * np.sin(u)
    y_vals = radius * np.sin(phi) * np.sin(u)
    z_vals = radius * np.cos(u)
    ax.plot(x_vals, y_vals, z_vals, color='k')

for j in range(1, num_segments_latitudinaux):
    theta = np.radians(j * 180 / num_segments_latitudinaux)
    x_vals = radius * np.cos(u) * np.sin(theta)
    y_vals = radius * np.sin(u) * np.sin(theta)
    z_vals = radius * np.cos(theta)
    ax.plot(x_vals, y_vals, z_vals, color='k')

'''
# Add numbering to sectors
for j in range(num_segments_latitudinaux):
    for i in range(num_segments_longitudinaux):
        phi = np.radians(i * 360 / num_segments_longitudinaux)
        theta = np.radians(j * 360 / num_segments_latitudinaux)
        x_text = 1.1 * radius * np.cos(phi+0.5) * np.sin(theta+0.5)  # Adjust position for text
        y_text = 1.1 * radius * np.sin(phi+0.5) * np.sin(theta+0.5)  # Adjust position for text
        z_text = 1.1 * radius * np.cos(theta)  # Adjust position for text
        ax.text(x_text, y_text, z_text, str(j+i + 1), color='r', fontsize=10)
'''

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
