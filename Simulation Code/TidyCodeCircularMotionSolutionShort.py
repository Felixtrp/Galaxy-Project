# ---- CODE THAT PRODUCES CIRCULAR MOTION IN GALAXIES WITH AN ORBITING STAR --------------------


#### BEFORE SIMULATION CODE, REPLACE ALL INSTANCES OF 'PARABOLIC_VELOCITIES' WITH 'CIRCULAR_VELOCITIES' AS FOLLOWS:

def circular_velocities(separation, mass_ratio):
    v2sq = 1/(separation*2)
    v2 = np.sqrt(v2sq)
    v1 = v2
    velocities = [-v1, v2]
    return velocities


#### AFTER SIMULATION CODE:

radii= [2]
densities=[1]

time_before_collision = 50
separation = 15
mass_ratio = 1
clockwise = 1
steps = 5001

particles = sum(densities)

particle_coordinates, galaxy1_coordinates, galaxy2_coordinates, time_list, particle_distances = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
plt.figure(figsize=(20, 10))
plt.plot(galaxy1_coordinates[:, 0], galaxy1_coordinates[:, 1], label = 'Galaxy 1 Path')
plt.plot(galaxy2_coordinates[:, 0], galaxy2_coordinates[:, 1], label= 'Galaxy 2 Path')
plt.plot(particle_coordinates[:, 0, 0], particle_coordinates[:, 0, 1], label = 'Particle Path')
plt.axis('equal')
plt.legend(fontsize=16)
plt.title("Galaxies in Circular Motion", fontsize = 20)
plt.xlabel('x-position', fontsize = 16)
plt.ylabel('y-position', fontsize = 16)
plt.savefig("AnalyticSolution.png")

plt.show()


