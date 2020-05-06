## TAIL SIZE ANALYSIS - 

# define function that counts the fraction of stars in the tail

# For CLOCKWISE interactions
def find_tail_size_cw(inner_limit, steps, particles, particle_coordinates, galaxy_coordinates, dr):
    tail_size = np.zeros(steps-2)
    for i in range(steps - 2):
        for j in range(particles):
            if dr[i, j] > inner_limit: # Star must be outside a certain inner radial limit from galaxy
                if particle_coordinates[i, j, 1]  > galaxy_coordinates[i, 1]: # Star must have y-coordinate greater than galaxy
                    tail_size[i] += 1
        tail_size[i] = tail_size[i]/particles
    return tail_size

# For ANTICLOCKWISE interactions
def find_tail_size_acw(inner_limit, steps, particles, particle_coordinates, galaxy_coordinates, dr):
    tail_size = np.zeros(steps-2)
    for i in range(steps - 2):
        for j in range(particles):
            if dr[i, j] > inner_limit: # Star must be outside a certain inner radial limit from galaxy
                if particle_coordinates[i, j, 0]  > galaxy_coordinates[i, 0]: # Star must have x-coordinate greater than galaxy
                    tail_size[i] += 1
        tail_size[i] = tail_size[i]/particles
    return tail_size


# Initial Conditions
radii = np.zeros(50)
densities = np.zeros(50)

for i in range(50):
    radii[i] = 2 + (0.1*i)
    densities[i] = 20 + i

densities = densities.astype(int)


time_before_collision = 100
separation = 8
mass_ratio = 2
steps = 3000
particles = sum(densities)

# define an inner limit
inner_limit = 10

# Simulate CLOCKWISE
clockwise = 1
particle_coordinates_cw, galaxy1_coordinates_cw, galaxy2_coordinates_cw, time_list_cw, particle_distances_cw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)

# find size of tail over tail
tail_size_cw = find_tail_size_cw(inner_limit, steps, particles, particle_coordinates_cw, galaxy1_coordinates_cw, particle_distances_cw)

# Simulate ANTICLOCKWISE
clockwise = -1
particle_coordinates_acw, galaxy1_coordinates_acw, galaxy2_coordinates_acw, time_list_acw, particle_distances_acw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)

# find size of tail over tail
tail_size_acw = find_tail_size_acw(inner_limit, steps, particles, particle_coordinates_acw, galaxy1_coordinates_acw, particle_distances_acw)


#plot results
plt.figure(figsize=(20, 10))
plt.title("Tail Size Vs. Time - Closest Approach = 8, Mass Ratio = 2", fontsize = 20)
plt.plot(time_list_cw[:steps-2], tail_size_cw, label = "Clcokwise Interaction")
plt.plot(time_list_acw[:steps-2], tail_size_acw, label = "Anticlockwise Interaction")
plt.xlabel("Time Since Closest Approach", fontsize = 16)
plt.ylabel("Fraction of Total Stars in Tail", fontsize = 16)
plt.legend(fontsize = 14)
plt.savefig("TailSizeVsTime.png")
plt.show()
