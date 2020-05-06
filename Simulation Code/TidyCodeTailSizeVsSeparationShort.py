##### Tail Fraction vs Closest Approach -----------
# Tail Size Functions 
def find_tail_size_cw(inner_limit, steps, particles, particle_coordinates, galaxy_coordinates, dr):
    tail_size = np.zeros(steps-2)
    for i in range(steps - 2):
        for j in range(particles):
            if dr[i, j] > inner_limit: # Star must be outside a certain inner radial limit from galaxy
                if particle_coordinates[i, j, 1]  > galaxy_coordinates[i, 1]: # Star must have y-coordinate greater than galaxy
                    tail_size[i] += 1
        tail_size[i] = tail_size[i]/particles
    return tail_size

def find_tail_size_acw(inner_limit, steps, particles, particle_coordinates, galaxy_coordinates, dr):
    tail_size = np.zeros(steps-2)
    for i in range(steps - 2):
        for j in range(particles):
            if dr[i, j] > inner_limit: # Star must be outside a certain inner radial limit from galaxy
                if particle_coordinates[i, j, 0]  > galaxy_coordinates[i, 0]: # Star must have x-coordinate greater than galaxy
                    tail_size[i] += 1
        tail_size[i] = tail_size[i]/particles
    return tail_size

# Initial conditions

radii = np.zeros(8)
densities = np.zeros(int(len(radii)))
densities = densities.astype(int)

for i in range(int(len(radii))):
    radii[i] = 2 + (5/int(len(radii)))*i
    densities[i] = radii[i]*6

particles = sum(densities)
time_before_collision = 20
mass_ratio = 1
steps = 2000

inner_limit = 10


# Produce a range of closest approach to test
list_of_separations = np.linspace(6, 15, 20)


# Repeat for error analysis
repeats = 10


# Tail Size Analysis
# CLOCKWISE
clockwise = 1
final_tail_size_cw = np.zeros([repeats, int(len(list_of_separations))])
for m in range(repeats):
    for n in range(int(len(list_of_separations))):
        separation = list_of_separations[n]
        particle_coordinates_cw, galaxy1_coordinates_cw, galaxy2_coordinates_cw, time_list_cw, particle_distances_cw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
        tail_size_cw = find_tail_size_cw(inner_limit, steps, particles, particle_coordinates_cw, galaxy1_coordinates_cw, particle_distances_cw)
        final_tail_size_cw[m, n] = tail_size_cw[steps-3]

end = time.time()


# ANTICLOCKWISE
clockwise = -1
final_tail_size_acw = np.zeros([repeats, int(len(list_of_separations))])
for m in range(repeats):
    for n in range(int(len(list_of_separations))):
        separation = list_of_separations[n]
        particle_coordinates_acw, galaxy1_coordinates_acw, galaxy2_coordinates_acw, time_list_acw, particle_distances_acw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
        tail_size_acw = find_tail_size_acw(inner_limit, steps, particles, particle_coordinates_acw, galaxy1_coordinates_acw, particle_distances_acw)
        final_tail_size_acw[m, n] = tail_size_acw[steps-3]

end = time.time()

# Produce average results and standard error
# Clcokwise
avg_cw = np.zeros(int(len(list_of_separations)))
error_cw = np.zeros(int(len(list_of_separations)))

for i in range((int(len(list_of_separations)))):
    avg_cw[i] = np.mean(final_tail_size_cw[:, i])
    error_cw[i] = np.std(final_tail_size_cw[:, i])

# Anticlockwise
avg_acw = np.zeros(int(len(list_of_separations)))
error_acw = np.zeros(int(len(list_of_separations)))

for i in range((int(len(list_of_separations)))):
    avg_acw[i] = np.mean(final_tail_size_acw[:, i])
    error_acw[i] = np.std(final_tail_size_acw[:, i])


# Plot Results
plt.figure(figsize=(20, 10))

plt.errorbar(list_of_separations, avg_cw, yerr=error_cw, xerr=None, fmt='none', label="Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.errorbar(list_of_separations, avg_acw, yerr=error_acw, xerr=None, fmt='none', label="Antilockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)

plt.title("Fraction of Stars in Tail at Different Closest Approaches - Mass Ratio = 1", fontsize=20)
plt.xlabel("Separation at Closest Approach", fontsize = 16)
plt.ylabel("Fraction of Stars in Tail", fontsize = 16)
plt.legend(fontsize=14)
plt.savefig("TailSizeVsSeparation.png")
plt.show()
