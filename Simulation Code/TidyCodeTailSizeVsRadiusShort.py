## -------------- FRACTION OF STARS THAT FORM TAIL AT DIFFERENT ORBITAL RADII ---------------
# define function that counts the fraction of stars in the tail for each radius

# CLOCKWISE
def radius_tail_fraction_cw(radii, densities, particle_coordinates, galaxy_coordinates, dr, steps):
    tail_fraction = np.zeros(int(len(radii)))
    index = 0
    for i in range(int(len(radii))):
        for j in range(densities[i]):
            if dr[steps-3, index] > inner_limit:
                if particle_coordinates[steps-3, index, 1]  > galaxy_coordinates[steps-3, 1]:
                    tail_fraction[i] += 1
                    index += 1
                else:
                    index += 1
            else:
                index += 1

        
        tail_fraction[i] = tail_fraction[i]/densities[i]

            
    return tail_fraction


#ANTICLOCKWISE
def radius_tail_fraction_acw(radii, densities, particle_coordinates, galaxy_coordinates, dr, steps):
    tail_fraction = np.zeros(int(len(radii)))
    index = 0
    for i in range(int(len(radii))):
        for j in range(densities[i]):
            if dr[steps-3, index] > inner_limit:
                if particle_coordinates[steps-3, index, 0]  > galaxy_coordinates[steps-3, 0]:
                    tail_fraction[i] += 1
                    index += 1
                else:
                    index += 1
            else:
                index += 1

        
        tail_fraction[i] = tail_fraction[i]/densities[i]

            
    return tail_fraction




# Initial Conditions
radii = np.zeros(50)
densities = np.zeros(int(len(radii)))

for i in range(int(len(radii))):
    radii[i] = 2 + (5*i)/int(len(radii))
    densities[i] = 15*radii[i]

densities = densities.astype(int)


time_before_collision = 100
separation = 8
mass_ratio = 2
steps = 3000
particles = sum(densities)
print(particles)
# define an inner limit
inner_limit = 10

# Record a few repeats of the simulations to obtain accurate results w/ errorbars
repeats = 5

tail_fraction_cw = np.zeros([repeats, int(len(radii))])
tail_fraction_acw = np.zeros([repeats, int(len(radii))])

for j in range(repeats):

    # Simulate CLOCKWISE
    clockwise = 1
    particle_coordinates_cw, galaxy1_coordinates_cw, galaxy2_coordinates_cw, time_list_cw, particle_distances_cw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)

    # find size of tail over tail
    tail_fraction_cw[j] = radius_tail_fraction_cw(radii, densities, particle_coordinates_cw, galaxy1_coordinates_cw, particle_distances_cw, steps)


    # Simulate ANTICLOCKWISE
    clockwise = -1
    particle_coordinates_acw, galaxy1_coordinates_acw, galaxy2_coordinates_acw, time_list_acw, particle_distances_acw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)

    # find size of tail over tail
    tail_fraction_acw[j] = radius_tail_fraction_acw(radii, densities, particle_coordinates_acw, galaxy1_coordinates_acw, particle_distances_acw, steps)


average_tail_fraction_cw = np.zeros(int(len(radii)))
average_tail_fraction_acw = np.zeros(int(len(radii)))

error_tail_fraction_cw = np.zeros(int(len(radii)))
error_tail_fraction_acw = np.zeros(int(len(radii)))


for i in range(int(len(radii))):
    average_tail_fraction_cw[i] = np.mean(tail_fraction_cw[:, i])
    average_tail_fraction_acw[i] = np.mean(tail_fraction_acw[:, i])

    error_tail_fraction_cw[i] = np.std(tail_fraction_cw[:, i])
    error_tail_fraction_acw[i] = np.std(tail_fraction_acw[:, i])


#plot results
plt.figure(figsize=(20, 10))
plt.errorbar(radii, average_tail_fraction_cw, error_tail_fraction_cw, fmt='none', label = "Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.errorbar(radii, average_tail_fraction_acw, error_tail_fraction_acw, fmt='none', label = "Anticlockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)

plt.xlabel("Initial Orbital Radius", fontsize = 16)
plt.ylabel("Fraction of Stars at given Radius in Tail", fontsize = 16)
plt.title("Fraction of Stars that form Tail at a Given Radius Vs. Radius - Closest Approach = 8, Mass Ratio = 2", fontsize=20)
plt.legend(fontsize = 14)
plt.savefig("RadialDistribution.png")


plt.show()
