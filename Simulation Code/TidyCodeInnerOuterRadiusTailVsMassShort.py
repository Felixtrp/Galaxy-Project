## -------------- FRACTION OF STARS AT A SPECIFIC RADIUS IN TAIL FOR DIFFERENT MASSES ----------

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

# Only Consider innermost and outermost stars
radii = [2, 7]
densities = np.zeros(int(len(radii)))

for i in range(int(len(radii))):
    densities[i] = 7*radii[i]

densities = densities.astype(int)


time_before_collision = 100
separation = 8
steps = 3000
particles = sum(densities)
print(particles)
# define an inner limit
inner_limit = 10

list_of_masses = np.logspace(-2, 2, 20)

inner_fraction_cw = np.zeros(int(len(list_of_masses)))
inner_fraction_acw = np.zeros(int(len(list_of_masses)))

outer_fraction_cw = np.zeros(int(len(list_of_masses)))
outer_fraction_acw = np.zeros(int(len(list_of_masses)))

error_inner_fraction_cw = np.zeros(int(len(list_of_masses)))
error_inner_fraction_acw = np.zeros(int(len(list_of_masses)))

error_outer_fraction_cw = np.zeros(int(len(list_of_masses)))
error_outer_fraction_acw = np.zeros(int(len(list_of_masses)))







# Repeat at different masses
for n in range(int(len(list_of_masses))):
    mass_ratio = list_of_masses[n]

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


    
    error_tail_fraction_cw = np.zeros(int(len(radii)))
    error_tail_fraction_acw = np.zeros(int(len(radii)))


    inner_fraction_cw[n] = np.mean(tail_fraction_cw[:, 0])
    inner_fraction_acw[n] = np.mean(tail_fraction_acw[:, 0])

    outer_fraction_cw[n] = np.mean(tail_fraction_cw[:, 1])
    outer_fraction_acw[n] = np.mean(tail_fraction_acw[:, 1])

    error_inner_fraction_cw[n] = np.std(tail_fraction_cw[:, 0])
    error_inner_fraction_acw[n] = np.std(tail_fraction_cw[:, 0])

    error_outer_fraction_cw[n] = np.std(tail_fraction_cw[:, 1])
    error_outer_fraction_acw[n] = np.std(tail_fraction_cw[:, 1])


log_masses = np.log10(list_of_masses)


#plot results

#Inner Radius
plt.figure(figsize=(20, 10))
plt.errorbar(log_masses, inner_fraction_cw, error_inner_fraction_cw, fmt='none', label = "Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.errorbar(log_masses, inner_fraction_acw, error_inner_fraction_acw, fmt='none', label = "Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.xlabel("log(Mass Ratio)", fontsize = 16)
plt.ylabel("Fraction of Stars Initially at Radius = 2 in Tail", fontsize = 16)
plt.title("Fraction of Stars Initially at Radius = 2 that form Tail Vs. log(Mass Ratio) - Closest Approach = 8", fontsize=20)
plt.legend(loc="upper left", fontsize=14)
plt.axis([-2.2, 2.2, -0.1, 1.1])
plt.savefig("InnerRadiusTailSizeVsMass")
plt.show()




# Outer Radius
plt.figure(figsize=(20, 10))
plt.errorbar(log_masses, outer_fraction_cw, error_outer_fraction_cw, fmt='none', label = "Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.errorbar(log_masses, outer_fraction_acw, error_outer_fraction_acw, fmt='none', label = "Clockwise Interaction", capsize=5, elinewidth=2, markeredgewidth=2)
plt.xlabel("log(Mass Ratio)", fontsize = 16)
plt.ylabel("Fraction of Stars Initially at Radius = 7 in Tail", fontsize = 16)
plt.title("Fraction of Stars Initially at Radius = 7 that form Tail Vs. log(Mass Ratio) - Closest Approach = 8", fontsize=20)
plt.legend(loc="upper left", fontsize=14)
plt.axis([-2.2, 2.2, -0.1, 1.1])
plt.savefig("OuterRadiusTailSizeVsMass")
plt.show()
