## CRITCAL MASS CALCULATIONS AT DIFFERENT CLOSEST APPROACHES

radii = [2, 7]
densities = [2, 7]

time_before_collision = 5
particles = sum(densities)
steps = 501
mass_ratio = 1
inner_limit = 10
repeats = 5
critical_mass_repeats = 5



list_of_separations = [7, 8, 9, 10, 11, 12, 13, 14, 15]
critical_mass_2_cw = np.zeros([int(len(list_of_separations)), critical_mass_repeats])
critical_mass_7_cw = np.zeros([int(len(list_of_separations)), critical_mass_repeats])

critical_mass_2_acw = np.zeros([int(len(list_of_separations)), critical_mass_repeats])
critical_mass_7_acw = np.zeros([int(len(list_of_separations)), critical_mass_repeats])


start = time.time()


# CLOCKWISE
clockwise = 1

list_of_mass_ratios = np.logspace(-1.2, 1.7, 10)

for m in range(critical_mass_repeats):
    for n in range(int(len(list_of_separations))):
        separation = list_of_separations[n]

        extreme_fraction_cw = np.zeros([repeats, int(len(list_of_mass_ratios)), 2])

        average_2_cw = np.zeros(int(len(list_of_mass_ratios)))
        average_7_cw = np.zeros(int(len(list_of_mass_ratios)))

        for j in range(repeats):
            for i in range(int(len(list_of_mass_ratios))):
                mass_ratio = list_of_mass_ratios[i]
                particle_coordinates_cw, galaxy1_coordinates_cw, galaxy2_coordinates_cw, time_list_cw, particle_distances_cw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
                extreme_fraction_cw[j, i] = extreme_radius_tail_fraction_cw(radii, densities, particle_coordinates_cw, galaxy1_coordinates_cw, particle_distances_cw, steps)
                
                
        for i in range(int(len(list_of_mass_ratios))):
            average_2_cw[i] = np.mean(extreme_fraction_cw[:, i, 0])
            average_7_cw[i] = np.mean(extreme_fraction_cw[:, i, 1])
                                   


        count = 0

        for i in range(int(len(average_2_cw))):
            if count < 1:
                if average_2_cw[i] > 0.1:
                    count += 1
                    y1 = average_2_cw[i-1]
                    x1 = np.log10(list_of_mass_ratios[i-1])

                    y2 = average_2_cw[i]
                    x2 = np.log10(list_of_mass_ratios[i])

                    dy = y2 - y1
                    dx = x2 - x1
                    
                    xcrit = x1 + (dx/dy)*(0.1 - y1)

                    critical_mass_2_cw[n, m] = 10**xcrit

        count = 0
        for i in range(int(len(average_7_cw))):
            if count < 1:
                if average_7_cw[i] > 0.1:
                    count += 1
                    y1 = average_7_cw[i-1]
                    x1 = np.log10(list_of_mass_ratios[i-1])

                    y2 = average_7_cw[i]
                    x2 = np.log10(list_of_mass_ratios[i])

                    dy = y2 - y1
                    dx = x2 - x1
                    
                    xcrit = x1 + (dx/dy)*(0.1 - y1)
                
                    critical_mass_7_cw[n, m] = 10**xcrit



        

average_critical_mass_2_cw = np.zeros(int(len(list_of_separations)))
average_critical_mass_7_cw = np.zeros(int(len(list_of_separations)))

error_critical_mass_2_cw = np.zeros(int(len(list_of_separations)))
error_critical_mass_7_cw = np.zeros(int(len(list_of_separations)))

for i in range(int(len(list_of_separations))):
    average_critical_mass_2_cw[i] = np.mean(critical_mass_2_cw[i])
    average_critical_mass_7_cw[i] = np.mean(critical_mass_7_cw[i])

    error_critical_mass_2_cw[i] = np.std(critical_mass_2_cw[i])
    error_critical_mass_7_cw[i] = np.std(critical_mass_7_cw[i])


#ANTICLOCKWISE

clockwise = -1

list_of_mass_ratios = np.logspace(-1.2, 1.7, 10)

for m in range(critical_mass_repeats):
    for n in range(int(len(list_of_separations))):
        separation = list_of_separations[n]

        extreme_fraction_acw = np.zeros([repeats, int(len(list_of_mass_ratios)), 2])

        average_2_acw = np.zeros(int(len(list_of_mass_ratios)))
        average_7_acw = np.zeros(int(len(list_of_mass_ratios)))

        for j in range(repeats):
            for i in range(int(len(list_of_mass_ratios))):
                mass_ratio = list_of_mass_ratios[i]
                particle_coordinates_acw, galaxy1_coordinates_acw, galaxy2_coordinates_acw, time_list_acw, particle_distances_acw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
                extreme_fraction_acw[j, i] = extreme_radius_tail_fraction_acw(radii, densities, particle_coordinates_acw, galaxy1_coordinates_acw, particle_distances_acw, steps)
                
        for i in range(int(len(list_of_mass_ratios))):
            average_2_acw[i] = np.mean(extreme_fraction_acw[:, i, 0])
            average_7_acw[i] = np.mean(extreme_fraction_acw[:, i, 1])
                                   


        count = 0

        for i in range(int(len(average_2_acw))):
            if count < 1:
                if average_2_acw[i] > 0.1:
                    count += 1
                    y1 = average_2_acw[i-1]
                    x1 = np.log10(list_of_mass_ratios[i-1])

                    y2 = average_2_acw[i]
                    x2 = np.log10(list_of_mass_ratios[i])

                    dy = y2 - y1
                    dx = x2 - x1
                    
                    xcrit = x1 + (dx/dy)*(0.1 - y1)

                    critical_mass_2_acw[n, m] = 10**xcrit

        count = 0
        for i in range(int(len(average_7_acw))):
            if count < 1:
                if average_7_acw[i] > 0.1:
                    count += 1
                    y1 = average_7_acw[i-1]
                    x1 = np.log10(list_of_mass_ratios[i-1])

                    y2 = average_7_acw[i]
                    x2 = np.log10(list_of_mass_ratios[i])

                    dy = y2 - y1
                    dx = x2 - x1
                    
                    xcrit = x1 + (dx/dy)*(0.1 - y1)
                    
                    critical_mass_7_acw[n, m] = 10**xcrit


average_critical_mass_2_acw = np.zeros(int(len(list_of_separations)))
average_critical_mass_7_acw = np.zeros(int(len(list_of_separations)))

error_critical_mass_2_acw = np.zeros(int(len(list_of_separations)))
error_critical_mass_7_acw = np.zeros(int(len(list_of_separations)))




for i in range(int(len(list_of_separations))):
    average_critical_mass_2_acw[i] = np.mean(critical_mass_2_acw[i])
    average_critical_mass_7_acw[i] = np.mean(critical_mass_7_acw[i])

    error_critical_mass_2_acw[i] = np.std(critical_mass_2_acw[i])
    error_critical_mass_7_acw[i] = np.std(critical_mass_7_acw[i])

end = time.time()

print("Time Taken:")
print(end-start)


# Plot Results
plt.figure(figsize=(20, 10))
plt.errorbar(list_of_separations, average_critical_mass_2_cw, error_critical_mass_2_cw, fmt='none', capsize=5, elinewidth=2, markeredgewidth=2, label="Clockwise Interaction")
plt.errorbar(list_of_separations, average_critical_mass_2_acw, error_critical_mass_2_acw, fmt='none', capsize=5, elinewidth=2, markeredgewidth=2, label="Anticlockwise Interaction")
plt.xlabel("Separation of Galaxies at Closest Approach", fontsize = 16)
plt.ylabel("Critical Mass", fontsize = 16)
plt.title("'Critical Mass' Vs. Separation at Closest Approach - Initial Radius = 2", fontsize = 20)
plt.legend(fontsize=14)
plt.savefig("Radius2CriticalMassVsClosestApproach.png")
plt.show()


plt.figure(figsize=(20, 10))
plt.errorbar(list_of_separations, average_critical_mass_7_cw, error_critical_mass_7_cw, fmt='none', capsize=5, elinewidth=2, markeredgewidth=2, label="Clockwise Interaction")
plt.errorbar(list_of_separations, average_critical_mass_7_acw, error_critical_mass_7_acw, fmt='none', capsize=5, elinewidth=2, markeredgewidth=2, label="Anticlockwise Interaction")
plt.xlabel("Separation of Galaxies at Closest Approach", fontsize = 16)
plt.ylabel("Critical Mass", fontsize = 16)
plt.title("'Critical Mass' Vs. Separation at Closest Approach - Initial Radius = 7, Clockwise Interaction", fontsize = 20)
plt.legend(fontsize=14)
plt.savefig("Radius7CriticalMassVsClosestApproach.png")
plt.show()
