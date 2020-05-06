# ---- CODE THAT PLOTS THE EVOLUTION OF INTERACTING GALAXIES OVER TIME IN, CLOCKWISE AND ANTICLOCKWISE -------------------

## Run Simulation, Plot results at some points along the evolution of the tail
# Initial Conditions
radii = np.zeros(50)
densities = np.zeros(int(len(radii)))

for i in range(int(len(radii))):
    radii[i] = 2 + (5/int(len(radii)))*i
    densities[i] = 20 * radii[i]

densities = densities.astype(int)




time_before_collision = 100
separation = 8
mass_ratio = 2
steps = 1000
particles = sum(densities)


# Simulate
# CLOCKWISE
clockwise = 1
particle_coordinates_cw, galaxy1_coordinates_cw, galaxy2_coordinates_cw, time_list_cw, particle_distances_cw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)

#ANTICLCKWISE
clockwise = -1
particle_coordinates_acw, galaxy1_coordinates_acw, galaxy2_coordinates_acw, time_list_acw, particle_distances_acw = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)


# Plot results
#CLOCKWISE
plt.figure(figsize=(20, 10))
plt.plot(galaxy1_coordinates_cw[:, 0], galaxy1_coordinates_cw[:, 1], linestyle = '--', label='Galaxy 1 Trajectory')
plt.plot(galaxy2_coordinates_cw[:, 0], galaxy2_coordinates_cw[:, 1], linestyle = '--', label='Galaxy 2 Coordinates')

def plot_cw(step):
    time = time_list_cw[step]
    plot_label = "Time = " + str(int(time))
    plt.scatter(particle_coordinates_cw[step, :, 0], particle_coordinates_cw[step, :, 1], s = 5, label = plot_label)

plot_cw(0)
plot_cw(200)
plot_cw(500)
plot_cw(800)

plt.title("The Evolution of a Tidal Tail Over Time (Clockwise Interaction) - Separation = 8, Mass Ratio = 2", fontsize = 20)
plt.xlabel("x-position", fontsize = 16)
plt.ylabel("y-position", fontsize = 16)
plt.axis([-100, 90, -35, 60])
plt.legend(fontsize = 16)
plt.savefig("ClockwiseExamplePlot.png")
plt.show()


# ANTICLOCKWISE
plt.figure(figsize=(20, 10))
plt.plot(galaxy1_coordinates_acw[:, 0], galaxy1_coordinates_acw[:, 1], linestyle = '--', label='Galaxy 1 Trajectory')
plt.plot(galaxy2_coordinates_acw[:, 0], galaxy2_coordinates_acw[:, 1], linestyle = '--', label='Galaxy 2 Coordinates')

def plot_acw(step):
    time = time_list_acw[step]
    plot_label = "Time = " + str(int(time))
    plt.scatter(particle_coordinates_acw[step, :, 0], particle_coordinates_acw[step, :, 1], s = 5, label = plot_label)

plot_acw(0)
plot_acw(200)
plot_acw(500)
plot_acw(800)


plt.title("The Evolution of a Tidal Tail Over Time (Anticlockwise Interaction) - Separation = 8, Mass Ratio = 2", fontsize = 20)
plt.xlabel("x-position", fontsize = 16)
plt.ylabel("y-position", fontsize = 16)
plt.axis([-25, 165, -35, 150])
plt.legend(fontsize = 16)
plt.savefig("AntiClockwiseExamplePlot.png")
plt.show()
