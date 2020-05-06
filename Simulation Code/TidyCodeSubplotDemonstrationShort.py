## -- CREATE SUBPLOTS TO ILLUSTRATE TWO GALAXY INTERACTION --------

radii=np.zeros(40)
densities=np.zeros(int(len(radii)))

for i in range(int(len(radii))):
    radii[i] = 2 + (5/int(len(radii)))*i
    densities[i] = 4 * radii[i]

densities= densities.astype(int)

time_before_collision = 50
separation = 8
mass_ratio = 1
clockwise = 1
steps = 801

particles = sum(densities)

start = time.time()
particle_coordinates, galaxy1_coordinates, galaxy2_coordinates, time_list, particle_distances = simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
particle_coordinates2, galaxy1_coordinates2, galaxy2_coordinates2, time_list2, particle_distances2 = simulate_galaxy2(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps)
finish = time.time()

time_taken = finish-start
print("Number of Particles:")
print(particles)
print()
print("Time:")
print(time_taken)
print()
print("Time Simulated:")
print(time_list[steps-1])


def plot(step):
    plt.plot(galaxy1_coordinates[:, 0], galaxy1_coordinates[:, 1])
    plt.plot(galaxy2_coordinates[:, 0], galaxy2_coordinates[:, 1])

    time = time_list[step]
    plot_time = "Time: " + str(int(time))
    plt.scatter(particle_coordinates[step, :, 0], particle_coordinates[step, :, 1], s =5, label= plot_time)
    plt.scatter(particle_coordinates2[step, :, 0], particle_coordinates2[step, :, 1], s = 5)
    plt.axis([-30, 30, -60, 60])
    plt.legend(loc="upper left", fontsize=14)
    plt.xlabel("x-position", fontsize=16)
    


plt.subplots(1, 5, sharey=True, figsize=(20, 10))
plt.suptitle("A Typical Interaction of Two Galaxies - Mass Ratio = 1, Closest Approach = 8, Both Rotating Clockwise", fontsize = 20)
plt.subplot(1, 5, 1)
plot(0)
plt.ylabel("y-position", fontsize=16)

plt.subplot(1, 5, 2)
plot(200)

plt.subplot(1, 5, 3)
plot(400)

plt.subplot(1, 5, 4)
plot(600)

plt.subplot(1, 5, 5)
plot(800)

# plt.savefig("TidyCodeSubplots.png")
plt.show()
