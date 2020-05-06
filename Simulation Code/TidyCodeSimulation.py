
## ------------------------------------SIMULATION CODE ----------------------------
## ---------------------------------------------------------------------------------
## This code produces on simulation given some input initial conditions



# Variable Timestep Included
# Verlet Program that uses variable timesteps to improve efficiency
# Verlet integrate galaxies

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.integrate import odeint
import time
from matplotlib.colors import LogNorm

# Derivatives that describe galaxy motion
def derivatives(variables, t, mass_ratio):
    x1, y1, vx1, vy1, x2, y2, vx2, vy2 = variables
    mr = mass_ratio
    dx = x2 - x1
    dy = y2 - y1
    drsq = (dx*dx) + (dy*dy)
    dr = np.sqrt(drsq)
    rcubed = dr*dr*dr 
    dydt = [vx1, vy1, mr*dx/rcubed, mr*dy/rcubed, vx2, vy2, -dx/rcubed, -dy/rcubed]
    return dydt

# Function that creates parabolic initial conditions
def parabolic_velocities(separation, mass_ratio):
    v2sq = 2/(separation*(mass_ratio + 1))
    v2 = np.sqrt(v2sq)
    v1 = mass_ratio * v2
    velocities = [-v1, v2]
    return velocities

# function that simulates the galaxies up to a specified time with a suitable timestep (1/10 000 of the time)
# The function should return the final coordinates with reversed velocities, which gives the new initial conditions for the parabolic orobit

def find_galaxy_initials(separation, mass_ratio, time_before_collision):
    timestep = time_before_collision/10000
    steps = int(time_before_collision/timestep)
    t = np.linspace(0, steps*timestep, steps)
    
    velocities = parabolic_velocities(separation, mass_ratio)

    # initial positions at closest approach
    initials = [-separation/2, 0, 0, velocities[0], separation/2, 0, 0, velocities[1]]

    # Generate solution
    sol = odeint(derivatives, initials, t, args=(mass_ratio,))

    #Find final coordinates
    coordinates = sol[steps - 1, :]

    # Reverse Velocities
    # Galaxy 1
    coordinates[2] = -coordinates[2]
    coordinates[3] = -coordinates[3]
    # Galaxy 2
    coordinates[6] = -coordinates[6]
    coordinates[7] = -coordinates[7]

    return coordinates



# Function that generates test particles around a central mass

def generate_particles(radii, densities, mass, clockwise):
    particles = sum(densities)
    particle_coordinates = np.zeros([particles, 4])
    number_of_radii = int(len(radii))

    index = 0
    for i in range(number_of_radii):
        number_at_radius = densities[i]
        theta_interval = 2*np.pi/number_at_radius
        R = radii[i]
        random_angle = 2* np.pi * np.random.rand() # eliminates possible artefacts from particles being alligned from theta = 0
        V = np.sqrt((mass/R))
        for j in range(number_at_radius):
            angle = random_angle + j*theta_interval
            particle_x = R*np.sin(angle)
            particle_y = R*np.cos(angle)
            particle_vx = clockwise*V*np.cos(angle)
            particle_vy = -clockwise*V*np.sin(angle)
            particle_coords = [particle_x, particle_y, particle_vx, particle_vy]
            particle_coordinates[index, :] = particle_coords

            index += 1

    return particle_coordinates



# Translate particle position and velocity to sit around local galaxy at initial conditions

def generate_translated_particles(radii, densities, separation, mass_ratio, time_before_collision, clockwise):
    central_coordinates = generate_particles(radii, densities, 1, clockwise)
    galaxy_coordinates = find_galaxy_initials(separation, mass_ratio, time_before_collision)

    particles = sum(densities)
    new_particle_initials = np.zeros([particles, 4])

    for i in range(particles):
        for j in range(4):
            new_particle_initials[i, j] = central_coordinates[i, j] + galaxy_coordinates[j]

    return new_particle_initials


# Acceleration of a particle given coordinates of galaxies and their mass ratios

def acceleration_particle(xp, galaxy1, galaxy2, mass_ratio):
    x_p = xp[0]
    y_p = xp[1]
    x1 = galaxy1[0]
    y1 = galaxy1[1]
    x2 = galaxy2[0]
    y2 = galaxy2[1]

    dx1 = x1 - x_p
    dy1 = y1 - y_p
    dr1sq = (dx1*dx1) + (dy1*dy1)
    dr1 = np.sqrt(dr1sq)
    dr1cubed = dr1*dr1*dr1
    
    dx2 = x2 - x_p
    dy2 = y2 - y_p
    dr2sq = (dx2*dx2) + (dy2*dy2)
    dr2 = np.sqrt(dr2sq)
    dr2cubed = dr2*dr2*dr2
    
    ax1 = dx1/dr1cubed
    ay1 = dy1/dr1cubed

    ax2 = mass_ratio*dx2/dr2cubed
    ay2 = mass_ratio*dy2/dr2cubed

    ax = ax1+ax2
    ay = ay1+ay2

    acceleration = [ax, ay]

    return acceleration, dr1

# Acceleration of galaxies due to each other (ignore mass ratio factor for now)

def acceleration_galaxy(galaxy1, galaxy2):
    x1 = galaxy1[0]
    y1 = galaxy1[1]
    x2 = galaxy2[0]
    y2 = galaxy2[1]

    dx = x2 - x1
    dy = y2 - y1
    drsq = (dx*dx) + (dy*dy)
    dr = np.sqrt(drsq)

    acceleration = np.zeros(2)
    acceleration[0] = dx/(dr*dr*dr)
    acceleration[1] = dy/(dr*dr*dr)

    return acceleration


# Function that returns particle and galaxy coordinates at every step with cooresponding time and distances of particls
# from the local galaxy - initial conditions are of two galaxies in parabolic orbits with a circular collection of stars around one of the galaxies.

# Input parameters for the orbit are the separation of the galaxy centres at closest approach, the mass ratio of the galaxies, the particular radial
# distribution of the stars in the galaxy.
#   whether the stars are orbitting clockwise or anticlockwise (prograde/retrograde), and the time before the closest approach at initialisation.

def simulate_galaxy(radii, densities, time_before_collision, separation, mass_ratio, clockwise, steps):
    min_timestep = 0.1
    theoretical_maximum_acceleration =0.5*((-0.01417*separation*mass_ratio) + (0.237*mass_ratio) + (0.00685*separation) + 0.11)
    time_list = np.zeros(steps)
    time_list[0] = -time_before_collision

    # Initializw coordinates

    particle_initial_coordinates = generate_translated_particles(radii, densities, separation, mass_ratio, time_before_collision, clockwise)
    galaxy_coordinates = find_galaxy_initials(separation, mass_ratio, time_before_collision)
    galaxy1_initial_coordinates = galaxy_coordinates[0:4]
    galaxy2_initial_coordinates = galaxy_coordinates[4:8]
    
    particle_coordinates = np.zeros([steps, particles, 2])
    particle_distances = np.zeros([steps, particles])
    galaxy1_coordinates = np.zeros([steps, 2])
    galaxy2_coordinates = np.zeros([steps, 2])


    # initial positions
    
    for i in range(particles):
        particle_coordinates[0, i, 0] = particle_initial_coordinates[i, 0]
        particle_coordinates[0, i, 1] = particle_initial_coordinates[i, 1]
    
    galaxy1_coordinates[0, 0] = galaxy1_initial_coordinates[0]
    galaxy1_coordinates[0, 1] = galaxy1_initial_coordinates[1]
    
    galaxy2_coordinates[0, 0] = galaxy2_initial_coordinates[0]
    galaxy2_coordinates[0, 1] = galaxy2_initial_coordinates[1]



    # initial velocities

    particle_initial_velocities = np.zeros([particles, 2])
    galaxy1_initial_velocity = np.zeros(2)
    galaxy2_initial_velocity = np.zeros(2)
    
    for i in range(particles):
        particle_initial_velocities[i, 0] = particle_initial_coordinates[i, 2]
        particle_initial_velocities[i, 1] = particle_initial_coordinates[i, 3]
        
    galaxy1_initial_velocity[0] = galaxy1_initial_coordinates[2]
    galaxy1_initial_velocity[1] = galaxy1_initial_coordinates[3]

    galaxy2_initial_velocity[0] = galaxy2_initial_coordinates[2]
    galaxy2_initial_velocity[1] = galaxy2_initial_coordinates[3]
    
    

    ## calculate initial accelerations
    # initialise vectors
    particle_accelerations = np.zeros([particles, 2])
    galaxy1_acceleration = np.zeros(2)
    galaxy2_acceleration = np.zeros(2)

    # assign galaxy positions and accelerations
    galaxy1_position = galaxy1_coordinates[0]
    galaxy2_position = galaxy2_coordinates[0]

    list_of_accelerations = np.zeros(2*particles)

    # assign particle positions and accelerations, record in list to check for maximum values later
    for i in range(particles):
        particle_position = particle_coordinates[0, i]
        particle_accelerations[i], particle_distances[0, i] = acceleration_particle(particle_position, galaxy1_position, galaxy2_position, mass_ratio)

        list_of_accelerations[(i*2)] = abs(particle_accelerations[i, 0])
        list_of_accelerations[(i*2) + 1] = abs(particle_accelerations[i, 1])

    # assign current galaxy acceleration value (without factor of mass)
    current_acceleration = acceleration_galaxy(galaxy1_position, galaxy2_position)

    # acceleration of each galaxy
    galaxy1_acceleration = mass_ratio*current_acceleration
    galaxy2_acceleration = -current_acceleration


    # Calculate timestep size using error analysis results, check for 95th percentile acceleration
    list_of_accelerations.sort()
    current_maximum_acceleration = list_of_accelerations[int(particles*1.9)]
    timestep_n = min_timestep + np.log(theoretical_maximum_acceleration/current_maximum_acceleration)

    if timestep_n < min_timestep:
        timestep_n = min_timestep # safety measure, if maximum goes above theoretical prediction, timestep does not go negative or ridiculous;y small that stops computation all together

        
    # Use verlet formula
    #x1 = x0 +v0*dt + 1/2*a0*dtt*dt

    # iterate particles
    for i in range(particles):
        particle_coordinates[1, i, 0] = particle_coordinates[0, i, 0] + (particle_initial_velocities[i, 0] * timestep_n) + (0.5*particle_accelerations[i, 0]*timestep_n*timestep_n)
        particle_coordinates[1, i, 1] = particle_coordinates[0, i, 1] + (particle_initial_velocities[i, 1] * timestep_n) + (0.5*particle_accelerations[i, 1]*timestep_n*timestep_n)

    # iterate galaxies
    galaxy1_coordinates[1, 0] = galaxy1_coordinates[0, 0] + (galaxy1_initial_velocity[0]*timestep_n) + (0.5*galaxy1_acceleration[0]*timestep_n*timestep_n)
    galaxy1_coordinates[1, 1] = galaxy1_coordinates[0, 1] + (galaxy1_initial_velocity[1]*timestep_n) + (0.5*galaxy1_acceleration[1]*timestep_n*timestep_n)

    galaxy2_coordinates[1, 0] = galaxy2_coordinates[0, 0] + (galaxy2_initial_velocity[0]*timestep_n) + (0.5*galaxy2_acceleration[0]*timestep_n*timestep_n)
    galaxy2_coordinates[1, 1] = galaxy2_coordinates[0, 1] + (galaxy2_initial_velocity[1]*timestep_n) + (0.5*galaxy2_acceleration[1]*timestep_n*timestep_n)

    # record timestep size
    time_list[1] = time_list[0] + timestep_n

    # reassign timestep_n
    timestep_n_minus_one = timestep_n

    # initialise list of times
    timestep_list = []



    ## Before full Verlet method, initialise important vectors

    #initialise current positions
    galaxy1_current_position = np.zeros(2)
    galaxy2_current_position = np.zeros(2)

    # Change since last step
    dg1 = np.zeros(2)
    dg2 = np.zeros(2)

    # changes in position for each particle
    dpos_list = np.zeros([particles, 2])

    # accelerations of each particle, in one list for ordering purposes
    list_of_accelerations = np.zeros(2*particles)
    
    ## Now implement full verlet integration method (n > 1)
    for i in range(1, steps - 1):

        ## and assign.....

        # galaxy positions
        galaxy1_current_position[0] = galaxy1_coordinates[i, 0]
        galaxy1_current_position[1] = galaxy1_coordinates[i, 1]

        
        galaxy2_current_position[0] = galaxy2_coordinates[i, 0]
        galaxy2_current_position[1] = galaxy2_coordinates[i, 1]
        
        #current galaxy accelerations
        current_acceleration = acceleration_galaxy(galaxy1_current_position, galaxy2_current_position)
        galaxy1_current_acceleration = mass_ratio*current_acceleration
        galaxy2_current_acceleration = -current_acceleration

        # changes in galaxy position
        dg1[0] = galaxy1_coordinates[i, 0] - galaxy1_coordinates[i-1, 0]
        dg1[1] = galaxy1_coordinates[i, 1] - galaxy1_coordinates[i-1, 1]

        dg2[0] = galaxy2_coordinates[i, 0] - galaxy2_coordinates[i-1, 0]
        dg2[1] = galaxy2_coordinates[i, 1] - galaxy2_coordinates[i-1, 1]

        # particle positions, accelerations, changes in position
        for j in range(particles):
            current_particle_position = particle_coordinates[i, j]
            # caluclate acceleration, also record distance from local galaxy centre (this was used in acceleration calculation and is analysed later)
            particle_acceleration, particle_distances[i, j] = acceleration_particle(current_particle_position, galaxy1_current_position, galaxy2_current_position, mass_ratio)
            dpos_list[j, 0] = particle_coordinates[i, j, 0] - particle_coordinates[i-1, j, 0]
            dpos_list[j, 1] = particle_coordinates[i, j, 1] - particle_coordinates[i-1, j, 1]

            # lists of accelerations
            list_of_accelerations[(j*2)] = (particle_acceleration[0])
            list_of_accelerations[(j*2) + 1] = (particle_acceleration[1])


        # Calculate timestep size
        ordered_list = abs(list_of_accelerations)
        ordered_list.sort() # order magnitude of accelerations, find 95th percentile
        current_maximum_acceleration = ordered_list[int(particles*1.9)]

        # find timestep using formula from errors
        timestep_n = min_timestep + np.log(theoretical_maximum_acceleration/current_maximum_acceleration)
        if timestep_n < min_timestep:
            timestep_n = min_timestep # minimum timestep safety feature to keep simulation running

        # calclulate ratio in current to previous timestep
        timestep_ratio = timestep_n/timestep_n_minus_one

        # record size of timestep
        timestep_list.append(timestep_n)
        
        # iterate galaxy
        galaxy1_coordinates[i+1, 0] = galaxy1_current_position[0] + (dg1[0]*timestep_ratio) + (galaxy1_current_acceleration[0] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)
        galaxy1_coordinates[i+1, 1] = galaxy1_current_position[1] + (dg1[1]*timestep_ratio) + (galaxy1_current_acceleration[1] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)

                
        galaxy2_coordinates[i+1, 0] = galaxy2_current_position[0] + (dg2[0]*timestep_ratio) + (galaxy2_current_acceleration[0] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)
        galaxy2_coordinates[i+1, 1] = galaxy2_current_position[1] + (dg2[1]*timestep_ratio) + (galaxy2_current_acceleration[1] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)

        # iterate particles
        for j in range(particles):
            dpos = dpos_list[j]
            current_particle_position = particle_coordinates[i, j]
            particle_acceleration[0] = list_of_accelerations[(j*2)]
            particle_acceleration[1] = list_of_accelerations[(j*2) + 1]
            dpos[0] = particle_coordinates[i, j, 0] - particle_coordinates[i-1, j, 0]
            dpos[1] = particle_coordinates[i, j, 1] - particle_coordinates[i-1, j, 1]
            particle_coordinates[i+1, j, 0] = current_particle_position[0] + (dpos[0]*timestep_ratio) + (particle_acceleration[0] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)
            particle_coordinates[i+1, j, 1] = current_particle_position[1] + (dpos[1]*timestep_ratio) + (particle_acceleration[1] * (timestep_n + timestep_n_minus_one) * 0.5 * timestep_n)

        #iterate timestep
        time_list[i+1] = time_list[i] + timestep_n

        # reassign previous timestep
        timestep_n_minus_one = timestep_n

    return particle_coordinates, galaxy1_coordinates, galaxy2_coordinates, time_list, particle_distances
