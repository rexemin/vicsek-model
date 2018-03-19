"""
This module implements the Vicsek simulation.

Author: Ivan A. Moreno Soto
Last updated: 16/March/2018

TODO:
1. Make animation of the system evolution.
"""

################################################################################

from math import sqrt, cos, sin, pi
from random import random, uniform

# For plotting.
from matplotlib import pyplot as plt

# For writing a csv file
import csv

################################################################################

class Particle:
    """
    Defines a Particle object used in the Vicsek simulation.
    """
    def __init__(self, iniPos, iniVel, iniHead):
        """
        Defines the starting parameters of a Particle.

        @param iniPos: Vector of initial position.
        @param iniVel: Vector of initial velocity.
        @param iniHead: Angle of initial heading.
        """
        self.position = iniPos
        self.velocity = iniVel
        self.heading = iniHead

################################################################################

def dist(point1, point2, l):
    """
    Computes the distance between two points in the torus.

    @param point1: Point where the difference vector begins.
    @param point2: Point where the difference vector ends.
    @param l: Length of the square's side.
    """

    # First, we compute the difference vector.
    diff = [point1[0] - point2[0], point1[1] - point2[1]]

    # Now, we compute the shortest distance, assuming that
    # opposite sides are identified.
    for dim in range(2):
        while diff[dim] > l:
            diff[dim] -= l
        while diff[dim] < -l:
            diff[dim] += l

    return sqrt( (diff[0])**2 + (diff[1])**2 )

################################################################################

def nextHeading(particle, system, r, p_int, l):
    """
    Computes the next heading a particle will have on time t + step.

    @param particle: Particle whose heading will be computed.
    @param system: List of Particles in the system.
    @param r: Radius of interaction.
    @param p_int: Interval for the noise.
    @param l: Length of the square's side.
    """

    # We get the number of neighbors.
    N = sum([1 for neig in system if dist(particle.position, neig.position, l) <= r])

    avg_vel = ( sum([neig.heading for neig in system
                     if dist(particle.position, neig.position, l) <= r])
              ) / N

    return avg_vel + uniform(-p_int/2, p_int/2)

################################################################################

def nextPosition(particle, nextVelocity, step):
    """
    Computes the position a particle will have on time t + step.

    @param particle: Particle whose position will be computed.
    @param nextVelocity: Velocity for the moment t + step.
    @param step: Difference from moment to moment.
    """
    return [particle.position[0] + (nextVelocity[0] * step),
            particle.position[1] + (nextVelocity[1] * step)]

################################################################################

def makeHeaders(headers, file_path):
    """
    Prints the headers of a csv file.

    @param headers: Headers of the file (as an iterable object).
    @param file_path: Path of the file to be written.
    """
    with open(file_path, 'w', newline='') as db:
        writer = csv.writer(db)
        writer.writerow(headers)

################################################################################

def printSystem(system, t, file_path = ''):
    """
    Prints the position, velocity, and heading of every particle in
    a system.

    @param: System to be printed.
    @param t: Current time of the system's state.
    @param file_path: File where the system will be printed.
    """
    with open(file_path, 'a', newline='') as db:
        writer = csv.writer(db)

        for p in system:
            writer.writerow([t,
                             p.position[0],
                             p.position[1],
                             p.velocity[0],
                             p.velocity[1]])

################################################################################

def plotOrderParameter(order, iterations, step):
    """
    Plots the history of values of an Order parameter.

    @param order: List containing the history of values of the Order parameter.
    @param iterations: Length of order.
    @param step: Step of time from iteration to iteration.
    """

    # We create the figure.
    plt.ioff()
    fig = plt.figure()
    ax = plt.axes()

    # We get the times when the parameter was computed.
    x = [t for t in range(0, int(iterations*step), int(step))]

    ax.plot(x, order)
    plt.show()

################################################################################

def plot_grid(ax, system):
    """
    Plots the current state of a system with the respective directions
    of movement of every particle.

    @param ax: Container for the plot.
    @param system: System of particles.
    """
    plotx = [p.position[0] for p in system]
    ploty = [p.position[1] for p in system]

    plotu = [p.velocity[0] for p in system]
    plotv = [p.velocity[1] for p in system]

    return ax.quiver(plotx, ploty, plotu, plotv)

################################################################################

def vicsek(N, l, r, iterations, p_int, step, v_0, file_path):
    """
    Simulates a flocking system using the Vicsek model.
    Returns a history of the order parameter value for each moment
    in time.

    @param N: Number of particles to be created in the system.
    @param l: Length of the square's side.
    @param r: Radius of interaction.
    @param iterations: Number of iterations of the system.
    @param p_int: Perturbance interval for the noise.
    @param step: Step of time from iteration to iteration.
    @param v_0: Velocity of self propelled particles.
    """
    # We initialize the order parameter history and variables
    # needed for plotting.
    order_history = [0]

    makeHeaders(['t', 'px', 'py', 'vx', 'vy'], file_path)

    # Make random starting position.
    angles = [uniform(-2*pi, 2*pi) for part in range(N)]
    velos_x = [cos(part)*v_0 for part in angles]
    velos_y = [sin(part)*v_0 for part in angles]

    # We create the system of particles.
    system = [Particle([uniform(0.0, l), uniform(0.0, l)],
                       [velos_x[p], velos_y[p]],
                        angles[p]) for p in range(N)]

    printSystem(system, 0, file_path)

    # Evolve the system.
    for iteration in range(1, iterations):
        # We sequentially compute the headings, velocities, and positions
        # for every particle at the same time.
        new_headings = [nextHeading(particle, system, r, p_int, l)
                        for particle in system]

        new_velocities = [(v_0 * cos(head),
                           v_0 * sin(head)) for head in new_headings]

        new_pos = [nextPosition(particle, v, step)
                   for particle, v in zip(system, new_velocities)]

        # Changes all the particles attributes.
        for (particle, pos, h, v) in zip(system, new_pos, new_headings, new_velocities):
            particle.heading = h

            particle.velocity[0] = v[0]
            particle.velocity[1] = v[1]

            particle.position[0] = pos[0]
            particle.position[1] = pos[1]

            # Checks if the particle loops around.
            if particle.position[0] >= l:
                particle.position[0] -= l
            elif particle.position[0] <= 0:
                particle.position[0] += l

            if particle.position[1] >= l:
                particle.position[1] -= l
            elif particle.position[1] <= 0:
                particle.position[1] += l

        # Calculate the order parameter for this time step.
        order_x = abs(sum([particle.velocity[0] for particle in system])) / (N * v_0)
        order_y = abs(sum([particle.velocity[1] for particle in system])) / (N * v_0)
        order = sqrt(order_x**2 + order_y**2)

        order_history.append(order)

        print("Order parameter for time " + str(iteration*step) + ": " + str(order))
        printSystem(system, iteration*step, file_path)

    return order_history

################################################################################

# Test of the Vicsek model defined here.
N = 300
l = 25
r = 1.0
steps = 1000
p_int = 0.1
dt = 1.0
v_0 = 0.03
order_history = vicsek(N, l, r, steps, p_int, dt, v_0, 'sim.csv')

plotOrderParameter(order_history, steps, dt)
