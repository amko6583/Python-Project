
#get_ipython().magic(u'pylab inline')
import pylab
import numpy as np

import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D



G = 6.673e-11
sunList = []

with open('initial_Sun_amanda.txt','r') as FILE:
    lines = FILE.readlines()
    for line in lines:
        if line[0] == '#':
            continue
        sunList.append(line.split())
for i in range(len(sunList)):
    for j in range(7):
        sunList[i][j] = float(sunList[i][j])
    """
    This opens the file initial_Sun_amanda.txt
    It reads the lines, and if the line starts with #, it skips it and moves on to the next line.
    For each value in the line, it adds it to sunList when there is a space between values.
    The for loop then changes the strings within sunList to floats so they can be used in evolve particles
    """
        
massList = [sunList[0][0],sunList[1][0]]    #Mass Sun, Mass Earth
posList = [sunList[0][1:4],sunList[1][1:4]] #X,Y,Z position of Sun, Earth
velList = [sunList[0][4:],sunList[1][4:]]   #X,Y,Z velocity of Sun, Earth

print sunList
print massList
print posList
print velList



masses = np.array(massList)
positions = np.array(posList)*149597870700
velocities = np.array(velList)
print positions


def gForce(mass1, mass2, radius, epsilon = 0.0):
    """
    Calculate the gravitational force between two bodies.

    Input Parameters
    ----------
    mass1 : int, float
        Mass of the first body, in kilograms.
    mass2 : int, float
        Mass of the second body, in kilograms.
    radius : int, float
        Separation of the bodies, in meters.

    Returns
    -------
    Force in Newtons : float
    """
    force = G * mass1 * mass2 / (radius**2)  # Calculate the force
    
    return force


def sepVector(point1, point2):
    """
    Calculates a separation vector (as a list) between two particles.

    The returned vector points from particle p1 to particle p2.

    Parameters
    ----------
    point1 : list
        List containing [x, y, z] position of particle #1, which is
        the starting point of the vector.
    point2 : list
        List containing [x, y, z] position of particle #2, which is 
        the end point of the vector.

    Examples
    --------
    Vector from the Earth to the Sun:
    >>> vector = sepVector([0, 0, 0], [1.5e11,0,0])

    Returns
    -------
    Separation vector (3-element list). Units are whatever the positions were
    supplied in : list of ints or floats
    """
    #  Create the separation vector
    vector = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
    
    return vector


from math import sqrt

def magOfVector(vector):
    """
    Calculate the magnitude of the vector.
    
    Parameters
    ----------
    vector : list of ints or floats
        List containing [x, y, z] values

    Returns
    -------
    Magnitude of the vector : float
    """
    
    #  Calculate the magnitude of the given vector
    magnitude = sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)
    
    return magnitude



def gForceSum(particleList, massList, epsilon=0.0):
    """
    Computes the sum of the forces on each object from all other 
    objects.

    The returned list is comprised of sublists with the x, y, and z
    forces on each object.

    Parameters
    ----------
    particleList : 2D list of floats
        Each list contains mass, x,y,z positions, and vx, vy, vz velocities
    massList : 1D list of floats
        A list of masses in kg

    Returns
    -------
    List of forces with dimensions masses x 3 : list with lists of floats
    """
    #  Create lists to keep track of the total x, y, and z forces
    xForceTotal = [0] * len(particleList)
    yForceTotal = [0] * len(particleList)
    zForceTotal = [0] * len(particleList)
    
    #  Go through each particle (this gets us particle 1)
    for index1 in range(len(particleList)):
        #  Go through each particle again (this gets us particle 2)
        for index2 in range(len(particleList)):
            #  Don't compare the particle to itself
            if index1 != index2:
                #  Get particle 1 and particle 2 masses
                mass1 = massList[index1]
                mass2 = massList[index2]
                
                #  Get particle 1 and particle 2 coordinate points
                point1 = particleList[index1]
                point2 = particleList[index2]
                
                #  Calculate the separation vector between point 1 and 2
                sVector = sepVector(point1, point2)
                
                #  Calculate the radius between the 2 points
                radius = magOfVector(sVector)
                
                #  Calculate the force from the x, y, and z components and add them to the totals
                xForceTotal[index1] += gForce(mass1, mass2, radius,epsilon) * (sVector[0]) / np.sqrt(radius**2+epsilon**2)
                yForceTotal[index1] += gForce(mass1, mass2, radius,epsilon) * (sVector[1]) / np.sqrt(radius**2+epsilon**2)
                zForceTotal[index1] += gForce(mass1, mass2, radius,epsilon) * (sVector[2]) / np.sqrt(radius**2+epsilon**2)
    
    #  Convert the x, y, and z total force lists into a particle list with force totals
    #
    #  From:
    #     [[particle1 xForceTot, particle2 xForceTot, etc],
    #      [particle1 yForceTot, particle2 yForceTot, etc],
    #      [particle1 zForceTot, particle2 zForceTot, etc]]
    #
    #  To:
    #     [[particle1 xForceTot, particle1 yForceTot, particle1 zForceTot],
    #      [particle2 xForceTot, particle2 yForceTot, particle2 zForceTot],
    #      etc...]
    forceVectors = [(xForceTotal[index], yForceTotal[index], zForceTotal[index]) for index in range(len(particleList))]
    
    return forceVectors


def evolveParticles(positions, velocities, masses, deltaT, epsilon=0.0):
    """ 
    Evolve particles in time via leap-frog integrator scheme.
    
    Parameters
    ----------
    positions : np.ndarray
        2D array containing (x, y, z) positions for all particles. 
        Shape is (N, 3) where N is the number of particles.
    velocities : np.ndarray
        2D array containing (x, y, z) velocities for all particles. 
        Shape is (N, 3) where N is the number of particles.
    masses : np.ndarray
        1D array containing masses for all particles, length N, where
        N is the number of particles.
    deltaT : float
        Evolve system for time deltaT
    epsilon : float
        Softening length

    Returns
    -------
    Updated particle positions and particle velocities, each being a 2D
    array with shape (N, 3), where N is the number of particles.

    """ 
    
    # Make copies of position/velocity arrays that we can modify in-place.
    positions = positions.copy()
    velocities = velocities.copy()
    
    N = len(masses)               # Number of particles in system
    dims = positions.shape[-1]    # Dimensionality of problem
    
    # Compute net force vectors on all particles.
    # If you use your own code, change "gForceSum" and its 
    # input variables to match your code.
    #forces = gForceSum(positions, masses)
    forces = gForceSum(positions, masses, epsilon = epsilon)
    
    # 2 half steps for leap-frog method (first dimension)
    acc = np.zeros([2,N,dims])
    
    # Leap-frog integrator takes two half-steps
    for step in xrange(2):
        # Loop over particles, compute acceleration,
        # update positions and velocities
        for k in xrange(N):
            
            # Rec-calculate acceleration at each half-step
            acc[step,k] = forces[k] / masses[k]
            
            # Update position on first half-step, velocity on second
            if step == 0:
                positions[k,:] = positions[k] + velocities[k] * deltaT + 0.5 * acc[0,k] * deltaT**2
            else:
                velocities[k,:] = velocities[k] + 0.5 * (acc[0,k] + acc[1,k]) * deltaT
    
    return positions, velocities



Torbit = 2*np.pi*149597870700/29779.5
deltaT = Torbit/1e4
tf = Torbit
dt = deltaT
t1 = 0


def evolve(positions, velocities, masses, deltaT, tf, epsilon=0.0):
    i = 1
    t1 = 0
    #List of zeros, with the length of masses(number of particles), 3 for x, y and z positions, and the number of steps in the loop
    pos_list = np.zeros([len(masses), 3, (tf/deltaT)])
    #List of zeros, with the length of masses(number of particles), 3 for x, y and z velocities, and the number of steps in the loop
    vel_list = np.zeros([len(masses), 3, (tf/deltaT)])
    pos_list[:,:,0] = positions.copy() #Puts initial positions into pos_list
    vel_list[:,:,0] = velocities.copy() #Puts initial velocities into vel_list
    current_positions = positions.copy()
    current_velocities = velocities.copy()

    
    while t1 < tf-deltaT:
        new_positions, new_velocities = evolveParticles(current_positions, current_velocities, masses, dt)
        pos_list[:,:,i] = new_positions 
        vel_list[:,:,i] = new_velocities 
        t1 += deltaT
        i+=1
        current_positions = new_positions 
        current_velocities = new_velocities
        """
        This loop evolves the particles.
        it plugs in the initial positions, then puts the new positions into pos_list as new_positions.
        Then another step is taken, and it plugs in the last output of new_positions as current_positions. 
        It does this until t1 is the same value as tf
        """
    return pos_list, vel_list    


positionsse, velocitiesse =evolve(positions, velocities, masses, deltaT, tf, epsilon = 0.0)
#evolves Sun and Earth from the initial positions, velocities, and masses to show their positions.

#creates scatter plot of the positions of the Sun and Earth for each time step. Since the time step is small, it appears as a line
pl.scatter(positionsse[0,0,:], positionsse[0,1,:], color = 'y')
pl.scatter(positionsse[1,0,:], positionsse[1,1,:], color = 'g')



G = 1
bodyList = []

with open('initial_conditions_Amanda.txt','r') as FILE:
    lines = FILE.readlines()
    for line in lines:
        if line[0] == '#':
            continue
    
         
        bodyList.append(line.split())
    """
    This opens the file, and skips the first line that starts with #.
    line.split seperates each value when the text document has a space. 
    The for loop converts the strings within boyList to floats so the values can be used in evolve particles.
    """

for i in range(len(bodyList)):
    for j in range(7):
        bodyList[i][j] = float(bodyList[i][j])

massList = [bodyList[0][0],bodyList[1][0], bodyList[2][0], bodyList[3][0], bodyList[4][0], bodyList[5][0], bodyList[6][0], bodyList[7][0], bodyList[8][0], bodyList[9][0]]    #Mass Sun, Mass Earth
posList = [bodyList[0][1:4],bodyList[1][1:4],bodyList[2][1:4],bodyList[3][1:4],bodyList[4][1:4],bodyList[5][1:4],bodyList[6][1:4],bodyList[7][1:4],bodyList[8][1:4],bodyList[9][1:4]] #X,Y,Z position of Sun, Earth
velList = [bodyList[0][4:],bodyList[1][4:],bodyList[2][4:],bodyList[3][4:],bodyList[4][4:],bodyList[5][4:],bodyList[6][4:],bodyList[7][4:],bodyList[8][4:],bodyList[9][4:]]   #X,Y,Z velocity of Sun, Earth

masses1 = np.array(massList)
positions1 = np.array(posList)
velocities1 = np.array(velList)

#X positions
posxList = [bodyList[0][1],bodyList[1][1],bodyList[2][1],bodyList[3][1],bodyList[4][1],bodyList[5][1],bodyList[6][1],bodyList[7][1],bodyList[8][1],bodyList[9][1]] #X position of Sun, Earth
#Y positions
posyList = [bodyList[0][2],bodyList[1][2],bodyList[2][2],bodyList[3][2],bodyList[4][2],bodyList[5][2],bodyList[6][2],bodyList[7][2],bodyList[8][2],bodyList[9][2]] #X position of Sun, Earth
#Z positions
poszList = [bodyList[0][3],bodyList[1][3],bodyList[2][3],bodyList[3][3],bodyList[4][3],bodyList[5][3],bodyList[6][3],bodyList[7][3],bodyList[8][3],bodyList[9][3]] #X position of Sun, Earth

#Plots the initial positions in a 3D plot
fig = pylab.figure()
ax = Axes3D(fig)

ax.scatter(posxList, posyList, poszList)
pl.xlabel('X Position')
pl.ylabel('Y Position')
pl.show()


#X initial velocity
velxList = [bodyList[0][4],bodyList[1][4],bodyList[2][4],bodyList[3][4],bodyList[4][4],bodyList[5][4],bodyList[6][4],bodyList[7][4],bodyList[8][4],bodyList[9][4]] #X position of Sun, Earth
#Y initial velocity
velyList = [bodyList[0][5],bodyList[1][5],bodyList[2][5],bodyList[3][5],bodyList[4][5],bodyList[5][5],bodyList[6][5],bodyList[7][5],bodyList[8][5],bodyList[9][5]] #X position of Sun, Earth
#Z initial velocity
velzList = [bodyList[0][6],bodyList[1][6],bodyList[2][6],bodyList[3][6],bodyList[4][6],bodyList[5][6],bodyList[6][6],bodyList[7][6],bodyList[8][6],bodyList[9][6]] #X position of Sun, Earth

#Plots initial velocities in a 3D plot
fig = pylab.figure()
ax = Axes3D(fig)

ax.scatter(velxList, velyList, velzList)
pl.xlabel('X Velocity')
pl.ylabel('Y Velocity')
pl.show()


deltaT = 1e-3 #Change in time
tf = 10 #final time
dt = deltaT
t1 = 0 #initial time



def gForce(mass1, mass2, radius, epsilon=0.0):
    G = 1
    """
    Calculate the gravitational force between two bodies.

    Input Parameters
    ----------
    mass1 : int, float
        Mass of the first body, in kilograms.
    mass2 : int, float
        Mass of the second body, in kilograms.
    radius : int, float
        Separation of the bodies, in meters.
    0.05 is epsilon, so that when r is very small, the force does not go to infinity.

    Returns
    -------
    Force in Newtons : float
    """
    force = G * mass1 * mass2 / (radius**2+epsilon**2) # Calculate the force
    
    return force


def evolve1(positions1, velocities1, masses1, deltaT, tf, epsilon=0.0):
    """
    This evolves the initial positions and velocities of the ten bodies.
    for each time step deltaT, and adds the new position into pos_list and new velocity into vel_list, which are arrays, then uses that
    Each particle interacts with eachother, and evolve1 takes all forces into consideration when evolving.
    Starts with i = 1, so it iis taking the initial positions, then adds 1 to i for each time step.
    """
    i = 1
   
    t1 = 0
    pos_list = np.zeros([len(masses1), 3, (tf/deltaT)+1])
    vel_list = np.zeros([len(masses1), 3, (tf/deltaT)+1])
    pos_list[:,:,0] = positions1.copy()
    vel_list[:,:,0] = velocities1.copy()
    current_positions = positions1.copy()
    current_velocities = velocities1.copy()

    
    while t1 < tf-deltaT:
        new_positions, new_velocities = evolveParticles(current_positions, current_velocities, masses1, deltaT,epsilon)
        pos_list[:,:,i] = new_positions.copy()
        vel_list[:,:,i] = new_velocities.copy()
        t1 += deltaT
        i+=1
        current_positions = new_positions
        current_velocities = new_velocities
    return pos_list, vel_list   

positionsevol, velocitiesevol = evolve1(positions1, velocities1, masses1, deltaT, 10, epsilon = 0.05)


fig = pylab.figure()
ax = Axes3D(fig)
ax.scatter(positionsevol[0,0,0], positionsevol[0,1,0], positionsevol[0,2,0])

ax.plot(positionsevol[0,0,:], positionsevol[0,1,:], positionsevol[0,2,:])
ax.plot(positionsevol[1,0,:], positionsevol[1,1,:], positionsevol[1,2,:])
ax.plot(positionsevol[2,0,:], positionsevol[2,1,:], positionsevol[2,2,:])
ax.plot(positionsevol[3,0,:], positionsevol[3,1,:], positionsevol[3,2,:])
ax.plot(positionsevol[4,0,:], positionsevol[4,1,:], positionsevol[4,2,:])
ax.plot(positionsevol[5,0,:], positionsevol[5,1,:], positionsevol[5,2,:])
ax.plot(positionsevol[6,0,:], positionsevol[6,1,:], positionsevol[6,2,:])
ax.plot(positionsevol[7,0,:], positionsevol[7,1,:], positionsevol[7,2,:])
ax.plot(positionsevol[8,0,:], positionsevol[8,1,:], positionsevol[8,2,:])
ax.plot(positionsevol[9,0,:], positionsevol[9,1,:], positionsevol[9,2,:])
pl.xlabel('X Position')
pl.ylabel('Y Position')
pl.show()


positionArray = positionsevol.copy()
velocityArray = velocitiesevol.copy()

def writeMe(positionArray, velocityArray):
    ''' Write the position array and velocity array to 2 text files.
        
        Input:  Position Array <3D numpy array>
                Velocity Array <3D numpy array>
        Output: None
        
        Example array
            [ [[5.0, 3.0, 0.0], [4.0, 6.0, 1.0], [3.0, 3.0, 3.0]],
              [[4.3, 2.3, 1.2], [9.8, 8.6, 4.6], [8.9, 7.8, 4.2]] ]

        Array Explained
               | x0   x1  x2 |  | y0   y1  y2 |  | z0   z1  z2 |
               -------------------------------------------------
            [ [[5.0, 3.0, 0.0], [4.0, 6.0, 1.0], [3.0, 3.0, 3.0]],     <-- Particle 0
              [[4.3, 2.3, 1.2], [9.8, 8.6, 4.6], [8.9, 7.8, 4.2]] ]    <-- Particle 1

            This example has 2 particles/bodies; each with x, y, z axis.
            Each axis has 3 points. The first coordinate is (x0, y0, z0),
            then (x1, y1, z1), and finally (x2, y2, z2).
    '''
    # Force numpy to print out the whole array
    np.set_printoptions(threshold = np.nan)
    
    positionArray = np.array(positionArray)
    velocityArray = np.array(velocityArray)
    
    #  Write the position array to a text file
    with open("positionsAmanda.txt", 'w') as WRITE:
        WRITE.write(repr(positionArray))
    
    #  Write the velocity array to a text file
    with open("velocitiesAmanda.txt", 'w') as WRITE:
        WRITE.write(repr(velocityArray))
    
    # Reset the threshold to 1000
    np.set_printoptions(threshold = 1000)

    return

#writeMe(positionArray, velocityArray)



