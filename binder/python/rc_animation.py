"""A module for generating animations for explaining fundamental rotation curves concepts.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani

###############################################

def CalculatePosition(radius,velocity,time,dt):
    """
    Calculates positions of multiple objects around a circle after some time with interval dt.

    Any set of units may be used, provided they match -- 
    for example, if using a radius of meters and time in seconds, velocity should be in meters per second and the time increments should be in seconds.

    :parameters:
        radius : [array]
            Radius of the objects from the center. Units to match velocity.
        velocity : [array]
            Velocity of the objects. Units to match radius and time. Length must match length of radius.
        time : [int]
            Maximum time. Units to match velocity. If not divisible by dt, the final time provided may exceed this value by a number less than dt.
        dt : [int]
            Time increment. Units to match time.

    :returns:
        Two [ndarray]s of x-position (first array) and y-position (second array) of all objects. These arrays are stored such that each row represents a single time 
        and each column represents a single object.  
        [1darray] of associated time.  

    :example:
        >>> import numpy as np #To generate the arrays
        >>> radius = np.array([1,2,3,4,5])              # in m
        >>> velocity = np.array([0.1,0.2,0.3,0.4,0.5])  # in m/s
        >>> time = 100                                  # in s
        >>> dt = 33                                     # in s
        >>> CalculatePosition(radius,velocity,time,dt)
        (array([[ 1.        ,  2.        ,  3.        ,  4.        ,  5.        ],
            [-0.98747977, -1.97495954, -2.96243931, -3.94991908, -4.93739885],
            [ 0.95023259,  1.90046518,  2.85069778,  3.80093037,  4.75116296],
            [-0.88919115, -1.77838231, -2.66757346, -3.55676461, -4.44595576],
            [ 0.80588396,  1.61176792,  2.41765187,  3.22353583,  4.02941979]]), 
        array([[ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ],
            [-0.15774569, -0.31549139, -0.47323708, -0.63098278, -0.78872847],
            [ 0.31154136,  0.62308273,  0.93462409,  1.24616545,  1.55770682],
            [-0.45753589, -0.91507179, -1.37260768, -1.83014358, -2.28767947],
            [ 0.59207351,  1.18414703,  1.77622054,  2.36829406,  2.96036757]]), 
        array([  0,  33,  66,  99, 132]))
    """
    
    # Initial conditions
    theta = 0
    xini = radius * np.cos(theta)
    yini = radius * np.sin(theta)
    t = 0
    
    # Store positions and time
    xposition = [xini]
    yposition = [yini]
    storedtime = [t]
    
    # Calculate positions
    while t < time:
        t += dt
        x = radius * np.cos((velocity/radius)*t)
        y = radius * np.sin((velocity/radius)*t)
        
        xposition.append(x)
        yposition.append(y)
        storedtime.append(t)
    
    # Make an array
    xposition = np.array(xposition)
    yposition = np.array(yposition)
    storedtime = np.array(storedtime)
    
    return xposition,yposition,storedtime

###############################################

def MultiplePositions(radius,velocity,time,dt):
    """
    Calculates transposed positions of multiple objects around a circle after some time with interval dt. Setup for the animation.

    Any set of units may be used, provided they match -- 
    for example, if using a radius of meters and time in seconds, velocity should be in meters per second and the time increments should be in seconds.

    :parameters:
        radius : [array]
            Radius of the objects from the center. Units to match velocity.
        velocity : [array]
            Velocity of the objects. Units to match radius and time. Length must match length of radius.
        time : [int]
            Maximum time. Units to match velocity. If not divisible by dt, the final time provided may exceed this value by a number less than dt.
        dt : [int]
            Time increment. Units to match time.

    :returns:
        Two [ndarray]s of x-position (first array) and y-position (second array) of all objects. These arrays are stored such that each row represents a single object 
        and each column represents a single time.  
        [1darray] of associated time.

    :example:
        >>> radius = np.array([1,2,3,4,5])              # in m
        >>> velocity = np.array([0.1,0.2,0.3,0.4,0.5])  # in m/s
        >>> time = 100                                  # in s
        >>> dt = 50                                     # in s
        >>> MultiplePositions(radius,velocity,time,dt)
        (array([[ 1.        ,  0.28366219, -0.83907153],
            [ 2.        ,  0.56732437, -1.67814306],
            [ 3.        ,  0.85098656, -2.51721459],
            [ 4.        ,  1.13464874, -3.35628612],
            [ 5.        ,  1.41831093, -4.19535765]]), 
        array([[ 0.        , -0.95892427, -0.54402111],
            [ 0.        , -1.91784855, -1.08804222],
            [ 0.        , -2.87677282, -1.63206333],
            [ 0.        , -3.8356971 , -2.17608444],
            [ 0.        , -4.79462137, -2.72010555]]), 
        array([ 0, 50]))
    """
    
    # Stop the calculation when the outermost point takes a whole revolution
    # Outermost point position
    outerposition = radius[len(radius)-1]
    
    # Calculate the positions of outermost point:
    xouter = CalculatePosition(radius[len(radius)-1],velocity[len(radius)-1],time,dt)[0]
    
    # Circumference of the outer circle
    circouter = 2*np.pi*outerposition
    
    # Distance the outer object traveled
    distance = 0
    istop = 0
        
    # New time
    storedtime = CalculatePosition(radius[0],velocity[0],time,dt)[2]
    for t in storedtime:
        if distance < circouter:
            distance = velocity[len(radius)-1]*t
        else:
            istop = np.where(storedtime == t)[0]     # find index of the numpy array 
            break
                   
    istop = int(istop)
    newstoredtime = storedtime[:istop-1]
            
    xmultiple = []
    ymultiple = []

    for i in range(len(radius)):
        x = CalculatePosition(radius[i],velocity[i],time,dt)[0]
        y = CalculatePosition(radius[i],velocity[i],time,dt)[1]
        xmultiple.append(x)
        ymultiple.append(y)

    xmultiple = np.array(xmultiple)
    ymultiple = np.array(ymultiple)
    
    return xmultiple, ymultiple, newstoredtime

###############################################

def PlotRotationCurve(radius,velocity,title,
                      xlabel='Radius (km)',ylabel='Velocity (km/s)',
                      xlim=1,
                      ylim=0.1):
    """
    Plot rotation curve, given the radius and velocity.

    Any set of units may be used, provided they match -- 
    for example, if using a radius of meters, velocity should be in meters per unit time and the labels should correctly indicate both.

    :parameters:
        radius : [array]
            Radius of the objects from the center. Units to match velocity.
        velocity : [array]
            Velocity of the objects. Units to match radius and time. Length must match length of radius.
        title : [string]
            Title of the plot.
        xlabel : [string]
            X-label of the plot. 
            
            :default: 'Radius (km)'
        ylabel : [string]
            Y-label of the plot. 
            
            :default: 'Velocity (km/s)'
        xlim : [float]
            X-limit of the plot. 
            
            :default: 1
        ylim : [float]
            Y-limit of the plot. 
            
            :default: 0.1

    :returns:
        None

    :example:
        >>> radius = np.array([1,2,3,4,5])              # in m
        >>> velocity = np.array([0.1,0.2,0.3,0.4,0.5])  # in m/s
        >>> PlotRotationCurve(radius,velocity,'Rigid Body Rotation Curve')

        .. image:: ../../images/PlotRotationCurveExample.png
            :alt: A plot, titled 'Rigid Body Rotation Curve', with the radii above on the x-axis and the velocities above on the y-axis
    """
    
    # Convert title to string
    title = str(title)
    
    # Plot
    fig = plt.figure(figsize=(6,6),dpi=100)
    ax = plt.axes()
    
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.set_xlabel(xlabel,color='white')
    ax.set_ylabel(ylabel,color='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.xaxis.label.set_color('white')
    ax.yaxis.label.set_color('white')
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    
    plt.title(title,color='white',fontsize='14')
    plt.scatter(radius,velocity,color='khaki')
    plt.plot(radius,velocity,color='white')
    plt.xlim(0,np.max(radius)+xlim)
    plt.ylim(0,np.max(velocity)+ylim)

    plt.show()
    
###############################################

def MakeAnimation(radius,velocity,time,dt,filename,title,
                 xlim=1,ylim=1,
                 size=False,masses=None):
    """
    Animation of rotating objects around a circle.

    Any set of units may be used, provided they match -- 
    for example, if using a radius of meters and time in seconds, velocity should be in meters per second and the time increments should be in seconds. 

    :parameters:
        radius : [array]
            Radius of the objects from the center. Units of km.
        velocity : [array]
            Velocity of the objects. Units to match time with radial units of km. Length must match length of radius.
        time : [int]
            Maximum time. Units to match velocity. If not divisible by dt, the final time provided may exceed this variable by a number less than dt.
        dt : [int]
            Time increment. Units to match time.
        filename : [string]
            File name to save animation.
        title : [string]
            Title of the plot.
        xlim : [float]
            X-limit of the plot. 
            
            :default: 1
        ylim : [float]
            Y-limit of the plot. 
            
            :default: 0.1
        size : [bool]
            Size of dots, based on masses. If True, the sizes of dots depend on masses. 
            
            :default: False
        masses : [array]
            Masses of objects, when needed. 
            
            :default: None

    :returns:
        None

    :example:
        >>> radius = np.array([1,2,3,4,5])              # in m
        >>> velocity = np.array([0.1,0.2,0.3,0.4,0.5])  # in m/s
        >>> time = 100
        >>> dt = 1
        >>> MakeAnimation(radius, velocity, time, dt, 
        ...     filename='images/solarsystem.gif', 
        ...     title='Planet-like Rotation')    

        .. image:: ../../images/solarsystem.gif
            :alt: An animation of planet-like rotation, in which the five planets are, coincidentally, moving like a rigid body. The x-axis and y-axis are both labeled with units of km with ranges of -6 to +6. Below the title 'Planet-like Rotation' is a counter indicating the current number of seconds, in scientific notation, resetting itself periodically. Other than this timer, the animation appears continuious.
    """
    
    # Extract x and y positions, and time
    xpositions = MultiplePositions(radius,velocity,time,dt)[0]
    ypositions = MultiplePositions(radius,velocity,time,dt)[1]
    #storedtimes = CalculatePosition(radius,velocity,time,dt)[2]
    storedtimes = MultiplePositions(radius,velocity,time,dt)[2]
    
    # Sizes of dots based on masses
    if size == True:   # use an array of masses as an input for sizes
        area = [s * 5e-25 for s in masses]
    if size == False:  # use a default size
        area = 100
                      
    # Create a movie write object, set frame rate
    writer = ani.FFMpegWriter(fps=25)

    # Create a figure, 8"x8" in size
    fig = plt.figure(figsize=(8,8))
    
    # Change background color of the plot
    fig.patch.set_facecolor('black')
    ax = plt.axes()
    ax.set_facecolor('black')
    ax.set_facecolor('black')
    ax.set_xlabel('x (km)',color='white')
    ax.set_ylabel('y (km)',color='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
  
    # Convert filename and title to string
    filename = str(filename)
    title = str(title)

    # Set things up to save frames to a movie:
    #   fig = the figure the writer will record from
    with writer.saving(fig, filename, 100):

    # Loop
        i = 0            # start counter
        for t in storedtimes:
            plt.cla()
            for r in radius:
                    circle = plt.Circle((0, 0), r, color='white', fill=False)
                    plt.gca().add_patch(circle)
            plt.scatter(xpositions[:,i],ypositions[:,i],s=area,color='khaki')
            plt.suptitle(title,color='white',fontsize='18')
            plt.title('{:.1e} seconds'.format(t),color='white')
            plt.xlim(-np.max(radius)-xlim,np.max(radius)+xlim)
            plt.ylim(-np.max(radius)-1,np.max(radius)+ylim)
            plt.xlabel('x (km)',color='white',fontsize='14')
            plt.ylabel('y (km)',color='white',fontsize='14')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            i += 1
            
            # Save the current plot as a movie frame
            writer.grab_frame()
        
        plt.close(fig)    # Do not display the image