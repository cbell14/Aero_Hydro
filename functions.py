import numpy as np
import math

class Source:
    """contains information about sources or sinks"""
    def __init__(self, strength, x, y):
        """Initialize the source or sink

        Params:
        -------
        strength  float, strength of the source or sink
        x, y      float, cartesian locations of the source or sink
        """

        self.strength = strength
        self.x, self.y = x, y

    def velocity(self, X, Y):
        """Computes the cartesian velocity field produced by
        a source or sink

        Params:
        -------
        X, Y     2D array of float, meshgrid
        """

        self.u = (self.strength/(2.0*math.pi)) * (X-self.x) /\
            ((X-self.x)**2 + (Y-self.y)**2)
        self.v = (self.strength/(2.0*math.pi)) * (Y-self.y) /\
            ((X-self.x)**2 + (Y-self.y)**2)

    def streamfunction(self, X, Y):
        """Computes the cartesian streamfunction produced by a
        source or sink

        Params:
        ------
        X, Y     2D array of float, meshgrid
        """

        self.psi = self.strength/(2.0*math.pi) * np.arctan2\
            ((Y-self.y),(X-self.x))
            
class Vortex:
    """Contains information about a vortex"""
    def __init__(self, strength, x, y):
        """Initialize the vortex

        Params:
        -------
        stength  float, strength of the vortex
        x, y     float, cartesian locations of the source or sink
        """

        self.strength = strength
        self.x, self.y = x, y

    def velocity(self, X, Y):
        """Computes the cartesian velocity field created by a vortex

        Params:
        -------
        X, Y     2D array of float, meshgrid
        """

        self.u = +self.strength/(2.0*math.pi) * (Y-self.y) / ((X-self.x)**2 + (Y-self.y)**2)
        self.v = -self.strength/(2.0*math.pi) * (X-self.x) / ((X-self.x)**2 + (Y-self.y)**2)

    def streamfunction(self, X, Y):
        """Computes the cartesian streamfunction created by a vortex

        Params:
        -------
        X, Y     2D array of float, meshgrid
        """

        self.psi = -self.strength/(4.0*math.pi) * np.log((X-self.x)**2 + (Y-self.y)**2)
        
class Doublet:
    """Contains information about a doublet"""
    def __init__(self, strength, x, y):
        """Initialize the doublet

        Params:
        -------
        strength  float, strength of the doublet
        x, y      float, cartesian locations of the doublet
        """

        self.strength = strength
        self.x, self.y = x, y

    def velocity(self, X, Y):
        """Computes the cartesian velocity field created by a
        doublet

        Params:
        ------
        X, Y     2D array of float, meshgrid
        """

        self.u = -self.strength/(2.0*math.pi) * ((X-self.x)**2\
            - (Y-self.y)**2)/((X-self.x)**2 + (Y-self.y)**2)**2
        self.v = -self.strength/(2.0*math.pi) * (2.0 * (X-self.x)\
            *(Y-self.y))/((X-self.x)**2 + (Y-self.y)**2)**2

    def streamfunction(self, X, Y):
        """Computes the cartesian streamfunction created by a
        doublet

        Params:
        ------
        X, Y    2D array of float, meshgrid
        """

        self.psi = -self.strength/(2.0*math.pi) * (Y-self.y)\
            /((X-self.x)**2 + (Y-self.y)**2)
            

class UniformFlow:
    """Contains information about uniform flow"""
    def __init__(self, u_inf, alpha, N, X, Y):
        """initialize the uniform flow

        Params:
        -------
        u_inf    float, freestream speed
        alpha    float, angle of attack
        N        float, number of points
        X, Y     2D array of float, meshgrid
        """

        self.u_inf = u_inf
        self.alpha = alpha
        self.N = N
        self.X, self.Y = X, Y

    def velocity(self):
        """Computes the uniform cartesian velocity field"""

        self.u = self.u_inf * (np.cos((self.alpha + np.zeros((self.N,self.N), dtype=float))))
        self.v = self.u_inf * (np.sin((self.alpha + np.zeros((self.N,self.N), dtype=float))))

    def streamfunction(self):
        """Computes the uniform cartesian streamfunction"""

        self.psi = self.u_inf * ((self.Y*np.cos(self.alpha)) - (self.X*np.sin(self.alpha)))




def gen_uform_grid(N, x_s, x_e, y_s, y_e):
    """Generates uniform cartesian grid of points for streamline calculations
    
    Params:
    -------
    N          float, number of points in each direction
    x_s, y_s   float, starting point of x and y
    x_e, y_e   float, ending point of x and y
    
    Returns:
    -------
    X,Y        2D array of float, meshgrid
    """
    x = np.linspace(x_s,x_e,N) #1-D array for x
    y = np.linspace(y_s,y_e,N) #1-D array for y
    X, Y = np.meshgrid(x,y)    #generates mesh grid
    
    return X, Y

def vel_uniform_flow_alpha(u_inf, alpha, N, X, Y):
    """Generates uniform cartesian flow velocity field
    
    Params:
    ------
    u_inf    float, free stream speed
    alpha    float, angle of attack
    N        float, number of points
    X,Y      2D array of float, meshgrid

    Returns:
    -------
    u,v      2D array of float, x and y velocities
    """

    # computes the freestream velocity field
    u = u_inf * (np.cos((alpha + np.zeros((N, N), dtype=float))))  #d(psi)/dy
    v = u_inf * (np.sin((alpha + np.zeros((N, N), dtype=float)))) #-d(psi)/dx
    return u, v

def sf_uniform_flow_alpha(u_inf, alpha, X, Y):
    """Generates uniform cartesian flow stream-function
    
    Params:
    ------
    u_inf    float, free stream speed
    alpha    float, angle of attack
    X,Y      2D array of float, meshgrid

    Returns:
    -------
    psi      2D array of float, streamfunction
    """

    # computes the stream-function
    psi = u_inf * ((Y*np.cos(alpha)) - (X*np.sin(alpha)))

    return psi
    
def vel_uniform_flow(u_inf, N, X, Y):
    """Generates uniform cartesian flow velocity field
    
    Params:
    ------
    u_inf    float, free stream speed
    N        float, number of points
    X,Y      2D array of float, meshgrid

    Returns:
    -------
    u,v      2D array of float, x and y velocities
    """

    # computes the freestream velocity field
    u = u_inf * np.ones((N, N), dtype=float)
    v = np.zeros((N, N), dtype=float)
    return u, v

def sf_uniform_flow(u_inf, X, Y):
    """Generates uniform cartesian flow stream-function
    
    Params:
    ------
    u_inf    float, free stream speed
    X,Y      2D array of float, meshgrid

    Returns:
    -------
    psi      2D array of float, streamfunction
    """

    # computes the stream-function
    psi = u_inf * Y

    return psi

def vel_source_sink(strength, xs, ys, X, Y):
    """Returns the cartesian velocity field generated by a source/sink.
    
    Params:
    ------
    strength  float, strength of the source/sink
    xs, ys    float, coordinates of the source/sink
    X, Y      2D array of float, mesh grid

    Returns:
    -------
    u, v      2D array of float, x and y velocities
    """
    u = strength/(2*np.pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*np.pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    
    return u, v

def sf_source_sink(strength, xs, ys, X, Y):
    """Returns the cartesian stream-function generated by a source/sink.
    
    Params
    ------
    strength   float, strength of the source/sink
    xs, ys     float, coordinates of the source/sink
    X, Y       2D array of float, mesh grid

    Returns
    -------
    psi        2D array of float, streamfunction
    """
    psi = strength/(2*math.pi)*np.arctan2((Y-ys), (X-xs))
    
    return psi

def vel_doublet(strength, xd, yd, X, Y):
    """Returns the cartesian velocity field generated by a doublet.
    
    Params:
    ---------
    strength    float, strength of the doublet
    xd, yd      float, coordinates of the doublet
    X, Y        2D array of float, mesh grid

    Returns:
    -------
    u, v        2D array of float, x and y velocities
    """
    u = - strength/(2*math.pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = - strength/(2*math.pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    
    return u, v

def sf_doublet(strength, xd, yd, X, Y):
    """Returns the cartesian stream-function generated by a doublet.
    
    Params
    ---------
    strength    float, strength of the doublet
    xd, yd      float, coordinates of the doublet
    X, Y        2D array of float, mesh grid

    Returns
    -------
    psi         2D array of float, streamfunction
    """
    psi = - strength/(2*math.pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    
    return psi

def vel_vortex(strength, xv, yv, X, Y):
    """Returns cartesian velocity field from vortex
    
    Params:
    ------
    strength   float, strength of vortex
    xv, yv     float, vortex coordinates
    X, Y       2D array of float, meshgrid
    
    Returns
    -------
    u,v        2D array of float, x and y velocities
    """

    u = + strength/(2*math.pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*math.pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    
    return u, v

def sf_vortex(strength, xv, yv, X, Y):
    """Returns cartesian stream-function from vortex
    
    Params:
    ------
    strength   float, strength of vortex
    x_v, y_v   float, vortex coordinates
    X, Y       2D array of float, meshgrid
    
    Returns
    -------
    psi         2D array of float, streamfunction
    """
    
    psi = strength/(2.0*math.pi)*np.log(np.sqrt((X-xv)**2+(Y-yv)**2))

    return psi

def vel_doublet_cylind(strength, xd, yd, X, Y):
    """Returns the cylindrical velocity field generated by a doublet.
    
    Params:
    ---------
    strength     float, strength of the doublet.
    xd, yd       float, coordinates of the doublet.
    X, Y         2D array of float, mesh grid.
    
    Returns
    --------
    v_r, v_theta 2D array of float, r and theta velocities
    """
    r = np.sqrt((X - xd)**2 + (Y - yd)**2)
    theta = np.arctan(Y/X)
    
    v_r = - strength/(2*math.pi)*((np.cos(theta))/r**2)
    v_theta = - strength/(2*math.pi)*((np.sin(theta))/r**2)
    
    return v_r, v_theta

def sf_doublet_cylind(strength, xd, yd, X, Y):
    """Returns the cylindrical stream-function generated by a doublet
    
    Params:
    ---------
    strength  float, strength of the doublet.
    xd, yd    float, coordinates of the doublet.
    X, Y      2D array of float, mesh grid.
    
    Returns:
    -------
    psi       2D array of float, streamfunction
    """
    r = np.sqrt((X - xd)**2 + (Y - yd)**2)
    theta = np.arctan(Y/X)
    
    psi = - strength/(2*math.pi)*((np.sin(theta))/r)
    
    return psi

def vel_uniform_flow_cylind(u_inf, xd, yd, alpha, N, X, Y):
    """Generates uniform cylindrical flow velocity field
    
    Params:
    ------
    u_inf        float, free stream speed
    xd, yd       float, coordinates of the doublet.
    alpha        float, angle of attack
    N            float, number of points
    X, Y         2D array of float, mesh grid.
    
    Returns
    -------
    v_r, v_theta 2D array of float, r and theta velocities
    """

    r = np.sqrt((X - xd)**2 + (Y - yd)**2)
    theta = np.arctan(Y/X)
    
    # computes the freestream velocity field
    v_r = u_inf * np.cos(theta - alpha)
    v_theta = - u_inf * np.sin(theta - alpha)
    
    return v_r, v_theta

def sf_uniform_flow_cylind(u_inf, xd, yd, alpha, X, Y):
    """Returns the cylindrical stream-function generated by uniform flow
    
    Params:
    ------
    u_inf    float, free stream speed
    xd, yd   float, coordinates of the doublet
    alpha    float, angle of attack
    X, Y     2D array of float, mesh grid

    Returns
    -------
    psi      2D array of float, streamfunction
    """
    r = np.sqrt((X - xd)**2 + (Y - yd)**2)
    theta = np.arctan(Y/X)
    
    psi = u_inf * (((r*np.sin(theta))*np.cos(alpha))\
                   - ((r*np.cos(theta))*np.sin(alpha)))
    
    return psi
    
def vel_combined_source_sink(strength, xs, ys, X, Y, N):
    """Calculates the velocities for several combined
    sources or sinks
    
    Params:
    ------
    strength   float, strength of sources/sinks
    xs, ys     float, coordinates of the sources/sinks
    X, Y       2D array of float, mesh grid
    N          float, number of points
    
    Returns:
    -------
    u_all, v_all      2D arrays of float, x and y velocities
    """
    
    #u_s, v_s = np.zeros((N,N)), np.zeros((N,N))
    u_all, v_all = np.zeros((N,N)), np.zeros((N,N))
    
    for i in range(0,len(source_strength)):
        u_s, v_s = vel_source_sink(strength[i], xs[i], ys[i], X, Y)
        u_all += u_s
        v_all += v_s
        
    return u_all, v_all
    
def sf_combined_source_sink(strength, xs, ys, X, Y, N):
    """Calculates the stream functions for several combined
    sources or sinks
    
    Params:
    ------
    strength   float, strength of sources/sinks
    xs, ys     float, coordinates of the sources/sinks
    X, Y       2D array of float, mesh grid
    N          float, number of points
    
    Returns:
    -------
    psi_s        2D array of float, streamfunction
    """
    
    #psi_s = np.zeros((N,N))
    psi_all = np.zeros((N,N))
    
    for i in range(0,len(source_strength)):
        psi_s = sf_source_sink(strength[i], xs[i], ys[i], X, Y)
        psi_all += psi_s
        
    return psi_all
