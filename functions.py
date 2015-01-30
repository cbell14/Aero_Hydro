import numpy as np
import math

def gen_uform_grid(N,x_s,x_e,y_s,y_e):
    """Generates uniform grid of points for streamline calculations
    
    Params
    -------
    N:     float, number of points in each direction
    x_s:   1D array of float, starting point of x
    x_e:   1D array of float, ending point of x
    y_s:   1D array of float, starting point of y
    y_e:   1D array of float, ending point of y
    
    Returns
    -------
    X,Y    2D array of float, meshgrid
    """
    x = np.linspace(x_s,x_e,N) #1-D array for x
    y = np.linspace(y_s,y_e,N) #1-D array for y
    X, Y = np.meshgrid(x,y)    #generates mesh grid
    
    return X,Y

def uniform_flow(u_inf, alpha, N, X, Y):
    """Generates uniform flow stream-function and velocity field
    
    Params:
    ------
    u_inf    float, free stream speed
    alpha    float, angle of attack
    N        float, number of points

    Returns
    -------
    u,v      1D array of float, x and y velocities
    psi      2D array of float, streamfunction
    """

    # computes the stream-function
    psi = u_inf * ((Y*np.cos(alpha)) - (X*np.sin(alpha)))

    # computes the freestream velocity field
    u = u_inf * (np.cos((alpha * np.ones((N, N), dtype=float))))  #dpsi/dy
    v = u_inf * (np.sin((alpha * np.ones((N, N), dtype=float)))) #-dpsi/dx
    return u, v, psi

def vel_source_sink(strength, xs, ys, X, Y):
    """Returns the velocity field generated by a source/sink.
    
    Params
    ------
    strength  float, strength of the source/sink
    xs, ys    1D array of float, coordinates of the source/sink
    X, Y      1D array of float, mesh grid

    Returns
    -------
    u,v       1D array of float, x and y velocities
    """
    u = strength/(2*np.pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*np.pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    
    return u, v

def sf_source_sink(strength, xs, ys, X, Y):
    """Returns the stream-function generated by a source/sink.
    
    Params
    ------
    strength   float, strength of the source/sink
    xs, ys     1D array of float, coordinates of the source/sink
    X, Y       1D array of float, mesh grid

    Returns
    -------
    psi        2D array of float, streamfunction
    """
    psi = strength/(2*math.pi)*np.arctan2((Y-ys), (X-xs))
    
    return psi

def vel_doublet(strength, xd, yd, X, Y):
    """Returns the velocity field generated by a doublet.
    
    Params
    ---------
    strength    float, strength of the doublet
    xd, yd      1D array of float, coordinates of the doublet
    X, Y        1D array of float, mesh grid

    Returns
    -------
    u,v         1D array of float, x and y velocities
    """
    u = - strength/(2*math.pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(Y-yd)**2)**2
    v = - strength/(2*math.pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    
    return u, v

def sf_doublet(strength, xd, yd, X, Y):
    """Returns the stream-function generated by a doublet.
    
    Params
    ---------
    strength    float, strength of the doublet
    xd, yd      1D array of float, coordinates of the doublet
    X, Y        1D array of float, mesh grid

    Returns
    -------
    psi         2D array of float, streamfunction
    """
    psi = - strength/(2*math.pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    
    return psi

def vel_vortex(strength, xv, yv, X, Y):
    """Returns velocity field from vortex
    
    Params:
    ------
    strength   float, strength of vortex
    x_v, y_v   1D array of float, vortex coordinates
    X, Y       1D array of float, meshgrid
    
    Returns
    -------
    u,v         1D array of float, x and y velocities
    """

    u = + strength/(2*math.pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*math.pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    
    return u, v

def sf_vortex(strength, xv, yv, X, Y):
    """Returns velocity field from vortex
    
    Params:
    ------
    strength   float, strength of vortex
    x_v, y_v   1D array of float, vortex coordinates
    X, Y       1D array of float, meshgrid
    
    Returns
    -------
    psi         2D array of float, streamfunction
    """
    
    psi = strength/(4*math.pi)*np.log((X-xv)**2+(Y-yv)**2)