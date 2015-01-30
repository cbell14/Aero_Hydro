def gen_uform_grid(N,x_start,x_end,y_start,y_end):
    """Generates grid of points for streamline calculations
    
    Params
    -------
    N:         float, number of points in each direction
    x_start:   float, starting point of x
    x_end:     float, ending point of x
    y_start:   float, starting point of y
    y_end:     float, ending point of y
    
    """
    x = np.linspace(x_start,x_end,N) #1-D array for x
    y = np.linspace(y_start,y_end,N) #1-D array for y
    X, Y = np.meshgrid(x,y)          #generates mesh grid
    
    return X,Y

def uniform_flow(u_inf, alpha, N):
    """Generates uniform flow stream-function and velocity field
    
    Params:
    ------
    u_inf    float, free stream speed
    alpha    float, angle of attack
    N        float, number of points
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
    strength  float, strength of the source/sink.
    xs, ys    float, coordinates of the source/sink.
    X, Y      array, mesh grid.
    """
    u = strength/(2*np.pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = strength/(2*np.pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    
    return u, v

def sf_source_sink(strength, xs, ys, X, Y):
    """Returns the stream-function generated by a source/sink.
    
    Params
    ------
    strength   float, strength of the source/sink.
    xs, ys     float, coordinates of the source/sink.
    X, Y       array, mesh grid.
    """
    psi = strength/(2*math.pi)*np.arctan2((Y-ys), (X-xs))
    
    return psi

