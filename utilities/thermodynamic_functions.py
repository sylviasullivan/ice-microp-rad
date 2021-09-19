import numpy as np

# Function to calculate the saturation vapor pressure over liquid.
# psatL [=] Pa; t_in [=] K
def satVapP_liq(t_in):
    R = 8.314             # J mol-1 K-1
    MWw = 18.015/1000     # kg mol-1
    rhoa = 1.395
    a1 = 54.842763
    a2 = -6763.22
    a3 = -4.21
    a4 = 0.000367
    a5 = 0.0415
    a6 = 218.8
    a7 = 53.878
    a8 = -1331.22
    a9 = -9.44523
    a10 = 0.014025
    factor = a7 + a8/t_in + a9*np.log(t_in) + a10*t_in
    psatL = a1 + a2/t_in + a3*np.log(t_in) + a4*t_in + np.arctan(a5*(t_in - a6))*factor
    psatL = np.exp(psatL)
    return psatL

# Function to calculate the saturation vapor pressure over ice.
# psatI [=] Pa and t_in [=] K
def satVapP_ice(t_in):
    a1 = 9.550426
    a2 = -5723.265
    a3 = 3.53068
    a4 = -0.00728332
    psatI = a1 + a2/t_in + a3*np.log(t_in) + a4*t_in
    psatI = np.exp(psatI)
    return psatI


# Function to calculate the potential temperature
# t_in [=] K, p_in [=] Pa
def calc_theta(t_in, p_in):
    p0 = 100000    # Reference pressure of 1000 hPa
    R_cp = 0.286   # Ratio of the gas constant to the heat capacity at constant pressure for air
    return t_in * ( p0 / p_in ) ** ( R_cp )


# Function to calculate the relative humidity with respect to ice
# t_in [=] K, p_in [=] Pa, qv [=] kg kg-1, return [=] %
def calc_RHi(t_in, p_in, qv_in):
    eps = 0.622  # ratio of molar mass of moist air to molar mass of dry air
    psat = satVapP_ice(t_in)
    qsat = eps * psat / (p_in - psat)
    return qv_in / qsat * 100


# Utility function to calculate running mean on a vector time series x over N time points.
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)

# Utility function to calculate running mean on a matrix of time series x over N time points along axis a.
# THIS FUNCTION DOES NOT YET WORK PROPERLY
def running_mean2(x, N, a):
    cumsum = np.cumsum(np.insert(x, 0, 0, axis=a)) 
    #print(cumsum.shape)
    #print(cumsum[:-N].shape)
 
    # Make a the first axis for the calculation and then reverse
    if a != 0:
        cumsum = np.swapaxes(cumsum, a, 0)
        smooth = cumsum[N:] - cumsum[:-N] / float(N)
        return np.swapaxes(smooth, a, 0)
    else:
        return (cumsum[N:] - cumsum[:-N]) / float(N)
