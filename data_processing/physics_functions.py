import numpy as np

# Functions to compute physical quantities from the simulation/obs data
# Resources for Goff-Gratch formula: http://cires1.colorado.edu/~voemel/vp.html

def get_esat(temp, type='liq'):
    """Compute saturation vapor pressure over water or ice using the formula from
    Goff-Gratch (1946)
    
    Parameters
    ----------
    temp : float
        Temperature in Kelvin
    type : str
        Type of formula to use: liq or ice

    Returns
    -------
    esat : float
        Saturation vapor pressure over water or ice in Pa
    """
    if type == 'liq':
        esat = 100 * 1013.246 * 10**(-7.90298*(373.16/temp-1) + 5.02808*np.log10(373.16/temp) - 
                                   1.3816e-7*(10**(11.344*(1-temp/373.16))-1) + 8.1328e-3 * 
                                   (10**(-3.49149*(373.16/temp-1))-1)) 
    elif type == 'ice':
        esat =  100 * 6.1071 * 10**( -9.09718*(273.16/temp - 1) - 3.56654*np.log10(273.16/ temp) +
               0.876793*(1 - temp/ 273.16))
    else:
        raise ValueError('type must be either liq or ice')
    
    return esat

def get_spec_hum_gg(temp, pres, rh):
    """Compute specific humidity based on rel humidity (0 to 1).
    
    Parameters
    ----------
    temp : float
        Temperature in Kelvin
    pres : float
        Pressure in Pa
    rh : float
        Relative humidity in percent
    
    Returns
    -------
    q : float
        Specific humidity in kg/kg
    """
    R_d = 287.04    # gas constant of dry air, in J kg^-1 K^-1
    R_v = 461.5     # gas constant of water vapour, in J kg^-1 K^-1
    M_dv = R_d / R_v # molar mass ratio , in ()

    esat_w = get_esat(temp, 'liq')
    e = rh * esat_w

    q = M_dv * e / (e*(M_dv - 1) + pres)

    return q

def get_abs_hum_gg(temp, rh):
    """Compute absolute humidity based on rel humidity (0 to 1).
    
    Parameters
    ----------
    temp : float
        Temperature in Kelvin
    rh : float
        Relative humidity in percent
    
    Returns
    -------
    ah : float
        Absolute humidity in kg/m^3
    """
    R_v = 461.5     # gas constant of water vapour, in J kg^-1 K^-1

    e_env = rh * get_esat(temp, 'liq')
    ah = e_env / (R_v * temp)
    return ah