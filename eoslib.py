
__author__ = "Guenter Schneider"

from numpy import ones, array, c_, newaxis
from scipy import linalg
from scipy.optimize import leastsq


def eos_harmonic(p,v):
    """Harmonic equation of state: e ~ v**2.
    p = [e0,v0,b0], with 
    e0 = energy at v0
    v0 = equilibrium volume
    b0 = bulk modulus
    """
    if len(p) < 3:
        raise Error
    e0 = p[0]
    v0 = p[1]
    b0 = p[2]
    return e0 + 0.5*b0*(v-v0)**2/v0


def eos_murnaghan(p,v):
    """Murnaghan equation of state: bp = db0/dp = const.
    p = [e0,v0,b0,bp], with
    e0 = energy at v0
    v0 = equilibrium volume
    b0 = bulk modulus
    bp = db0/dp with p := pressure
    """
    if len(p) < 4:
        raise Error
    e0 = p[0]
    v0 = p[1]
    b0 = p[2]
    bp = p[3]
    return e0 + b0*v/bp*((v0/v)**bp/(bp-1)+1.0) - b0*v0/(bp-1)


def eos_birch(p,v):
    """Birch equation of state: 3rd order eos
    p = [e0,v0,b0,bp], with
    e0 = energy at v0
    v0 = equilibrium volume
    b0 = bulk modulus
    bp = db0/dp with p := pressure
    """
    if len(p) < 4:
        raise Error
    e0 = p[0]
    v0 = p[1]
    b0 = p[2]
    bp = p[3]
    return (e0 + 9.0/16.0*b0*v0*(bp*((v0/v)**(2.0/3.0)-1.0)**3
               + (6.0-4.0*(v0/v)**(2.0/3.0))*((v0/v)**(2.0/3.0)-1.0)**2))

def pv_birch(p,v):
    """Birch equation of state: 3rd order eos
    p = [e0,v0,b0,bp], with
    e0 = energy at v0
    v0 = equilibrium volume
    b0 = bulk modulus
    bp = db0/dp with p := pressure
    """
    if len(p) < 4:
        raise Error
    e0 = p[0]
    v0 = p[1]
    b0 = p[2]
    bp = p[3]
    vr = v0/v
    return  ( 1.5*b0*(vr**(7./3.)-vr**(5./3.))*
                             (1.0+0.75*(bp-4.0)*(vr**(2./3.)-1.0)) )



class FitError(Exception):
    def __str__(self):
        return repr(self.value)


def parabolicLeastSquaresFit(xdata,ydata):
    """Linear least squares fit of data (xdata,ydata) to c[0]+c[1]*x+c[2]*x**2
    xdata -- x values of data
    ydata -- y values of data. len(ydata) = len(xdata)
    Return values:
    c -- coefficient of parabola
    res -- residual
    """
    one = ones(len(xdata))
    x = array(xdata)
    y = array(ydata)
    A = c_[one[:,newaxis],x[:,newaxis],(x**2)[:,newaxis]]
    c,res,rank,sigma = linalg.lstsq(A,y)
    return c, res


def nonlinearLeastSquaresFit(f,p_guess,xdata,ydata,ftol=1e-8,xtol=1e-8):
    """Nonlinear least squares fit of data (xdata,ydata) with function f.
    f -- function f(p,x) with is a len(p) parameters.
    p_guess -- initial guess for parameters p.
    xdata -- x values of data
    ydata -- y values of data. len(ydata) = len(xdata)
    Optional inputs:
    ftol -- Relative error desired in the sum of squares.
    xtol -- Relative error desired in the approximate solution.
    Return values:
    p -- parameters
    res -- residual
    """
    def residual(p):
        r = []
        for x,y in zip(xdata,ydata):
            r.append(y-f(p,x))
        return array(r)
    
    p,cov_x,infodict,mesg,status = leastsq(residual,p_guess,ftol=ftol,
                                           xtol=xtol,full_output=True)
#    p,cov_x,infodict,mesg,status = leastsq(residual,p_guess,ftol=ftol,xtol=xtol,
#                                           full_output=True,warning=True)
    if status < 1 or status > 4:
        raise FitError(mesg)
    res  = sum(residual(p)**2)
    return p[0],p[1],p[2],p[3],res


def eos_fit(eos,vol,ene,units=('Ang','eV')):
    """Fit of volume, energy data to equation of state.
    eos -- Equation of state = {'harmonic','Murnaghan','Birch'} (only first letter counts)
    vol -- volume data
    ene -- energy data
    units -- Units of volume and energy data, used to convert Bulk modulus to GPa
             units = {'Ang','bohr','eV','Htr','Ryd'} (only first letter counts)
    Return values:
    v0 -- equilibrium volume
    b0 -- bulk modulus in GPa
    bp -- pressure derivative of bulk modulus
    e0 -- energy at v0
    res -- residual
    """
    c,res = parabolicLeastSquaresFit(vol,ene)
    b0 = -c[1]
    v0 = -0.5*c[1]/c[2]
    e0 = c[0]+c[1]*v0+c[2]*v0**2
    if eos[0].lower()=='h':
        f = eos_harmonic
        bp = 0.0
    else:
        bp = 4.0
        if eos[0].lower()=='m':
            f = eos_murnaghan
        elif eos[0].lower()=='b':
            f = eos_birch
            e0,v0,b0,bp,res = nonlinearLeastSquaresFit(eos_murnaghan,(e0,v0,b0,bp),vol,ene)
        else:
            raise ArgumentError
        e0,v0,b0,bp,res = nonlinearLeastSquaresFit(f,(e0,v0,b0,bp),vol,ene)
    u = []
    if len(units)>2:
        raise ArgumentError
    for unit in units:
        u.append(unit[0].lower())
    bmod = b0
    if 'h' in u:
        bmod = bmod*27.2
    if 'r' in u:
        bmod = bmod*13.2
    if 'b' in u:
        bmod = bmod/0.529177**3
    bmod = bmod*1e-9*1.6e-19/1e-10**3  # bmod in GPa
    return v0,bmod,(e0,v0,b0,bp),res


if __name__=='__main__':

    import matplotlib.pyplot as plt
    import numpy as np
    
    def plot_f(f,a,b,N,symbol='o'):
        x = np.linspace(a,b,N)
        y = []
        for x_val in x:
            y.append(f(x_val))
        plt.plot(x,y,symbol)
    
    # diamond
    v= [3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10]

    e= [-4.621424, -5.804482, -6.70809, -7.398549, -7.922563, -8.316091, -8.606262, -8.813852, -8.954821, -9.042461, -9.08705, -9.0977, -9.07705, -9.03505, -8.9743, -8.8984, -8.8103, -8.7122, -8.60665, -8.4949, -8.3785, -8.2585, -8.1362, -8.0118, -7.88685, -7.7615, -7.63665, -7.51265, -7.3898]

    # graphite
    v1 = [5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5, 10.75, 11, 11.25, 11.5, 11.75, 12, 12.25, 12.5, 12.75, 13, 13.25, 13.5, 13.75, 14]

    e1 = [-6.984607, -7.419216, -7.784433, -8.090376, -8.345294, -8.556517, -8.729426, -8.869525, -8.980762, -9.067014, -9.131185, -9.176347, -9.204615, -9.218072, -9.218746, -9.207845, -9.187136, -9.157585, -9.120264, -9.076312, -9.02624, -8.971119, -8.911418, -8.847637, -8.780511, -8.710437, -8.637559, -8.562698, -8.48598, -8.407513, -8.327904, -8.24726, -8.16564, -8.083573, -8.00095]

    print(v)
    print(e)
   
    plt.figure(1)
    plt.plot(v,e,'ro')
    plt.xlabel('volume [Ang^3]')
    plt.ylabel('energy [eV]')
   
    vol,bmod,p,res = eos_fit('h',v,e,units=('eV','Ang'))
    bder = p[3]
    ene = p[0]
    print("                     volume      B modulus      B'    energy    residual")
    print("Harmonic  EOS fit :", vol,bmod,bder,ene,res)
    plot_f(lambda x: eos_harmonic(p,x), 2.7, 11., 121,symbol='-')
    vol,bmod,p,res = eos_fit('m',v,e,units=('eV','Ang'))
    bder = p[3]
    ene = p[0]
    print("Murnaghan EOS fit :", vol,bmod,bder,ene,res)
    plot_f(lambda x: eos_murnaghan(p,x), 2.7, 11.0, 121,symbol='-')
    vol,bmod,p,res = eos_fit('birch',v,e,units=('eV','Ang'))
    bder = p[3]
    ene = p[0]
    plot_f(lambda x: eos_birch(p,x), 2.7, 11.0, 121,symbol='--')
    print("Birch     EOS fit :", vol,bmod,bder,ene,res)

    v1=5.0
    e1=eos_birch(p,v1)
    p1=pv_birch(p,v1)
    p1GPa = p1*1e-9*1.6e-19/1e-10**3  # pressure in GPa


#    print(bmod, p[2])

    bmod_eVAng3 = p[2]
    bprime = p[3]
    vol0_Ang3 = p[1]
    ene_eV = p[0]

    print("B0 in eV/Ang^3 = ", bmod_eVAng3)
    print("B0'            = ", bprime)
    print("V0 in Ang^3    = ", vol0_Ang3)
    print("E(V) in eV     = ", ene_eV)
    print()
    print("Volume in Ang^3 = ",v1)
    print("Energy in eV    = ",e1)
    print("pressure in eV/Ang^3 = ",p1)
    print("pressure in GPa = ",p1GPa)
    print("p*V in eV    = ",p1*v1)
    print("Enthalpy in eV    = ",e1+p1*v1)
  
    plt.show()
