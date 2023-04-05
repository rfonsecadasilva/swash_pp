# -*- coding: utf-8 -*-
"""
Linear Wave Theory functions
function disper imported from openearthtools (https://github.com/openearth)
"""
import numpy as np

def disper(w, h, g=9.81):
    """Source (https://github.com/openearth)
    Calculate wave number k using  linear dispersion relation.
    (absolute error in k*h < 5.0e-16 for all k*h)
    Args:
        w (float): angular frequency (radians)
        h (float): depth (m)
        g (float, optional): gravity acceleration (m2/s2). Default to 9.81.
    """
    if not type(w) == list:
        w=[w]
    if not type(h) == list:
        h=[h]
    w2 = [iw**2 * ih / g for (iw,ih) in zip(w,h)]
    q = [iw2 / (1 - np.exp(-(iw2**(5/4))))**(2/5) for iw2 in w2]

    thq = np.tanh(q)
    thq2 = 1 - thq**2

    a = (1 - q * thq) * thq2
    b = thq + q * thq2
    c = q * thq - w2

    D = b**2 - 4 * a * c
    arg = (-b + np.sqrt(D)) / (2 * a)
    iq = np.where(D < 0)[0]
    if iq:
        print(iq)
        arg[iq] = -c[iq] / b[iq]
    q = q + arg

    k = np.sign(w) * q / h
    if np.isnan(k).any():
        k = np.array(k)
        k[np.isnan(k)] = 0

    return k
    
def shoal(T,h):
    """
    Return shoaling coefficient Ks based on linear theory
    Args:
        T (float): wave period (s)
        h (float): depth (m)
    """
    if (type(h) == list) or (type(h) == np.ndarray):
        h = np.array(h)
        T = np.array([T for i in range(len(h))])
        H = np.array([(2*(np.cosh(disper(2*np.pi/T[key],h[key])[0]*h[key]))**2/(np.sinh(2*disper(2*np.pi/T[key],h[key])[0]*h[key])+2*disper(2*np.pi/T[key],h[key])[0]*h[key]))**(0.5) for key,value in enumerate(h)])
    else:
        H = (2*(np.cosh(disper(2*np.pi/T,h)*h))**2/(np.sinh(2*disper(2*np.pi/T,h)*h)+2*disper(2*np.pi/T,h)*h))**(0.5)
    return H


def etau_separation(t,eta,u,h,hab,fcutoff):
    """Separate incoming and outgoing water levels and velocities.
    Return etai, etar, ui, and ur
    (eta for water level and u for velocity;
    i for incoming and r for reflected)
    Args:
        t (np array): time in s
        eta (np array): water level time series in m
        u (np array): incoming velocity in m/s (in wave propagation direction) 
        h (float): mean water depth in m (including setup)
        hab (float): height above bed where velocity is taken (positive) in m
        fcutoff (float): maximum frequency to avoid blowing up noise
    
    Ref: Dynamics of Wave Setup over a Steeply Sloping Fringing Reef, 2005,
    J. Phys. Oceanogr., 45, 3005-3023, doi:10.1175/Jpo-D-15-0067.1.
    """
    from numpy.fft import fft,ifft
    import numpy as np
    nantozero =lambda x: np.where(~np.isnan(x),x,0)
    g=9.81
    N = len(eta) # number of samples
    dt=np.mean(np.diff(t))
    f=np.linspace(0,1,int(N/2)+1)/(2*dt) # frequency in Hz
    omega=2*np.pi*f #angular freq
    #omega_2D = np.vstack([omega]*hab.shape[0]) #2D mode (time series)
    k=disper(omega,h,g).squeeze()
    #k_2D = np.vstack([k]*hab.shape[0]) # 2D mode (time series)
    #hab_2D=np.transpose(np.vstack([hab]*omega.shape[0]))
    Ku=np.cosh(k*hab)/np.cosh(k*h) # velocity response factor
    Ku[f>=fcutoff]=0
    rhs=g*k/omega*Ku # rhs of B1A, Buckley et al., 2015
    #water levels
    Feta,Fu=fft(eta)[:N//2+1],fft(u)[:N//2+1] # one-sided fourier for total eta and u
    Feta[0],Fu[0]=0,0 # detrend
    Fetai,Fetar = (rhs*Feta+Fu)/(2*rhs),(rhs*Feta-Fu)/(2*rhs) # one-sided fourier for etai and etar
    Fetai[rhs==0],Fetar[rhs==0]=0,0
    Fetai[f>=fcutoff],Fetar[f>=fcutoff]=0,0
    Fetai,Fetar=np.concatenate((Fetai,np.conj(Fetai[1:])[::-1]),axis=0),np.concatenate((Fetar,np.conj(Fetar[1:])[::-1]),axis=0) # add reversed complex conjugate
    Fetai,Fetar=nantozero(Fetai),nantozero(Fetar)
    #velocities
    Fui,Fur = rhs*Fetai[:N//2+1],-rhs*Fetar[:N//2+1] # one-sided fourier for etai and etar
    Fui[rhs==0],Fur[rhs==0]=0,0
    Fui[f>=fcutoff],Fur[f>=fcutoff]=0,0
    Fui,Fur=np.concatenate((Fui,np.conj(Fui[1:])[::-1]),axis=0),np.concatenate((Fur,np.conj(Fur[1:])[::-1]),axis=0) # add reversed complex conjugate
    Fui,Fur=nantozero(Fui),nantozero(Fur)
    # inverse fft
    etai,etar=ifft(Fetai),ifft(Fetar)
    ui,ur=ifft(Fui),ifft(Fur)
    return etai,etar,ui,ur