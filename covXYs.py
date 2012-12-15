"""
covXYs module for JRO MI rotation calculations
establish a namespace with JRO parameters and subroutines
for magnetic geometry (ionosphere independent) and covariance
(ionosphere and geometry dependent) calculations

usage:
import covXYs as covXY
geom=covXY.geom
covs=covXY.covs
px,py,uB,Bo,theta,uphi,uthe=geom(tx,ty,dr)
covs(px,py,uB,Bo,theta,uphi,uthe,Ne,h)
"""
from pylab import *

import jropack.igrf.igrf11 as igrf11 #imports the entire *namespace* of igrf11
igrf_B=igrf11.igrf_B                 # IGRF subroutine that returning Bv
a_igrf=igrf11.a                      # [km] mean Earth radius used in IGRF model
a_igrf

import scipy.constants as pc     #physical constants
q=pc.e                           # [C] electric charge
m=pc.m_e                         # [kg] mass of an electron
eps0=pc.epsilon_0                # [F/m] free space electrical permittivity
mu0=pc.mu_0                      # [H/m] free space magnetic permeability
c=pc.c                           # [m/s] light speed

print q,m,eps0,mu0,c

#wgs84 model of spheroid Earth
a=6378.137                                # [km] - equatorial radius WGS 84
f=1/298.257223563                         # flattening f = (a-b)/a
b=a*(1-f)                                 # [km] - polar radius
eccentricity=sqrt(a**2-b**2)/a            # eccentricity

#spacetime of JRO
year=2004.99     # experiment time needed in IGRF model
latR = deg2rad(-11.947917); lonR = deg2rad(-76.872306)  #geodedic --- like GPS indicates
hR =0.463                                  # [km] - local height above reference ellipsoid
nR=a/sqrt(1-(f*(2-f)*sin(latR)**2))        # radius of curvature
xR=(nR+hR)*cos(latR)*cos(lonR)
yR=(nR+hR)*cos(latR)*sin(lonR)
zR=(nR*(1-eccentricity**2)+hR)*sin(latR)

#north, east and zenith unit vectors at JRO
pR=sqrt(xR**2+yR**2)
uRz=array([cos(latR)*cos(lonR),cos(latR)*sin(lonR),sin(latR)])    # Zenith unit vector
uRe=array([-yR,xR,0])/pR                                          # East  unit vector
uRn=cross(uRz,uRe)                                                # North unit vector

#orthonormal basis vectors including the JRO on-axis direction
dec=deg2rad(-12.88); ha=deg2rad(-(4+37/60))                         
uo=array([cos(dec)*cos(ha/4+lonR),cos(dec)*sin(ha/4+lonR),sin(dec)])    # on axis
ux=cross(uRz,uo); ux=ux/norm(ux)                                        # along building to the right
uy=cross(uo,ux) 
ux,uy,uo

def geom(tx,ty,dr):
    """
sets geom
    """
    tz=sqrt(1-(tx**2 + ty**2))
    uk=tx*ux+ty*uy+tz*uo  # propagation direction to the target 
    #polarization vectors at ground
    px=cross(uk,cross(uk,ux))
    px=px/norm(px)
    py=cross(uk,cross(uk,uy))
    py=py/norm(py)

    xyzT=[xR,yR,zR] #initialize target vector
    uBs=[];Bos=[];thetas=[];uphis=[];uthes=[]; #initialize magneto-ionic parameter lists
    for i in range(200): #loop to targets at different ranges 0, dr, 2*dr, ... 199*dr
        pT=sqrt(xyzT[0]**2+xyzT[1]**2)
        latT=arctan2(xyzT[2],pT)
        lonT=arctan2(xyzT[1],xyzT[0]) 
        rT=norm(xyzT)
        
        XB,YB,ZB,B=igrf_B(year,rT-a_igrf,rad2deg(lonT),rad2deg(latT))   # X = North, Y = East, Z = Down (nT)

        uzT=xyzT/rT                          #zenith unit vector at target location
        ueT=array([-xyzT[1],xyzT[0],0])/pT   #east unit vector = -x*sin(lonR) + y*cos(lonR)
        unT=-cross(ueT,uzT)                  #north unit vector
        
        Bv=XB*unT+YB*ueT-ZB*uzT     # B-field vector
        Bo=norm(Bv)                 # [nT] where it is necessary multiply it by 10^-09 (1 nT = 10^-05 Gauss and 1 Gauss = 10^-04 kg/(C.s)
        uB=Bv/Bo                    # unit vector for magnetic field

        theta=arccos(dot(uB,uk))   # [rad] the angle theta between uB and uk            
        uphi=cross(uB,uk)          # phi unit vector
        uphi=uphi/norm(uphi)
        uthe=cross(uphi,uk)        # theta unit vector
        uthe=uthe/norm(uthe)

        uBs=uBs+[uB]; Bos=Bos+[Bo]; thetas=thetas+[theta]; uphis=uphis+[uphi]; uthes=uthes+[uthe] #add to lists
        xyzT=xyzT+dr*uk #increment the target position vector
    uB=array(uBs); Bo=array(Bos); theta=array(thetas); uphi=array(uphis); uthe=array(uthes) #convert lists to arrays
    return px,py,uB,Bo,theta,uphi,uthe

def covs(px,py,uB,Bo,theta,uphi,uthe,Ne,h):
    f0=49.98e06      # [Hz] operating frequency of radar
    w0=2*pi*f0       # [rad/s] angular operating frequency of radar
    k0=w0/c          # free space propagation constant
    wp=sqrt((Ne*q**2)/(eps0*m))              # plasma frequency
    X=wp**2/w0**2 
    ge=(Bo*1.0e-09)*q/m                      # electron gyrofrequency. Bo in [kg/(C*s)]
    Y=ge/w0;
    YL=Y*cos(theta); YT=Y*sin(theta)
    F1=(0.5*YT**2)/(1-X);
    F2=sqrt(F1**2+YL**2)
    Fo=F1-F2; Fx=F1+F2
    no=sqrt(1-X/(1-Fo)); nx=sqrt(1-X/(1-Fx))
    deln=0.5*(no-nx)                         # difference between the refractive indices
    avn=0.5*(no+nx)                          # avg refractive indices
    ar=Fo/YL                                 # axial ratio showing the polarization changes

    dr=h[1]-h[0]
    denom=1/(1+ar**2)
    phase=k0*deln*dr*1000
    A=(exp(-1j*phase)+ar**2*exp(1j*phase))*denom
    B=2*ar*sin(phase)*denom
    C=exp(1j*k0*dr*1000*avn)
    covXX=[];covXY=[];T=identity(3);R=T
    for i in range(200):    
        S=array([uthe[i],uphi[i]])                #transforms field vector from JROxyz- to th/phi-system
        M=array([[A[i],B[i]],[-B[i],conj(A[i])]]) #propagator of vector [Eth;Eph] across slab i
        M=dot(S.T,dot(M,S))                       #propagator across slab i in xyz-system        
        T=dot(M,T)                                #forward propagator from slab 1 to slab i+1
        R=dot(R,M)                                #back propagator from slab i+1 to slab 1
        RT=dot(R,T)                               #backscatter matrix for slab i+1 
        xp=dot(RT,px.T)                           #backscattered field vector per unit ionosphere    
        vx=dot(px,xp)                             #x dipole voltage response to x-pol tx
        vxx=abs(vx)**2
        vy=dot(py,xp)                             ##x dipole voltage response to x-pol tx 
        vxy=abs(vy)**2    
        covXX=covXX+[X[i]*vxx]
        covXY=covXY+[X[i]*vxy]
    covXX=array(covXX)
    covXY=array(covXY)
    return covXX,covXY

if __name__=='__main__':   
    dr=5. # [km] stepping size
    tx=-0.0245;ty=-tx # direction cosines
    px,py,uB,Bo,theta,uphi,uthe=geom(2*tx,2*ty,dr)
    #ionosphere specification
    h=arange(200)*dr                                # [km] height grid vector
    Ne=8.0e12*exp(-h/180.-50.*exp(-h/60.))     
    covXX,covXY=covs(px,py,uB,Bo,theta,uphi,uthe,Ne,h)
    plot(h,covXX); xlabel('Height (km)'); plot(h,covXY,color='r');
    grid('on'); show()
