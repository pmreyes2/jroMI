from pylab import *
import covXYs as covXY
geom=covXY.geom
covs=covXY.covs

dr=5.
h=arange(200)*dr
for ir in range(12):
    tx=-0.0245*(ir-4);ty=-tx 
    px,py,uB,Bo,theta,uphi,uthe=geom(tx,ty,dr) 
    Ne=8.0e12*exp(-h/180.-50.*exp(-h/60.))     
    covXX,covXY=covs(px,py,uB,Bo,theta,uphi,uthe,Ne,h)
    subplot(12,1,ir)
    plot(h,covXX); xlabel('Height (km)'); plot(h,covXY,color='r');

show()
