"""
Grid2D establishes a namespace to display
the covariances for each tx-ty in a 2-D grid

Usage :
import Grid2D as G2D
Gr2D = G2D.grid_2D
FigGrid = G2D.FigGrid
CovFigGrid = G2D.CovFigGrid

grid,tx,ty = Gr2D(dci,ns)
FigGrid(grid,dci,ns)   # Optional
CovFigGrid(tx,ty,dr,Ne,h)

"""
from pylab import *

import covXYs as covXY

geom = covXY.geom
covs = covXY.covs

def grid_2D(dci,ns):

    """
        This definition sets up the grid composed of direction cosines (tx-ty).
    """

    grid = mgrid[-dci*ns:dci*ns:dci, dci*ns:-dci*ns:-dci] 
    tx = grid[0].T[0]; print 'tx : ',tx    #  direction cosines specified in a way that loop will begin from the upper left corner
    ty = grid[1][0];     print 'ty : ', ty

    return grid, tx, ty

def FigGrid(grid,dci,ns):

    """
        This definition displays the figure of the grid set up by definition 'grid2D' to give the user an idea how the grid looks like.
    """
    
    figure(1)

    mag = dci*ns + dci; mig = -mag

    scatter(grid[0],grid[1]);

    xlabel('tx'); ylabel('ty'); title('Grid Setup')

    scatter(diag(grid[0]),diag(grid[1]),color='r'); # Diagonal shows the N-S line (N is the upper left point)

    xlim([mig,mag]); ylim([mig,mag]);

    text(-dci*ns - dci/2,dci*ns, 'N', fontsize=12)

    show()

    return

def CovFigGrid(tx,ty,dr,Ne,h):


    """
        This definition plots the covariance values for all the tx-ty pairs defined by definition 'grid_2D'.
        N in figure shows the north direction.
    """
    gran = tx.shape[0];    # grid sample range 

    sf = 0;                # counter of subfigure numbers

    fig = figure(figsize = (20,15))

    fig.subplots_adjust(hspace=.7,wspace = 0.3)

    for ny in range(gran):
    
        for nx in range(gran):

            ax = fig.add_subplot(gran, gran, sf)
        
            px, py, uB, Bo, theta, uphi, uthe = geom(tx[nx],ty[ny],dr)
        
            covXX, covXY = covs(px, py, uB, Bo, theta, uphi, uthe, Ne, h)
        
            ax.plot(h,covXX); ax.plot(h,covXY,color='r')

            ax.set_xticks((0,500,1000)); 
        
            ax.set_yticks((0.000,0.015,0.035));  
        
            sf +=1
        
            if nx != 1:
            
                setp(ax.get_yticklabels(), visible=False)

    show()

    return


if __name__=='__main__':

    grid,tx,ty = grid_2D(0.0245,2)

    FigGrid(grid,0.025,2)

    dr = 5.
    h  = arange(200)*dr
    Ne=8.0e12*exp(-h/180.-50.*exp(-h/60.))

    CovFigGrid(tx,ty,5.,Ne,h)

#    h  = arange(200)*dr
    
#    Ne=8.0e12*exp(-h/180.-50.*exp(-h/60.))
