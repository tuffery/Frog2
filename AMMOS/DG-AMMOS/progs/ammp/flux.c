/* this is the initial prototype of a solver for pb-electrostatics
*  the eventual idea is to output the interesting parts
* of the  solution as well as the field.
*
*  i.e. we need maps of the debye charge density
*       and of the "dielectric" charge.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "numeric.h"

main()
{
    float fd_sweep();
    float (*grid)[],(*field)[];
    float delta;
    int i,j,k;
    int nx,ny,nz;

    nx = 10; ny = 10; nz = 10;
    grid = malloc( nx*ny*nz*sizeof(float));
    field = malloc(nx*ny*nz*sizeof(float));

    (*grid)[4*nx*ny + 3*nx +3] = 0.1;
    (*grid)[4*nx*ny + 7*nx +7] = -0.1;

    for( k=0; k< 80; k++)
    {
        delta = fd_sweep(grid,field,nx,ny,nz,1.0,1.0,1.0);
        /*
        	delta = fd_sweep(grid,field,nx,ny,nz,0.5,0.5,0.5);  
        	delta = fd_sweep(grid,field,nx,ny,nz,1.0,1.0,1.0); 

        	delta = fd_sweep(grid,field,nx,ny,nz,2.0,2.0,2.0);  
        	delta = fd_sweep(grid,field,nx,ny,nz,4.0,4.0,4.0); 
         */
        if( k == 79){

            printf("%f\n",delta);
            for( i=0; i< nx; i++)
            {
                for( j=0; j<ny; j++)
                {
                    printf("%10.5e ",(*field)[4*nx*ny + j*nx + i]);
                }
                printf("\n");
            }
        }
    }
}/* end of main */

float fd_sweep(grid,field, nx,ny,nz,dx,dy,dz)
float (*grid)[]; /* the charges*/
float (*field)[]; /* the field */
int nx,ny,nz; /* the numbers of grid points */
float dx,dy,dz; /* the size of the grids in Angstroms */
{
    int i,j,k;
    float delta,deltamax;
    float flux;
    float flx,fly,flz,dv;
    int izp,iyp,ixp,inxy;


    flx = dy*dz/dx;
    fly = dx*dz/dy;
    flz = dx*dy/dz;
    dv = 2.*(dx*dy + dy*dz + dx*dz);
    dv = 1./dv;
    dv = dv/(4.*PI);
    inxy = nx*ny;
    deltamax = 0.;

    for( k=1; k< nz-1; k++)
    {
        izp = k*inxy;
        for( j=1; j< ny-1; j++)
        {
            iyp  = j*nx;
            for( i=1; i< nx-1; i++)
            {
                ixp = i+iyp + izp;
                flux = ((*field)[ixp] -(*field)[ixp-1])*flx;
                flux += ((*field)[ixp] -(*field)[ixp+1])*flx;
                flux += ((*field)[ixp] -(*field)[ixp-nx])*fly;
                flux += ((*field)[ixp] -(*field)[ixp+nx])*fly;
                flux += ((*field)[ixp] -(*field)[ixp-inxy])*flz;
                flux += ((*field)[ixp] -(*field)[ixp+inxy])*flz;
                delta = (*grid)[ixp] - flux*dv;
                if( fabs(delta) > deltamax) deltamax = fabs(delta);
                (*field)[ixp] += delta*0.5;
            }/* i */
        }/* j */
    }/* k */

    return deltamax;
}/* end of fd_sweep */


