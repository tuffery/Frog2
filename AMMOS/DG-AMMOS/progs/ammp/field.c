
/* AMMP version of finite difference support
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ammp.h"
#include "numeric.h"
/* fd_field is the driver of the whole thing */
int fd_field(op1,op2,D,Z,T,guard,dr,deby)
FILE *op1,*op2; /* output files if needed */
float D,Z,T,guard,dr; /* dielectric, Ion concentration, Temperature, guard radius, grid size */
float deby;  /* radius for D and Z effects */
{
    ATOM *ap,*a_next();
    int numatm, a_number();
    float xmin,xmax,xbar;
    float ymin,ymax,ybar;
    float zmin,zmax,zbar;
    float (*grid)[],(*dielectric)[],(*induced)[],(*field)[];
    int nx,ny,nz,i;
    float dx,dy,dz, din,dout;
    int use_debye, use_dielec;
    float fd_sweep();

    use_debye = 1==0;
    use_dielec = 1==0;
    if( Z > 0) use_debye= 1==1;
    if( D >= 1.) use_dielec = 1==1;
    din = 1.;
    dout = D; if( !use_dielec) dout = 10.;
    if( deby > guard*0.5) deby = guard*0.5;

    if( op1 == NULL ) op1 = stdout;
    /* if( op2 == NULL ) op2 = stdout; */
    grid = NULL;
    dielectric = NULL;
    induced = NULL;

    numatm = a_number();
    if( numatm < 1) return 1==0;

    xmin = 10.e10; ymin = 10.e10; zmin = 10.e10;
    xmax = -xmin; ymax = -ymin; zmax = -zmin;

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->x > xmax) xmax = ap->x;
        if( ap->y > ymax) ymax = ap->y;
        if( ap->z > zmax) zmax = ap->z;
        if( ap->x < xmin) xmin = ap->x;
        if( ap->y < ymin) ymin = ap->y;
        if( ap->z < zmin) zmin = ap->z;
    }
    xmax += guard; ymax += guard; zmax += guard;
    xmin -= guard; ymin -= guard; zmin -= guard;

    xbar =  0.5*(xmax+xmin);
    ybar =  0.5*(ymax+ymin);
    zbar =  0.5*(zmax+zmin);

    nx = (xmax-xmin)/dr + 1;
    ny = (ymax-ymin)/dr + 1;
    nz = (zmax-zmin)/dr + 1;
    grid = malloc( nx*ny*nz*sizeof(float));
    field = malloc( nx*ny*nz*sizeof(float));
    if( grid == NULL || field == NULL)
    {aaerror("cannot allocate memory in fd_field");
        if(field!= NULL) free(field);
        if(grid != NULL) free(grid);
        return 1==0;}
    if( Z > 0.)
    {
        dielectric = malloc( nx*ny*nz*sizeof(float));
        induced  = malloc( nx*ny*nz*sizeof(float));
        if( dielectric == NULL || induced == NULL)
        {aaerror("cannot allocate memory in fd_field");
            if(induced  != NULL) free( induced);
            if( dielectric != NULL) free( dielectric);
            if(field!= NULL) free(field);
            if(grid != NULL) free(grid);
            return 1==0;}
    } else if( D >=1.) {
        dielectric = malloc( nx*ny*nz*sizeof(float));
        if( dielectric == NULL )
        {aaerror("cannot allocate memory in fd_field");
            if( dielectric != NULL) free( dielectric);
            if(field!= NULL) free(field);
            if(grid != NULL) free(grid);
            return 1==0;}
    }
    /* now we have memory, let's set up the data */
    dx = dr; dy = dr; dz = dr;

    for(  i=0; i< nx*ny*nz; i++)
    {
        (*field)[i] = 0.;
    }

    fd_build_charge_map( grid,nx,ny,nz,dx,dy,dz,-xmin,-ymin,-zmin);
    if( use_debye)
    {
        fd_build_dielectric( dielectric,nx,ny,nz,dx,dy,dz,-xmin,-ymin,-zmin,deby,din,dout);
        for(  i=0; i< nx*ny*nz; i++)
        {
            (*induced)[i] = 0.;
        }
    } else if( use_dielec)
    {
        fd_build_dielectric( dielectric,nx,ny,nz,dx,dy,dz,-xmin,-ymin,-zmin,deby,din,dout);
    }
    /* now do the work (fixed iteration number, is this bogus?) */
    for( i=0; i< 1000; i++)
    {
        if( use_dielec && use_debye)
            /*float fd_sweep(grid,induced,dielectric,field, nx,ny,nz,dx,dy,dz) */
            fd_sweep( grid,induced,dielectric,field,nx,ny,nz,dx,dy,dz);
        else if( use_dielec )
            fd_sweep( grid,NULL,dielectric,field,nx,ny,nz,dx,dy,dz);
        else if( use_debye)
            fd_sweep( grid,induced,NULL,field,nx,ny,nz,dx,dy,dz);
        else
            fd_sweep( grid,NULL,NULL,field,nx,ny,nz,dx,dy,dz);

        if( i%5==4 &&  use_debye)
            fd_ionize( field,dielectric,induced,nx,ny,nz,din,dout,Z,T*1.987e-3);
    }
    /* output the maps */
    /*int dump_xmap( where,what,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz) */
    dump_xmap(op1,field,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
    if( use_debye && op2 != NULL)
        dump_xmap(op2,induced,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
    /* for debugging */
    /*
    	if( !use_debye && op2 != NULL)
    	dump_xmap(op2,grid,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz);
    */
    fd_interpolate( stdout,field, nx,ny,nz,dx,dy,dz,xmin,ymin,zmin);

    /* clean up of the memory we may use */
    if( induced != NULL) free(induced);
    if( dielectric != NULL) free(dielectric);
    free(field);
    free(grid);
    return 1==1;
}/* end of fd_field */

int fd_ionize( grid,dielectric,induced,nx,ny,nz,in,out, Z,kt)
float (*grid)[],(*dielectric)[],(*induced)[];
float in ,out,Z,kt;  /* inner, outer dielectric, Ionic strength, kT*/
{
    int i;
    float x,xmax,xmin;
    kt = 1./kt;
    xmax = 0.;
    xmin = 0.;

    for( i=0; i< nx*ny*nz; i++)
    {
        if( (*dielectric)[i] != in)
        {
#define linear
#ifdef linear
            x  = -(*grid)[i]*Z*kt; /* magically now in kcal/mol */
            (*induced)[i] = 2.*x/(*dielectric)[i]*4*PI ;
#else
            x  = -(*grid)[i]*332.17752*Z*kt*4*PI; /* magically now in kcal/mol */

            x = Z*( exp(x) -exp(-x));
            (*induced)[i] = x;
#endif

            if( x > xmax) xmax = x;
            if( x < xmin) xmin = x;

        }
    }
    printf("%f %f\n",xmax,xmin);
    return 1==1;
}/* end of fd_ionize */

int fd_build_dielectric( grid,nx,ny,nz,dx,dy,dz,ox,oy,oz,rad,in,out)
float (*grid)[],dx,dy,dz,ox,oy,oz;
float rad,in,out; /* radius, dielectric within radius, dielectric outside*/
int nx,ny,nz;
{

    int i,j,k,inxy;
    int ix,iy,iz;
    int irx,iry,irz;
    int jrx,jry,jrz;
    int ii,jj,kk;
    float fx,fy,fz;
    float r2, rx,ry,rz;
    ATOM *ap,*a_next();
    int iatm,numatm, a_number();
    int (*mask)[];

    numatm = a_number();
    if( numatm < 1) return 1==0 ;
    inxy = nx*ny;

    fx = 1./dx; fy = 1./dy; fz = 1./dz;
    irx = rad*fx; iry = rad*fy; irz = rad*fz;
    jrx = 2*irx+1; jry = 2*iry + 1; jrz = 2*irz + 1;
    mask  = malloc( jrx*jry*jrz*sizeof(int));
    if( mask == NULL)
    {aaerror("cannot allocate space in fd_build_dielectric"); return 1==0;}

    /*
    	in = 1./in;
    	out = 1./out;
    */

    r2 = rad*rad;
    for( k=0; k< jrz; k++)
        for( j=0; j< jry; j++)
            for( i=0; i< jrx; i++)
            {
                rx = (float)(i- irx)*dx;
                ry = (float)(j- iry)*dy;
                rz = (float)(k- irz)*dz;
                rx = rx*rx + ry*ry + rz*rz;
                if( rx < rad)
                    (*mask)[(k*jry+j)*jrx + i] = 1==1;
                else
                    (*mask)[(k*jry+j)*jrx + i] = 1==0;
            }


    for( i=0; i< nx*ny*nz; i++)
        (*grid)[i] = out;
    for( iatm = 0; iatm < numatm; iatm++)
    {  	ap = a_next(iatm);
        ix = fx*(ap->x + ox) +0.5;
        iy = fy*(ap->y + oy) +0.5;
        iz = fz*(ap->z + oz) +0.5;
        if( ix > nx-irx) {free(mask); return 1==0;}
        if( iy > ny-irx) {free(mask); return 1==0;}
        if( iz > nz-irx) {free(mask); return 1==0;}
        if( ix < irx) {free(mask); return 1==0;}
        if( iy < iry) {free(mask); return 1==0;}
        if( iz < irz) {free(mask); return 1==0;}
        for(k=0; k< jrz; k++)
        {
            kk = k+iz-irz;
            for(j=0; j< jry; j++)
            {
                jj = j+iy - iry;
                for(i=0; i< jrx; i++)
                {
                    ii = i+ix -irx;
                    ii = kk*inxy+jj*nx + ii;
                    if((*mask)[(k*jry+j)*jrx+i]) (*grid)[ii] = in;
                } } }
    }
    free(mask);
    return 1==1;
}/* end of fd_build_dielectric */

int fd_build_charge_map( grid,nx,ny,nz,dx,dy,dz,ox,oy,oz)
float (*grid)[],dx,dy,dz,ox,oy,oz;
int nx,ny,nz;
{

    int i,j,k,inxy;
    int ix,iy,iz;
    float fx,fy,fz;
    float dv;
    ATOM *ap,*a_next();
    int iatm,numatm, a_number();

    numatm = a_number();
    if( numatm < 1) return 1==0 ;
    inxy = nx*ny;

    fx = 1./dx; fy = 1./dy; fz = 1./dz;
    dv = fx*fy*fz;
    for( i=0; i< nx*ny*nz; i++)
        (*grid)[i] = 0.;

    for( iatm = 0; iatm < numatm; iatm++)
    {  	ap = a_next(iatm);
        if( ap->active){
            ix = fx*(ap->x + ox) +0.5;
            iy = fy*(ap->y + oy) +0.5;
            iz = fz*(ap->z + oz) +0.5;
            if( ix > nx) return 1==0;
            if( iy > ny) return 1==0;
            if( iz > nz) return 1==0;
            if( ix < 0) return 1==0;
            if( iy < 0) return 1==0;
            if( iz < 0) return 1==0;
            i = iz*inxy + iy*nx + ix;
            (*grid)[i] += ap->q*dv;
        }
    }
    return 1==1;
}/* end of fd_build_charge_map */

/* fd_sweep
*  finite difference sweep program
*
* Use Gauss's theorm to solve the Possion equations
*  this is one sweep and you call it a whole bunch of times to get the
*  the feild
*
*  solve the differential equation (del) dot D = 4 pi Rho
*  where D = Dielectric (Del) phi
*
*  but use integration shells
*  much more numerically stable and ALMOST finite elements
*
* if induced or dielectric is NULL the effects are ignored
*
*/
float fd_sweep(grid,induced,dielectric,field, nx,ny,nz,dx,dy,dz)
float (*grid)[]; /* the charges*/
float (*induced)[]; /* induced charges*/
float (*dielectric)[]; /* 3-dimensional map of the dielectric */
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
    /* the next lines give the numerical integral of the flux from an element
    *   of unit charge
    */
    dv = 2.*(dx*dy + dy*dz + dx*dz);
    dv = 1./dv;
    dv = dv/(4.*PI);
    inxy = nx*ny;
    deltamax = 0.;

    if( induced == NULL && dielectric == NULL){
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
        }/* k */} /* if not done */
    if( induced != NULL && dielectric == NULL){
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
                    delta = (*grid)[ixp] - flux*dv + (*induced)[ixp];
                    if( fabs(delta) > deltamax) deltamax = fabs(delta);
                    (*field)[ixp] += delta*0.5;
                }/* i */
            }/* j */
        }/* k */}
    if( induced == NULL && dielectric != NULL){
        for( k=1; k< nz-1; k++)
        {
            izp = k*inxy;
            for( j=1; j< ny-1; j++)
            {
                iyp  = j*nx;
                for( i=1; i< nx-1; i++)
                {
                    ixp = i+iyp + izp;
                    flux = ((*field)[ixp]*(*dielectric)[ixp]
                            -(*field)[ixp-1]*(*dielectric)[ixp-1])*flx;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+1]*(*dielectric)[ixp+1])*flx;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp-nx]*(*dielectric)[ixp-nx])*fly;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+nx]*(*dielectric)[ixp+nx])*fly;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp-inxy]*(*dielectric)[ixp-inxy])*flz;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+inxy]*(*dielectric)[ixp+inxy])*flz;
                    delta = (*grid)[ixp] - flux*dv;
                    delta /= (*dielectric)[ixp];
                    if( fabs(delta) > deltamax) deltamax = fabs(delta);
                    (*field)[ixp] += delta*0.5;
                }/* i */
            }/* j */
        }/* k */}

    if( induced != NULL && dielectric != NULL){
        for( k=1; k< nz-1; k++)
        {
            izp = k*inxy;
            for( j=1; j< ny-1; j++)
            {
                iyp  = j*nx;
                for( i=1; i< nx-1; i++)
                {
                    ixp = i+iyp + izp;
                    flux = ((*field)[ixp]*(*dielectric)[ixp]
                            -(*field)[ixp-1]*(*dielectric)[ixp-1])*flx;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+1]*(*dielectric)[ixp+1])*flx;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp-nx]*(*dielectric)[ixp-nx])*fly;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+nx]*(*dielectric)[ixp+nx])*fly;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp-inxy]*(*dielectric)[ixp-inxy])*flz;
                    flux += ((*field)[ixp]*(*dielectric)[ixp]
                             -(*field)[ixp+inxy]*(*dielectric)[ixp+inxy])*flz;
                    delta = (*grid)[ixp] +(*induced)[ixp] - flux*dv;
                    delta /= (*dielectric)[ixp];
                    if( fabs(delta) > deltamax) deltamax = fabs(delta);
                    (*field)[ixp] += delta*0.5;
                }/* i */
            }/* j */
        }/* k */}


    return deltamax;
}/* end of fd_sweep */



/* map.c
*
* routine to  output XPLOR MAPS
*
*/
/*
*  copyright 1993,1994,1995 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/


/* do the map */

int dump_xmap( where,what,xmin,xmax,ymin,ymax,zmin,zmax,nx,ny,nz)
FILE *where;
float (*what)[];
float xmin,xmax,ymin,ymax,zmin,zmax;
int nx,ny,nz;
{
    int irow;
    int imin,imax,jmin,jmax,kmin,kmax;
    int i,j,k;

    imin = xmin*nx/(xmax-xmin); imax = (xmax*nx/(xmax-xmin));
    jmin = ymin*ny/(ymax-ymin); jmax = (ymax*ny/(ymax-ymin));
    kmin = zmin*nz/(zmax-zmin); kmax = (zmax*nz/(zmax-zmin));
    if( imin >= 0) imin += 1;
    if( jmin >= 0) jmin += 1;
    if( kmin >= 0) kmin += 1;

    /* write the header */
    fprintf(where,"\n       2 !NTITLE\nREMARKS AMMP xplor map                                                                                 \nREMARKS                                                                                                 \n");
    fprintf(where,"%8i%8i%8i%8i%8i%8i%8i%8i%8i\n",
            nx,imin,imax,ny,jmin,jmax,nz,kmin,kmax);
    fprintf(where,"%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\nZXY",
            xmax-xmin,ymax-ymin,zmax-zmin,90.,90.,90.);
    /* now do the work we loop with z outer then x then y
    * and fake a FORTRAN implied do */
    for( k=0; k< nz; k++)
    {
        fprintf(where,"\n%8d\n",k);
        irow = 0;
        for( j=0; j< ny; j++)
            for( i=0; i< nx; i++)
            {
                if(irow ==6) {fprintf(where,"\n"); irow = 0;}
                fprintf(where,"%12.5e",(*what)[(k*ny+j)*nx + i]);
                irow ++;
            }}
}/* end of routine */

/*
	fd_interpolate( stdout,field, nx,ny,nz,dx,dy,dz,xmin,ymin,zmin);
*/
int fd_interpolate( op, what,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin)
FILE *op;
float (*what)[];
int nx,ny,nz; /* map size */
float dx,dy,dz; /* map grid steps */
float xmin,ymin,zmin; /* origin of the map in Angstroms */
{
    ATOM *ap, *a_next();
    int numatom, a_number();
    int i,j,k;
    int ix,iy,iz;
    int iatom;
    float x,y,z;
    float xx,yy;
    float gx,gy,gz;
    float tensor[2][2][2];
    float accum;

    numatom = a_number();
    if( numatom == 0 ) return 0;

    for( iatom = 0; iatom < numatom ; iatom ++)
    {
        ap = a_next(iatom);
        if( ap->active){/* skip me if i'm active */}
        else{ /* if i'm inactive you need to interpolate */
            x = ap->x;
            y = ap->y;
            z = ap->z;
            gx = (x-xmin)/dx;
            gy = (y-ymin)/dy;
            gz = (z-zmin)/dz;
            ix = gx; /* should round down */
            iy = gy;
            iz = gz;

            gx = gx - ix;
            gy = gy - iy;
            gz = gz - iz;
            /* now gx should be the fractional offset */
            accum = 0.;
            i = ( iz*ny + iy)*nx + ix;
            /* test on x only */
            /*
            		x = (*what)[i+1]*gx + (*what)[i]*(1.-gx);
            		y = (*what)[i+nx+1]*gx + (*what)[i+nx]*(1.-gx);
            		xx = (*what)[i+nx*ny+1]*gx + (*what)[i]*(1.-gx);
            		yy = (*what)[i+nx*ny+1]*gx + (*what)[i+nx*ny]*(1.-gx);
            		x = y*gy +x*(1.-gy);
            		y = yy*gy + xx*(1.-gy);
            		accum = y*gz + x*(1.-gz);
            */
            tensor[0][0][0] = (*what)[i];
            tensor[1][0][0] = (*what)[i+1];
            tensor[0][1][0] = (*what)[i+nx];
            tensor[1][1][0] = (*what)[i+1+nx];
            tensor[0][0][1] = (*what)[i+nx*ny];
            tensor[1][0][1] = (*what)[i+1+nx*ny];
            tensor[0][1][1] = (*what)[i+nx*ny+nx];
            tensor[1][1][1] = (*what)[i+1+nx*ny+nx];
     /*       printf("%f \n", tensor[0][0][0]);
*/

            tensor[0][0][0] = tensor[0][0][1]*gz + tensor[0][0][0]*(1.-gz);
            tensor[1][0][0] = tensor[1][0][1]*gz + tensor[1][0][0]*(1.-gz);
            tensor[0][1][0] = tensor[0][1][1]*gz + tensor[0][1][0]*(1.-gz);
            tensor[1][1][0] = tensor[1][1][1]*gz + tensor[1][1][0]*(1.-gz);

            tensor[0][0][0] = tensor[0][1][0]*gy + tensor[0][0][0]*(1.-gy);
            tensor[1][0][0] = tensor[1][1][0]*gy + tensor[1][0][0]*(1.-gy);

            accum = tensor[1][0][0]*gx + tensor[0][0][0]*(1.-gx);
            fprintf(op," Atom %d at %f %f %f interpolated %f %f \n",
                    ap->serial, ap->x, ap->y,ap->z,accum,accum*331.17752);

        }
    }/* iatom */

    fflush(op);

    return 0;
} /* end of fd_interpolate */

