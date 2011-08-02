/* tether.c
*
* collection of routines to service tether  potentials
*  these fix an atom to a location 
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of tethers, with the associated position and force constant
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
*/
/*
*  copyright 1992 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-tethered potentials
*/
typedef struct{
    ATOM *atom1;
    float k,x,y,z;
    void *next;
}  TETHER;
#define tetherLONG sizeof(TETHER)

TETHER *tether_first = NULL;
TETHER *tether_last = NULL;
/* function tether adds a tether to the tether list
* returns 1 if ok
* returns 0 if not
*  is passed the atom serial number, position and constant
* allocates the new memory, initializes it and
* returns
*/
int tether( p1,fk,x,y,z)
int p1;
float fk ,x,y,z;
{
    ATOM *ap1,*a_m_serial();
    TETHER *new;
    char line[BUFSIZ];
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    if( (ap1 == NULL)  )
    {
        sprintf( line,"undefined atom in tether %d \0",p1);
        aaerror( line );
        return 0;
    }
    /* check to see if it exists and update if so */
    if( tether_first != NULL)
    {
        new = tether_first;
        while(1)
        {
            if( new ->next == NULL) break;
            if( new->atom1 == ap1)
            {
                new -> k = fk;
                new -> x = x;
                new -> y = y;
                new -> z = z;
                return 1;
            }
            if( new->next == new ) break;
            new = new->next;
        }
    }
    if( ( new = malloc( tetherLONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( tether_first == NULL) tether_first = new;
    if( tether_last == NULL) tether_last = new;
    new -> atom1 = ap1;
    new -> k = fk;
    new -> x = x;
    new -> y = y;
    new -> z = z;
    new -> next = new;
    tether_last -> next = new;
    tether_last = new;
    return 1;
}

int inactive_tether()
{
    TETHER *tp;
    ATOM *ap;
    tp = tether_first;
    if( tp == NULL) return 0;

    while(1==1)
    {
        ap = tp->atom1;
        if( ap != NULL) ap->active = (1==0);
        if( tp == tp->next) break;
        tp = tp->next;
        if( tp == NULL) break;
    }
    return 0;
}

int alltether( fk)
float fk ;
{
    ATOM *ap1 ;
    ATOM *a_next();
    int i,numatm,a_number();
    numatm = a_number();
    if( numatm < 1 ) return 0;

    /* */
    /* get the first atom and then loop over them all */
    /* write a tether for each with the given force constant */
    for( i=0; i< numatm; i++)
    {
        ap1 = a_next(i);
        if( fabs(ap1->x ) > 1.e-8)
            if( fabs(ap1->y ) > 1.e-8)
                if( fabs(ap1->z ) > 1.e-8)
                {
                    if( tether( ap1->serial, fk, ap1->x,ap1->y,ap1->z) )
                    { } else {return 0 ;}
                }
    }
    return numatm;
}


/* v_tether()
* this function sums up the potentials
* for the atoms defined in the TETHER data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_tether( V, lambda )
float *V,lambda;
{
    TETHER *bp;
    float r,xt,yt,zt;
    ATOM *a1;


    bp = tether_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1;
        if( a1->active ){
            if( lambda == 0.)
            {
                r =     (a1->x - bp->x)*(a1->x - bp->x);
                r = r + (a1->y - bp->y)*(a1->y - bp->y);
                r = r + (a1->z - bp->z)*(a1->z - bp->z);
            } else
            {
                xt = (a1->x -bp->x +lambda*(a1->dx));
                yt = (a1->y -bp->y +lambda*(a1->dy));
                zt = (a1->z -bp->z +lambda*(a1->dz));
                r = xt*xt+yt*yt+zt*zt;
            }
            *V += bp->k*r ;
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_tether()
*
* f_tether increments the forces in the atom structures by the force
* due to the tether components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_tether(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    TETHER *bp;
    float r,k,ux,uy,uz;
    ATOM *a1;


    bp = tether_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1;
        if( a1->active){
            if( lambda == 0.)
            {
                ux = (a1->x - bp->x );
                uy = (a1->y - bp->y );
                uz = (a1->z - bp->z );
            }else{
                ux = (a1->x -bp->x +lambda*(a1->dx));
                uy = (a1->y -bp->y +lambda*(a1->dy));
                uz = (a1->z -bp->z +lambda*(a1->dz));
            }
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { goto SKIP; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            ux = 2*k*(r)*ux;
            uy = 2*k*(r)*uy;
            uz = 2*k*(r)*uz;
            a1->fx -= ux;
            a1->fy -= uy;
            a1->fz -= uz;
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* routine dump_tethers
* this function outputs the tether parameters
* and does it in a simple form
* tether ser1,k,x,y,z
* the rest is just free format
*/
void dump_tethers( where )
FILE *where;
{
    TETHER *b;
    ATOM *a1;
    b = tether_first;
    if( b == NULL ) return;
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1;
        fprintf( where,"tether %d %f %f %f %f ;\n",a1->serial,
                 b->k,b->x,b->y,b->z);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1;
    fprintf( where,"tether %d %f %f %f %f ;\n",a1->serial,
             b->k,b->x,b->y,b->z);
}
/* a_tether()
* this function sums up the potentials
* for the atoms defined in the TETHER data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_tether( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    TETHER *bp;
    float r,xt,yt,zt;
    float rms, rmax;
    int tried,imax;
    ATOM *a1;
    float bstrot();
    float (*x)[],(*y)[],(*z)[],(*xx)[],(*yy)[],(*zz)[];
    float matrix[3][3],delta[3];
    int numatm,a_number();


    rms = 0.; rmax = -1.;
    tried = 0;
    bp = tether_first;
    if( bp == NULL ) return 1;
    numatm = a_number();
    x = malloc( numatm*sizeof(float));
    y = malloc( numatm*sizeof(float));
    z = malloc( numatm*sizeof(float));
    xx = malloc( numatm*sizeof(float));
    yy = malloc( numatm*sizeof(float));
    zz = malloc( numatm*sizeof(float));

    if( x == NULL || y == NULL || z== NULL || xx == NULL ||
            yy == NULL || zz== NULL)
    {aaerror("cannot allocate memory"); return 1;}

    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1;
        if( a1->serial >= ilow && a1->serial <= ihigh)
        {
            if( lambda == 0.)
            {
                r =     (a1->x - bp->x)*(a1->x - bp->x);
                r = r + (a1->y - bp->y)*(a1->y - bp->y);
                r = r + (a1->z - bp->z)*(a1->z - bp->z);
            } else
            {
                xt = (a1->x -bp->x +lambda*(a1->dx));
                yt = (a1->y -bp->y +lambda*(a1->dy));
                zt = (a1->z -bp->z +lambda*(a1->dz));
                r = xt*xt+yt*yt+zt*zt;
            }

            (*x)[tried] = a1->x+lambda*a1->dx;
            (*y)[tried] = a1->y+lambda*a1->dy;
            (*z)[tried] = a1->z+lambda*a1->dz;
            (*xx)[tried] = bp->x;
            (*yy)[tried] = bp->y;
            (*zz)[tried] = bp->z;

            tried += 1;
            rms = rms + r;
            if( r > rmax ) {rmax = r; imax = a1->serial;}

            zt = bp->k*r ;
            *V += zt;
            fprintf(op,"Tether %d E %f error %f\n",a1->serial,zt,sqrt(r));
        }
        if( bp == bp->next ) break;
        bp = bp->next;
    }
    if( tried > 0){
        rms = sqrt(rms/tried); rmax = sqrt(rmax);
        fprintf(op," RMSD %f Maximum Deviation %f on atom %d\n",rms,rmax,imax);
        fprintf(op," RMSD after superposition %f\n",
                bstrot( &(*x)[0], &(*y)[0], &(*z)[0], &(*xx)[0], &(*yy)[0],
                        &(*zz)[0],tried,matrix,delta) );
    }
    free(zz); free(yy); free(xx); free(z); free(y); free(x);
}
/* routine to superimpose two sets of coordinates
*  
* cadged from Larry Andrews' Newvector.f library
*
C  The algorithm is from Ferro and Hermans, Acta Cryst., A33,345-347
C  (1977).  Other algorithms give the same result, such as:
C
C  Kearsley, Acta Cryst., A45, 208-210 (1989)
C  Diamond,  Acta Cryst., A44, 211-216 (1988)
C  Kabsch,   Acta Cryst., A34, 827-828 (1978)
C  Kabsch,   Acta Cryst., A32, 922-923 (1976)
C
C  The algorithm of Ferro and Hermans has the advantage that the
C  rotational components are removed iteratively, so that residual
C  numerical errors are correctly removed.  This leaves no
C  possibility of poorly orthogonalized matrices.  The explicit
C  iteration is not a problem.  The other methods simply hide the
C  iteration in the eigenvector solver.  Ferro and Hermans
C  simply build what is a simple, special purpose version of the QR
C  method.
*/

#include <math.h>
/*void matmul(float[],float[],float[],int,int);
*/

/* return the RMS error, -1 on internal error */
float bstrot( x,y,z,xx,yy,zz, na, matrix,delta)
int na;
float x[],y[],z[];
float xx[],yy[],zz[];
float matrix[3][3], delta[3];
{

    float tensor[3][3],cx,cy,cz,cxx,cyy,czz;
    float tx,ty,tz,txx,tyy,tzz;
    float rms;
    float sx[3][3],sy[3][3],sz[3][3];
    float sq[3][3];
    void cpyvec();
    void matmul();
    int i,j,ipass;


    if( na < 1) return -1.;
    /* find the centers of mass */
    cx = 0.;
    cy = 0.;
    cz = 0.;
    cxx = 0.;
    cyy = 0.;
    czz = 0.;
    for( i=0; i< na; i++)
    {
        cx += x[i];
        cy += y[i];
        cz += z[i];
        cxx += xx[i];
        cyy += yy[i];
        czz += zz[i];
    }
    cx /= na;
    cy /= na;
    cz /= na;
    cxx /= na;
    cyy /= na;
    czz /= na;
    /* make the metric tensor */
    for( i=0; i< 3; i++)
        for( j=0; j< 3; j++)
        {
            tensor[i][j] = 0.;
            matrix[i][j] = 0.;
            sx[i][j] = 0.;
            sy[i][j] = 0.;
            sz[i][j] = 0.;
        }
    matrix[0][0] = 1.;
    matrix[1][1] = 1.;
    matrix[2][2] = 1.;
    sx[0][0] = 1.;
    sy[1][1] = 1.;
    sz[2][2] = 1.;
    for( i=0; i<na; i++)
    {
        tx = x[i] - cx;
        ty = y[i] - cy;
        tz = z[i] - cz;
        txx = xx[i] - cxx;
        tyy = yy[i] - cyy;
        tzz = zz[i] - czz;
        tensor[0][0] += tx*txx;
        tensor[0][1] += tx*tyy;
        tensor[0][2] += tx*tzz;
        tensor[1][0] += ty*txx;
        tensor[1][1] += ty*tyy;
        tensor[1][2] += ty*tzz;
        tensor[2][0] += tz*txx;
        tensor[2][1] += tz*tyy;
        tensor[2][2] += tz*tzz;
    }
    /* now find the linear orthogonal transformation which symetrizes
    * the metric tensor */

    for( ipass = 0; ipass < 20; ipass ++)
    {
        rms = 0.;
        /* x */
        tx = atan2( tensor[2][1]-tensor[1][2],
                    tensor[1][1]+tensor[2][2]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sx[1][1] = ty;
        sx[2][1] = -tz;
        sx[1][2] = tz;
        sx[2][2] = ty;
        matmul( sx,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sx,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /* y */
        tx = atan2( tensor[2][0]-tensor[0][2],
                    tensor[0][0]+tensor[2][2]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sy[0][0] = ty;
        sy[2][0] = -tz;
        sy[0][2] = tz;
        sy[2][2] = ty;
        matmul( sy,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sy,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /* z */
        tx = atan2( tensor[0][1]-tensor[1][0],
                    tensor[1][1]+tensor[0][0]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sz[1][1] = ty;
        sz[0][1] = -tz;
        sz[1][0] = tz;
        sz[0][0] = ty;
        matmul( sz,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sz,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /*  termination critereon here */
        if( rms < 1.e-7) break;
    }
    rms = 0.;
    for(i=0; i< na; i++)
    {
        txx = xx[i] - cxx;
        tyy = yy[i] - cyy;
        tzz = zz[i] - czz;
        tx = matrix[0][0]*txx + matrix[1][0]*tyy + matrix[2][0]*tzz;
        ty = matrix[0][1]*txx + matrix[1][1]*tyy + matrix[2][1]*tzz;
        tz = matrix[0][2]*txx + matrix[1][2]*tyy + matrix[2][2]*tzz;
        tx += cx - x[i];
        ty += cy - y[i];
        tz += cz - z[i];
        rms += tx*tx + ty*ty + tz*tz;
    }
    tx = matrix[0][0]*cxx + matrix[1][0]*cyy + matrix[2][0]*czz;
    ty = matrix[0][1]*cxx + matrix[1][1]*cyy + matrix[2][1]*czz;
    tz = matrix[0][2]*cxx + matrix[1][2]*cyy + matrix[2][2]*czz;
    delta[0] = cx - tx;
    delta[1] = cy - ty;
    delta[2] = cz - tz;
    return sqrt(rms/na);


}/* end of routine */
/* copy a vector into another */
void cpyvec(orig,copy,n)
float orig[],copy[];
int n;
{
    int i;
    for( i=0; i< n; i++)
        copy[i] = orig[i];
}
/*
multiply c[n][n] =  a[n][m] b[m][n];
*/
void matmul( a,b,c,n,m)
float a[],b[],c[];
int n,m;
{
    int i,j,k,ioff,koff;

    for( i=0; i< n*n; i++)
        c[i] = 0.;
    for( i=0; i< n; i++)
    {
        ioff = i*n;
        for( j=0; j< n; j++)
        {
            koff = 0.;
            for( k=0; k<m; k++)
            {
                c[ ioff +j] += a[ioff + k] *b[ j +koff];
                koff += m;
            }
        }
    }

}
/* v_ho_tether()
* this function sums up the potentials
* for the atoms defined in the TETHER data structure.
* homotropy version !!!
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_ho_tether( V, lambda )
float *V,lambda;
{
    TETHER *bp;
    float r,xt,yt,zt;
    ATOM *a1;
    float hol, get_f_variable();

    hol = get_f_variable( "lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;
    if( hol == 1.) return;

    bp = tether_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1;
        if( a1->active ){
            if( lambda == 0.)
            {
                xt = (a1->x -bp->x );
                yt = (a1->y -bp->y );
                zt = (a1->z -bp->z );
            } else
            {
                xt = (a1->x + lambda*(a1->dx) -bp->x);
                yt = (a1->y + lambda*(a1->dy) -bp->y);
                zt = (a1->z + lambda*(a1->dz) -bp->z);
            }
            r =(xt*xt+yt*yt+zt*zt)*(1.-hol)*(1.-hol);
            *V += bp->k*r ;
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* f_ho_tether()
*
* f_tether increments the forces in the atom structures by the force
* due to the tether components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*
*homotropy version
*/
int f_ho_tether(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    TETHER *bp;
    float r,k,ux,uy,uz;
    ATOM *a1;
    float hol, get_f_variable();

    hol = get_f_variable( "lambda");
    if( hol < 0. ) hol = 0.;
    if( hol > 1. ) hol = 1.;
    if( hol == 1.) return;



    bp = tether_first;
    if( bp == NULL ) return 1;
    while(1)
    {
        if( bp == NULL) return 0;
        k = bp->k;
        a1 = bp->atom1;
        if( a1->active){
            if( lambda == 0.)
            {
                ux = (a1->x - bp->x );
                uy = (a1->y - bp->y );
                uz = (a1->z - bp->z );
            }else{
                ux = (a1->x -bp->x +lambda*(a1->dx));
                uy = (a1->y -bp->y +lambda*(a1->dy));
                uz = (a1->z -bp->z +lambda*(a1->dz));
            }
            r = ux*ux + uy*uy + uz*uz;
            /* watch for FP errors*/
            if( r <= 1.e-5)
            { goto SKIP; }else{
                r = sqrt(r); ux = ux/r; uy = uy/r; uz = uz/r;
            }
            ux = 2*k*(r)*ux*(1.-hol);
            uy = 2*k*(r)*uy*(1.-hol);
            uz = 2*k*(r)*uz*(1.-hol);
            a1->fx -= ux;
            a1->fy -= uy;
            a1->fz -= uz;
        }
SKIP:
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
