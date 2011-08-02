/* abc.c
*
* collection of routines to service bond abc  potentials
*
* abc potentials are Angle Bond Correlations
*
* abc changes k and target bond length in the bond potentials
*
*
*   v(bond) = k(theta)( length - target(theta) )**2
*   k(theta) = k0 + dk/dtheta *delta theta  (two dk's per abc )
*  target(theta) = target0 + dr/dtheta *delta theta (one per abc)
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of abc, with the associated length and force constant
*
* (this could be table driven but what the hell memories cheap)
*
* the routines for potential value, force and (eventually) second
* derivatives are here also
*
* force and 2nd derivative routines assume zero'd arrays for output
* this allows for parralellization if needed (on a PC?)
*
* forces are abc wise symmetric - so we don't have to fuck around with
* s matrices and the like.
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
* but otherwise is self-contained. Note the hooks for Non-abced potentials
*/
typedef struct{
    ATOM *atom1,*atom2;
    float length,k;
    float dlength,dk; /* for abc terms */
    int ndiff;
    void *next;
}  BOND;

extern BOND *bond_first;

typedef struct{
    ATOM *atom1,*atom2,*atom3;
    float theta,thetazero,dr,dk1,dk2;
    void *b1,*b2 ; /* pointers to the bonds for atom1 atom 2 */
    void *next;
}  ABC;
#define ALONG sizeof(ABC)

ABC *abc_first = NULL;
ABC *abc_last = NULL;
/* function abc adds a abc to the abc list
* returns 1 if ok
* returns 0 if not
*  is passed the array pointers, length and constant
* allocates the new memory, initializes it and
* returns
*/
int abc( p1,p2,p3,theta,thetazero,dr,dk1,dk2)
int p1,p2,p3;
float theta ,thetazero,dr,dk1,dk2;
{
    ABC *new;
    ATOM *ap1,*ap2,*ap3,*a_m_serial();
    char line[BUFSIZ];
    int i;
    /* get the atom pointers for the two serial numbers */
    ap1 = a_m_serial( p1 );
    ap2 = a_m_serial( p2 );
    ap3 = a_m_serial( p3 );
    if( (ap1 == NULL) || (ap2 == NULL) || (ap3==NULL) )
    {
        sprintf( line,"undefined atom in abc %d %d %d \0",p1,p2,p3);
        aaerror( line );
        return 0;
    }

    if( ( new = malloc( ALONG ) ) == NULL)
    {
        return 0;
    }
    /* initialize the pointers */
    if( abc_first == NULL) abc_first = new;
    if( abc_last == NULL) abc_last = new;
    new -> atom1 = ap1;
    new -> atom2 = ap2;
    new -> atom3 = ap3;
    new -> dr = dr;
    new -> theta = theta;
    new -> thetazero = thetazero;
    new -> dk1 =  dk1;
    new -> dk2 =  dk2;
    new ->b1 = NULL ; new->b2 = NULL ; /* look up on first call */
    new -> next = new;
    abc_last -> next = new;
    abc_last = new;
    return 1;
}


/* v_abc()
* this function sums up the potentials
* for the atoms defined in the abc data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int v_abc( V, lambda )
float *V,lambda;
{
    ABC *abcp;
    BOND *bp, *get_bond_pointer();
    float r,x1,y1,z1;
    float keff,length_eff;
    ATOM *a1,*a2,*a3;


    abcp = abc_first;
    if( abcp == NULL ){ return v_bond( V,lambda); }
    bp = bond_first;
    if( bp == NULL ) return 1;

    do_abc( lambda);  /* setup the increments */

    /* loop through the BONDS and update the potentials */
    while(1)
    {
        if(bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active ) {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            r = x1*x1+ y1*y1 + z1*z1;
            if( bp->ndiff > 0 )
            {
                keff = bp->k + bp->dk/bp->ndiff;
                length_eff = bp->length + bp->dlength/bp->ndiff;
            }else{
                keff = bp->k ;
                length_eff = bp->length ;
            }
            r = sqrt(r);
            r = r - length_eff;
            *V += keff*r*r;
        } /* if active */
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
/* do_abc()
*
*  do the abc corrections, this is common to both v and f_abc
* since only called from v_abc and f_abc the existance
* and completness of the abc and bond list has been checked
*
*/
do_abc(lambda)
float lambda;
{
    ABC *abcp;
    BOND *bp,*get_bond_pointer();
    float r,x1,y1,z1,x2,y2,z2,dp;
    ATOM *a1,*a2,*a3;



    bp = bond_first;
    /* initialize the bond list */
    while(1)
    {
        if( bp == NULL ) break;
        bp->ndiff = 0;
        bp->dlength = 0.;
        bp->dk = 0.;
        if( bp == bp->next) break;
        bp = bp->next;
    }
    /* then walk through the abc list and adjust the increments */
    abcp = abc_first;
    while(1)
    {
        if( abcp == NULL ) break;

        if( abcp->b1 == NULL )
        { abcp->b1 = get_bond_pointer(abcp->atom1,abcp->atom2);
            if( abcp->b1 == NULL ) goto SKIP; }
        if( abcp->b2 == NULL )
        { abcp->b2 = get_bond_pointer(abcp->atom3,abcp->atom2);
            if( abcp->b2 == NULL ) goto SKIP;}
        /* if here everything is defined for a given abc so lets do it */

        a1 = abcp->atom1; a2 = abcp->atom2; a3 = abcp->atom3;
        x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
        y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
        z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
        x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
        y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
        z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));
        dp = x1*x2+y1*y2+z1*z2;
        r = (( x1*x1+y1*y1+z1*z1)*(x2*x2+y2*y2+z2*z2));
        r = sqrt(r);
        dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
        dp = acos(dp);

        r = dp - abcp->theta;
        r = r*(  abcp->thetazero - r );
        bp = abcp->b1;
        bp->ndiff += 1;
        bp->dlength += abcp->dr*r;
        bp->dk += abcp->dk1*r;
        bp = abcp->b2;
        bp->ndiff += 1;
        bp->dlength += abcp->dr*r;
        bp->dk += abcp->dk2*r;

SKIP: ;
        if( abcp == abcp->next) break;
        abcp = abcp->next;
    }
}



/* f_abc()
*
* f_abc increments the forces in the atom structures by the force
* due to the abc components.  NOTE THE WORD increment.
* the forces should first be zero'd.
* if not then this code will be invalid.  THIS IS DELIBERATE.
* on bigger (and better?) machines the different potential terms
* may be updated at random or in parrellel, if we assume that this routine
* will initialize the forces then we can't do this.
*/
int f_abc(lambda)
float lambda;
/*  returns 0 if error, 1 if OK */
{
    ABC *abcp;
    BOND *bp;
    float r,k,ux1,uy1,uz1,ux2,uy2,uz2;
    float rb,ubx,uby,ubz;
    ATOM *a1,*a2,*a3;
    float x1,y1,z1,x2,y2,z2;
    float r1,r2,dtheta,dp,drb;
    float r11,r22,sdth;
    float length_eff;
    float ddtheta;

    abcp = abc_first;
    if( abcp == NULL ){ return f_bond( lambda); }
    bp = bond_first;
    if( bp == NULL ) return 1;

    do_abc( lambda);  /* setup the increments */

    /*first loop over the bonds !!
    * iff the bp->ndiff == 0 then do the bond force 
    * just as in bonds.c 
    * otherwise don't do anything and we will do the forces
    * from the abc loop
     */
    while(1)
    {
        if( bp == NULL) break;
        if( bp->ndiff == 0 ){
            k = bp->k;
            length_eff = bp->length;
        }else{
            k = bp->k + bp->dk/bp->ndiff;
            length_eff = bp->length+ bp->dlength/bp->ndiff;
        }
        a1 = bp->atom1; a2 = bp->atom2;
        if( a1->active || a2->active ) {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            dp = x1*x1+y1*y1+z1*z1;
            if( dp < 1.e-15 )
            { rb = 0; ubx = 1.; uby = 0.; ubz = 0.;
            } else{
                rb = sqrt(dp);
                ubx = x1/rb;
                uby = y1/rb;
                ubz = z1/rb;
            }
            ubx =  2*k*(rb - length_eff)*ubx;
            uby =  2*k*(rb - length_eff)*uby;
            ubz =  2*k*(rb - length_eff)*ubz;

        } /* end of active check */
        if( a1->active){
            a1->fx -= ubx;
            a1->fy -= uby;
            a1->fz -= ubz;
        }

        if( a2->active){
            a2->fx += ubx;
            a2->fy += uby;
            a2->fz += ubz;
        }

SKIP:
        if( bp == bp->next ) break ;
        bp = bp->next;
    }

    abcp = abc_first;
    while(1)
    {
        if( abcp == NULL ) break;
        a1 = abcp->atom1; a2 = abcp->atom2; a3 = abcp->atom3;

        if( a1->active|| a2->active || a3->active )
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));
            x2 = (a3->x -a2->x +lambda*(a3->dx-a2->dx));
            y2 = (a3->y -a2->y +lambda*(a3->dy-a2->dy));
            z2 = (a3->z -a2->z +lambda*(a3->dz-a2->dz));

            dp = x1*x2+y1*y2+z1*z2;
            r1 = sqrt(x1*x1+y1*y1+z1*z1);
            r2 = sqrt(x2*x2+y2*y2+z2*z2);
            if( r1 < 1.e-5 || r2 < 1.e-5) goto SKIP_ANGLE;
            r = r1*r2;
            if( r > 1.e-8){

                dp = dp/r;  if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;
                dtheta = acos(dp);
                ddtheta = dtheta - abcp->theta;
                ddtheta =  abcp->thetazero -ddtheta -ddtheta;
                sdth = sin(dtheta); if( sdth < 1.e-3) sdth = 1.e-3;
                r11 = r1*sdth; r22 = r2*sdth;
                dtheta = 0.;
                bp = abcp->b1;
                if( bp->ndiff > 0 ){
                    rb = r1 - bp->length - bp->dlength/bp->ndiff;
                    dtheta = abcp->dk1*ddtheta*rb*rb/bp->ndiff;
                    dtheta += 2*abcp->dr*ddtheta *rb *(bp->k + bp->dk/bp->ndiff)/bp->ndiff;
                }
                bp = abcp->b2;
                if( bp->ndiff > 0 ){
                    rb = r2 - bp->length - bp->dlength/bp->ndiff;
                    dtheta += abcp->dk2*ddtheta*rb*rb/bp->ndiff;
                    dtheta += 2*abcp->dr*ddtheta *rb *(bp->k + bp->dk/bp->ndiff)/bp->ndiff;
                }
                /* these lines calculate the x1,x2,x3 dependence
                * of d/dtheta and are the same for all derivatives */
                ux2 = -(x1/r1 - dp*x2/r2)/r22*dtheta;
                uy2 = -(y1/r1 - dp*y2/r2)/r22*dtheta;
                uz2 = -(z1/r1 - dp*z2/r2)/r22*dtheta;
                ux1 = -(x2/r2 - dp*x1/r1)/r11*dtheta;
                uy1 = -(y2/r2 - dp*y1/r1)/r11*dtheta;
                uz1 = -(z2/r2 - dp*z1/r1)/r11*dtheta;
                if( a1->active)
                {
                    a1->fx += ux1;
                    a1->fy += uy1;
                    a1->fz += uz1;
                }

                if( a2->active)
                {
                    a2->fx += -ux1 - ux2;
                    a2->fy += -uy1 - uy2;
                    a2->fz += -uz1 - uz2;
                }

                if( a3->active)
                {
                    a3->fx += ux2;
                    a3->fy += uy2;
                    a3->fz += uz2;
                }
            }
        }

SKIP_ANGLE:
        if( abcp == abcp->next) break;
        abcp = abcp->next;
    }




}
/* function get_abc( a1,bonded,10,inbond);
* check the ABC list for atoms 1-3 ed to a1
*/
void get_abc( a1,bonded,mbond,inbond)
ATOM *a1, *bonded[];
int mbond,*inbond ;
{
    ABC *mine;
    mine = abc_first;
    *inbond = 0;
    while(1)
    {
        if( (mine == NULL) )
        {
            return;
        }
        if( mine->atom1 == a1)
        {
            bonded[(*inbond)++] = mine->atom3;
        }
        if( mine->atom3 == a1)
        {
            bonded[(*inbond)++] = mine->atom1;
        }
        if( mine == mine->next) return;
        mine = mine->next;
        if( *inbond == mbond ) return;
    }
}
/* routine dump_abcs
* this function outputs the abc parameters
* and does it in a simple form
* abc ser1,ser2,ser3,k,r,theta (in degrees )
* the rest is just free format
*/
void dump_abcs( where )
FILE *where;
{
    ABC *b;
    ATOM *a1,*a2,*a3;
    float rtodeg;
    b = abc_first;
    if( b == NULL ) return;
    rtodeg = 180./acos(-1.);
    while( (b->next != b) )
    {
        if( b->next == NULL) return;
        a1 = b->atom1; a2 = b->atom2;a3 = b->atom3;
        fprintf( where,"abc %d %d %d %f %f %f %f %f \;\n",a1->serial,a2->serial,
                 a3-> serial,b->theta*rtodeg,b->thetazero*rtodeg,
                 b->dr, b->dk1, b->dk2);
        b = b->next;
    }
    if( b->next == NULL) return;
    a1 = b->atom1; a2 = b->atom2;a3 = b->atom3;
    fprintf( where,"abc %d %d %d %f %f %f %f %f\;\n",a1->serial,a2->serial,
             a3-> serial,b->theta*rtodeg,b->thetazero*rtodeg,
             b->dr, b->dk1, b->dk2);
}

/* a_abc()
* this function sums up the potentials
* for the atoms defined in the abc data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int a_abc( V, lambda,ilow,ihigh,op )
float *V,lambda;
int ilow,ihigh;
FILE *op;
{
    BOND *bp;
    float r,x1,y1,z1;
    float keff,length_eff;
    ATOM *a1,*a2;


    bp = bond_first;
    if( bp == NULL ) return 1;
    if( abc_first == NULL ) return 1;
    /* modify the bonds as needed  */
    do_abc(lambda);
    while(1)
    {
        if( bp == NULL) return 0;
        a1 = bp->atom1; a2 = bp->atom2;
        if( (a1->serial >= ilow && a1->serial <= ihigh)
                ||  (a2->serial >= ilow && a2->serial <= ihigh))
        {
            x1 = (a1->x -a2->x +lambda*(a1->dx-a2->dx));
            y1 = (a1->y -a2->y +lambda*(a1->dy-a2->dy));
            z1 = (a1->z -a2->z +lambda*(a1->dz-a2->dz));

            r = x1*x1 + y1*y1 + z1*z1;
            r = sqrt(r);
            if( bp->ndiff == 0 )
            {
                keff = bp->k;
                length_eff = bp->length;
            } else {
                keff = bp->k + bp->dk/bp->ndiff;
                length_eff = bp->length + bp->dlength/bp->ndiff;
            }
            *V += keff*(r-length_eff)*(r-length_eff);

            fprintf(op,"abc bond %d %d  E %f r %f k %f error %f\n"
                    ,a1->serial,a2->serial, keff*(r-length_eff)*(r-length_eff),
                    r,keff,r-length_eff);
        }
        if( bp == bp->next ) return 1;
        bp = bp->next;
    }
}
