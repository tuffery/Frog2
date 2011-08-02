/* tgroup.c
*
* collection of routines to service  torsion searching
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
*/
/*
*  copyright 1993 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-torsioned potentials
*/
#define MAXTGROUP 100
#define MAXTGSTEP 24
typedef struct { ATOM *context,*b1,*b2,*b3 ;
    int which ,ntry, ingroup ;
    float base;
    ATOM *group[MAXTGROUP];
    void *next;
}  TGROUP ;
TGROUP *tg_first = NULL;
#define tg_long sizeof(TGROUP) 

int tgroup( which,context,b1,b2,b3,base,ntry )
int which ,context,b1,b2,b3,ntry ; /* context b1 b2 b3 are atom serial numbers */
float base;
{

    TGROUP *tgp,*tgo;
    ATOM *bonded[20],*a_m_serial();
    ATOM *newest[MAXTGROUP];
    ATOM *newer[MAXTGROUP];
    int in_newest,in_newer;
    int  i,j,k,l,ll;
    char line[80];
    void get_bond();
    if( ntry == 0 )
    { base = 0. ; ntry = 6 ;}

    tgp = tg_first;
    tgo = tgp;
    if( which == 0)
    { aaerror("warning tg_group 0 will never be accessed\n"); return 0;}
    while(tgp != NULL)
    {
        if( tgp->which == which) goto found;
        tgo = tgp;
        tgp = tgp->next;
    }
    if( (tgp = malloc( tg_long)) == NULL)
    { aaerror(" cannot allocate memory for a tg_group\n"); exit(0);}
    if( tgo != NULL) tgo->next = tgp;
    if( tgo == NULL) tg_first = tgp;
    tgp->next = NULL;
found:
    if( context < 0 || b1 < 0 || b2 < 0 || b3 < 0) /* free the group iff it is not definable */
    {
        if( tgo != NULL ) tgo->next = tgp->next; /* delete the link if its not the first one */
        if( tgo == NULL ) tg_first = tgp->next;  /* delete the first link */
        free( tgp );
        return 1;
    } /* else reinitialize the group */
    tgp->which = which;
    tgp->context = a_m_serial(context);
    tgp->b1 = a_m_serial(b1);
    tgp->b2 = a_m_serial(b2);
    tgp->b3 = a_m_serial(b3);
    if( tgp->context == NULL )
    {
        sprintf(line," tgroup %d  %d atom not defined cannot define tgroup ",which,context);
        aaerror( line);
        return 1;}
    if( tgp->b1 == NULL )
    {
        sprintf(line," tgroup %d  %d atom not defined cannot define tgroup ",which,b1);
        aaerror( line);
        return 1;}
    if( tgp->b2 == NULL )
    {
        sprintf(line," tgroup %d  %d atom not defined cannot define tgroup ",which,b2);
        aaerror( line);
        return 1;}
    if( tgp->b3 == NULL )
    {
        sprintf(line," tgroup %d  %d atom not defined cannot define tgroup ",which,b3);
        aaerror( line);
        return 1;}
    tgp->base = base *3.1415926589793 /180.;
    for( k=0; k< MAXTGROUP; k++)
        tgp->group[k] = NULL;
    if( ntry > MAXTGSTEP) ntry = MAXTGSTEP;
    tgp->ntry = ntry;
    get_bond( tgp->b2, bonded,20,&j);
    for( k=0; k< j; k++)
    {
        if( bonded[k] != tgp->context
                && bonded[k] != tgp->b3
                && bonded[k] != tgp->b1 ) break;
    }
    newest[0] = tgp->b3;
    tgp->group[0] = tgp->b3;
    tgp->ingroup = 1;
    in_newest = 1;
    if( k != j ){
        newest[1] = bonded[k];
        tgp->group[1] = bonded[k];
        tgp->ingroup = 2;
        in_newest = 2;
        for( i = 0; i < j; i++)
        {
            l = 1;
            for( k=0; k< in_newest; k++)
            { if(newest[k] == bonded[i]) { l = 0; break;} }
            if( bonded[i] == tgp->context) l = 0;
            if( bonded[i] == tgp->b1) l = 0;
            if( bonded[i] == tgp->b2) l = 0;
            if( bonded[i] == tgp->b3) l = 0;

            if( l == 1) {
                tgp->group[tgp->ingroup++ ] = bonded[i];
                newest[in_newest++] = bonded[i];
            }
        }
    }
    while( in_newest > 0 )
    {
        in_newer = 0;
        for( l=0; l < in_newest; l++)
        {
            get_bond( newest[l],bonded,20,&j);
            for( i =0; i< j; i++)
            {
                ll = 1;
                for( k=0; k< tgp->ingroup ; k++)
                    { if( tgp->group[k] == bonded[i])
                    { ll = 0; break; }
                }
                if( bonded[i] == tgp->context) ll = 0;
                if( bonded[i] == tgp->b1) ll = 0;
                if( bonded[i] == tgp->b2) ll = 0;
                if( bonded[i] == tgp->b3) ll = 0;
                if( ll == 1 )
                {
                    tgp->group[ tgp->ingroup++ ] = bonded[i];
                    newer[in_newer++] = bonded[i];
                    if( tgp->ingroup > MAXTGROUP )
                    {aaerror(" too many atoms in a tgroup - must exit\n"); exit(0);}
                }
            } /* end of loop i */
        }/* end of loop l */
        for( i= 0; i< in_newer; i++)
        { newest[i] = newer[i]; }

        in_newest = in_newer;
    }
    if( tgp->ingroup > MAXTGROUP )
    {aaerror(" too many atoms in a tgroup - must exit\n"); exit(0);}
    /*
    	for( i = 0; i < tgp->ingroup; i++)
    	printf(" atom %d in group %d\n",tgp->group[ i ]->serial,which); */
    return 1;
}/* end of routine tgroup */

/* tsearch()
*  this function performs the search over a bunch of torsion groups
*
*  tsearch ( n1 n2 n3 ...n8 )
* 
*  recurses n8 ...n3 n2 n1  from the last (deepest) non-zero tgroup 
*
*
* makes use of various registers 
*  best is returned in xyz registers
*/
tsearch( t1 ,t2 ,t3,t4,t5,t6,t7,t8 )
int t1,t2,t3,t4,t5,t6,t7,t8;
{

    TGROUP  *grouplist[8],*match_tgroup();
    int ngroup,i,j;
    int bestlist[8];
    float V;
    void tg_do_search();
    void tg_gen_con(), tg_init(),tg_apply();

    ngroup = 8;
    for( i=0; i< ngroup ; i++)
        grouplist[i] = NULL;

    if( t8 <= 0 ) ngroup = 7;
    if( t7 <= 0 ) ngroup = 6;
    if( t6 <= 0 ) ngroup = 5;
    if( t5 <= 0 ) ngroup = 4;
    if( t4 <= 0 ) ngroup = 3;
    if( t3 <= 0 ) ngroup = 2;
    if( t2 <= 0 ) ngroup = 1;
    if( t1 <= 0 ) ngroup = 0;
    if( ngroup == 0 ) return 1;
    grouplist[7] = match_tgroup( t8 );
    grouplist[6] = match_tgroup( t7 );
    grouplist[5] = match_tgroup( t6 );
    grouplist[4] = match_tgroup( t5 );
    grouplist[3] = match_tgroup( t4 );
    grouplist[2] = match_tgroup( t3 );
    grouplist[1] = match_tgroup( t2 );
    grouplist[0] = match_tgroup( t1 );
    /* check active atoms
    * if an  movable atom is inactive make the tgroup NULL 
    * movable atoms include >bp4 and everyone in the list 
    */
    for( i=0; i< ngroup; i++)
    {
        if( grouplist[i] != NULL)
            { for( j=0; j< grouplist[i]->ingroup; j++)
                if( !(grouplist[i]->group[j]->active) )
                {grouplist[i] = NULL; break;}
            if( grouplist[i] != NULL )
                if(!(grouplist[i]->b3->active) )
                {grouplist[i] = NULL; }
        }
    }
    /* now check for NULL tgroups */
    for( i = 0; i< ngroup; i++)
    {
        if( grouplist[i] == NULL)
        {
            for( j=i; j< ngroup-1; j++)
            { grouplist[j] = grouplist[j+1]; }
            ngroup -= 1;
        }
    }


    /* set all the torsions to the ideal values in x,y,z */
    tg_init( grouplist,&ngroup,0 );
    /* do the search  (recursive) */
    tg_do_search(&V, grouplist,bestlist,0,ngroup );
    /*
    	for( i=0; i< ngroup; i++)
    	printf(" group %d best %d\n",i,bestlist[i]);
    */
    /* fixup the results */
    /*	tg_init( grouplist,&ngroup ); */
    tg_gen_con( grouplist,bestlist,ngroup );

} /* end of tsearch */

/*	tg_gen_con( grouplist,bestlist,ngroup );
*
*  given the list of the best offsets
*  generate the best coordinates and
* return them to  xyz from dxdydz
*/
void tg_gen_con( gl,bl,ngl)
TGROUP *gl[];
int bl[], ngl;
{
    ATOM *a1;
    int i,j;
    void tg_apply();
    for( i=0 ; i < ngl ; i++)
    {
        j = ngl -i -1;
        tg_apply( gl[j],bl[j]) ;
    }
    for( i=ngl-1; i> -1 ; i -- )
    {
        for( j=0; j< gl[i]->ingroup ; j++)
        {
            a1 = gl[i]->group[j];
            a1->x = a1->dx;
            a1->y = a1->dy;
            a1->z = a1->dz;
        }
    }
}/* end of tg_gen_con */

/* tg_init( gl,ngl )
* go through all of the tgroups specified
*  find the torsion value and then rotate so
* that the torsion is  at the base (gl[i]->base)
* value. this prepares the tgroup set for
* a recursive search 
*/
void tg_init( gl,ng ,deep)
TGROUP *gl[];
int *ng,deep;
{
    ATOM *a1,*a2,*a3,*a4;
    float x1,y1,z1,x2,y2,z2,x3,y3,z3;
    float cx1,cy1,cz1,cx2,cy2,cz2;
    float dp,r;
    int i ,j;
    int ngl;
    void tg_d_apply();
    ATOM *a_next();
    int  a_number();
    if( deep == 0 )
    {
        j = a_number();
        for( i=0; i< j; i++)
        {
            a1 = a_next(i);
            a1->dx = a1->x;
            a1->dy = a1->y;
            a1->dz = a1->z;
        }
    }
    ngl = *ng;
    for( i=ngl-1; i> -1 ; i -- )
    {
        if( gl[i] == NULL) {if( i == 0) {*ng = 0; return;}
            *ng = i; tg_init( gl,ng,deep+1); return ;}
        for( j=0; j< gl[i]->ingroup ; j++)
        {
            a1 = gl[i]->group[j];
            if( a1 == NULL) return;
            a1->dx = a1->x;
            a1->dy = a1->y;
            a1->dz = a1->z;
        }
        a1 = gl[i]->context;
        a2 = gl[i]->b1;
        a3 = gl[i]->b2;
        a4 = gl[i]->b3;

        a1->dx = a1->x;
        a1->dy = a1->y;
        a1->dz = a1->z;
        a2->dx = a2->x;
        a2->dy = a2->y;
        a2->dz = a2->z;
        a3->dx = a3->x;
        a3->dy = a3->y;
        a3->dz = a3->z;
        a4->dx = a4->x;
        a4->dy = a4->y;
        a4->dz = a4->z;

        x1 = (a1->x -a2->x );
        y1 = (a1->y -a2->y );
        z1 = (a1->z -a2->z );
        x2 = (a3->x -a2->x );
        y2 = (a3->y -a2->y );
        z2 = (a3->z -a2->z );
        x3 = (a4->x -a3->x );
        y3 = (a4->y -a3->y );
        z3 = (a4->z -a3->z );
        /* 1 cross 2 */
        cx1 = y1*z2 - y2*z1;
        cy1 = -x1*z2 + x2*z1;
        cz1 = x1*y2 - x2*y1;
        r = cx1*cx1 + cy1*cy1 + cz1*cz1;
        if( r < 1.e-4) goto SKIP;
        r = sqrt(r);
        cx1 = cx1/r;
        cy1 = cy1/r;
        cz1 = cz1/r;
        /* 3 cross 2 */
        cx2 = y3*z2 - y2*z3;
        cy2 = -x3*z2 + x2*z3;
        cz2 = x3*y2 - x2*y3;
        r = cx2*cx2 + cy2*cy2 + cz2*cz2;
        if( r < 1.e-4) goto SKIP;
        r = sqrt(r);
        cx2 = cx2/r;
        cy2 = cy2/r;
        cz2 = cz2/r;
        /* if here everything is well determined */
        dp = cx1*cx2 + cy1*cy2 + cz1*cz2; /* cos( abs(theta)) */
        if( dp > 1.) dp = 1.; if( dp < -1.) dp = -1.;

        dp = acos(dp);
        /* determine the sign by triple product */
        r = cx1*x3 + cy1*y3 + cz1*z3;
        if( r > 0 ) dp =  -dp ;
        r =  gl[i]->base - dp ;
        tg_d_apply( gl[i], r );
        for( j=0; j< gl[i]->ingroup ; j++)
        {
            a1 = gl[i]->group[j];
            a1->x = a1->dx;
            a1->y = a1->dy;
            a1->z = a1->dz;
        }
SKIP:   i = i;
    } /* end of for( i */


} /* end of tg_init */

/* tg_do_search,
*
* recursively, (depth first) evaluate the tgroup tree in
* grouplist 
*/
void tg_do_search( Vp, gl,bl,igl, ngl )
float *Vp;
int igl,ngl,bl[];
TGROUP *gl[];
{
    float vl[MAXTGSTEP],vb;
    int i,ibest;
    int bestlist[MAXTGSTEP][8];
    int lbl[8];
    void tg_apply();
    /* safety first - these  should never happen */
    if( igl == 8 ) return;
    if( gl[igl] == NULL ) return;
    if( igl > ngl) return;

    for( i=0; i < gl[igl]->ntry ; i++)
    {
        tg_do_search(&vl[i], gl,lbl,igl+1,ngl );
        /* evaluate my energy */
        /* int tg_nonbon( V,tgp ) */
        tg_nonbon( &vl[i], gl[igl]);
        for( ibest=0; ibest<8; ibest++)
            bestlist[i][ibest] = lbl[ibest];
        /* generate the next step */
        if( gl[igl]->ntry > 1)
        {
            tg_apply( gl[igl],1);
        }/* end of if worth generating next try */
    } /* end of for each step (i) loop */
    vb = 10e20;
    ibest = 0;
    for( i=0; i< gl[igl]->ntry; i++)
    {
        /*	printf(" tg_do_search: %d %f\n",i,vl[i]); */
        if( vl[i] < vb ) { vb = vl[i]; ibest = i ;}
    }
    bestlist[ibest][igl] = ibest;
    for( i=0; i<8; i++)
    { bl[i] = bestlist[ibest][i]; }
    *Vp = vl[ibest];

}/* end of tg_do_search */

/* tg_apply
*  apply the torsion to the torsion group
* tg_apply( TGROUP *tgp, int num )
*  pointer and how many steps to use
*/
void tg_apply( tgp, num )
TGROUP *tgp;
int num;
{

    float nx,ny,nz;
    float phi,cphi,sphi;
    float rx,ry,rz, nnrx,nnry,nnrz, rnx,rny,rnz;
    ATOM *b1,*b2;
    int i;

    b1 = tgp->b1;
    b2 = tgp->b2;
    nx = b2->dx - b1->dx;
    ny = b2->dy - b1->dy;
    nz = b2->dz - b1->dz;
    rx = sqrt(nx*nx + ny*ny + nz*nz);
    if( rx < 1.e-6)
    {aaerror(" bad torsion radius in tg_apply \n"); return ;}
    nx = nx/rx;
    ny = ny/rx;
    nz = nz/rx;
    phi = 2.*3.141592653589793 /(float)tgp->ntry * (float)num;
    cphi = cos(phi); sphi = sin(phi);
    for( i=0; i< tgp->ingroup; i++)
    {
        rx = (tgp->group[i])->dx - b1->dx;
        ry = (tgp->group[i])->dy - b1->dy;
        rz = (tgp->group[i])->dz - b1->dz;
        phi = nx*rx + ny*ry + nz*rz;
        nnrx = phi*nx;
        nnry = phi*ny;
        nnrz = phi*nz;
        rnx = ny*rz - nz*ry;
        rny = -nx*rz + nz*rx;
        rnz = nx*ry - ny*rx;
        phi = (1.-cphi);
        rx = cphi*rx + phi*nnrx + sphi*rnx;
        ry = cphi*ry + phi*nnry + sphi*rny;
        rz = cphi*rz + phi*nnrz + sphi*rnz;
        (tgp->group[i])->dx = rx + b1->dx;
        (tgp->group[i])->dy = ry + b1->dy;
        (tgp->group[i])->dz = rz + b1->dz;
    }
}/* end of tg_apply */

/* tg_d_apply
*  apply the torsion to the torsion group
* tg_d_apply( TGROUP *tgp, float off )
* a pointer and the amount (radians) to rotate 
*/
void tg_d_apply( tgp, off )
TGROUP *tgp;
float off;
{

    float nx,ny,nz;
    float phi,cphi,sphi;
    float rx,ry,rz, nnrx,nnry,nnrz, rnx,rny,rnz;
    ATOM *b1,*b2;
    int i;

    b1 = tgp->b1;
    b2 = tgp->b2;
    nx = b2->dx - b1->dx;
    ny = b2->dy - b1->dy;
    nz = b2->dz - b1->dz;
    rx = sqrt(nx*nx + ny*ny + nz*nz);
    if( rx < 1.e-6)
    {aaerror(" bad torsion radius in tg_apply \n"); return ;}
    nx = nx/rx;
    ny = ny/rx;
    nz = nz/rx;
    /*	phi = 2.*3.141592653589793 /tgp->ntry * num;
    */
    phi = off;
    cphi = cos(phi); sphi = sin(phi);
    for( i=0; i< tgp->ingroup; i++)
    {
        rx = (tgp->group[i])->dx - b1->dx;
        ry = (tgp->group[i])->dy - b1->dy;
        rz = (tgp->group[i])->dz - b1->dz;
        phi = nx*rx + ny*ry + nz*rz;
        nnrx = phi*nx;
        nnry = phi*ny;
        nnrz = phi*nz;
        rnx = ny*rz - nz*ry;
        rny = -nx*rz + nz*rx;
        rnz = nx*ry - ny*rx;
        phi = (1.-cphi);
        rx = cphi*rx + phi*nnrx + sphi*rnx;
        ry = cphi*ry + phi*nnry + sphi*rny;
        rz = cphi*rz + phi*nnrz + sphi*rnz;
        (tgp->group[i])->dx = rx + b1->dx;
        (tgp->group[i])->dy = ry + b1->dy;
        (tgp->group[i])->dz = rz + b1->dz;
    }
}/* end of tg_d_apply */

/* dump_tgroup
*  outputs the tgroup stuff to the supplied output FILE pointer
*/
void dump_tgroup( where )
FILE *where;
{
    TGROUP *tgp;

    tgp = tg_first;
    while (tgp != NULL)
    {
        fprintf( where," tgroup %d %d %d %d %d %f %d ;\n",
                 tgp->which, tgp->context->serial,
                 tgp->b1->serial,
                 tgp->b2->serial,
                 tgp->b3->serial,
                 tgp->base*180./3.141592653589793, tgp->ntry);
        tgp = tgp->next;
    }
}

/* match_tgroup finds the tgroup which matches a group id
*  returns its address  or NULL 
*/
TGROUP *match_tgroup( i )
int i;
{
    TGROUP *tgp ;
    tgp = tg_first;
    while( tgp != NULL )
    {
        if( tgp ->which == i )
        {
            if( tgp->context == NULL) return NULL;
            if( tgp->b1 == NULL) return NULL;
            if( tgp->b2 == NULL) return NULL;
            if( tgp->b3 == NULL) return NULL;
            return tgp;
        }
        tgp = tgp->next;
    }
    return NULL;
}



/* tg_nonbon()
* this function sums up the potentials
* for the atoms defined in the nonbon data structure.
*/
/* standard returns 0 if error (any) 1 if ok
* V is the potential */
int tg_nonbon( V ,tgp )
float *V;
TGROUP *tgp ;
{
    float r,r0,xt,yt,zt;
    float lcutoff,cutoff,get_f_variable();
    int inbond,inangle,i,ii;
    ATOM *a1,*a2;
    ATOM *a_next( ); /* returns first ATOM when called with -1 */
    float dielectric,ve,va,vh;
    ATOM *a_m_serial();

    /* nonbonded potentials
    * do a double loop starting from the first atom to the 
    * last 
    * then from the second to the last 
    * etc
    *
    * also check to avoid bonded and 1-3 bonded atoms
    */
    dielectric = get_f_variable("dielec");
    if( dielectric < 1.) dielectric = 1.;
    dielectric = 332.17752/dielectric;
    cutoff = get_f_variable("cutoff");
    if( cutoff < 1.) cutoff = 1.e10;
    lcutoff = -cutoff;
    *V = 0.;
    if( tgp == NULL ) return 1;
    for( ii=0; ii<= tgp->ingroup; ii++)
    {
        a1 = tgp->group[ii];
        if( a1 == NULL ) goto NOTANATOM;
        ve = 0.; va = 0.; vh = 0.;
        a2 = a_next(-1);
        /*
        *	for(i = 0; i< a1->dontuse; i++)
        *	printf("%d ",a1->excluded[i]->serial);
        *	printf("\n");
        */
        /*
        	while(  (a2->next != a2) && (a2->next != NULL))
        	*/
        while(  (a2 != NULL) && (a2->next != NULL) && a2->next != a2)
        {
            /* goto SKIP is used because this is one case where it makes sense */
            /*	if( a2 == a1) break;  */
            if( a2 == a1) goto SKIP;
            for(i = 0; i< a1->dontuse; i++)
                if( a2 == a1->excluded[i]) goto SKIP;
            /* non - bonded are only used when the atoms arent bonded */

            xt = (a1->dx - a2->dx);
            if( (xt > cutoff) || (xt < lcutoff) ) goto SKIP;
            yt =  (a1->dy - a2->dy);
            if( (yt > cutoff) || (yt < lcutoff) ) goto SKIP;
            zt =  (a1->dz - a2->dz);
            if( (zt > cutoff) || (zt < lcutoff) ) goto SKIP;

            r = xt*xt+yt*yt+zt*zt;
            if( r < 1.) r = 1.;

            r0 = sqrt(r); r = r*r*r ;
            /* debugging
            *	printf(" %d %d %f %f %f \n", a1->serial,a2->serial,a1->q,a2->q,
            *	332.17752*a1->q*a2->q/r0);
            */
            ve += dielectric*a1->q*a2->q/r0;
            va -= a1->a*a2->a/r;
            vh += a1->b*a2->b/r/r;

SKIP:
            /*	if( a2->next == a1) break; */
            if( a2->next == a2) break;
            a2 = a2->next;
        }
        *V += ve + va + vh;
NOTANATOM:
        i = i;
    }
    return 1;

}

