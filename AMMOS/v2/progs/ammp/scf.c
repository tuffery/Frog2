/* v_scf,f_fscf routines
*  
*  after orbitals have been refined track the solutions in 
*  3-d space using the dscf suite
*
*
*/
/*
*  copyright 1996 Robert W. Harrison
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
#include "orbit.h"

#define DSTEP 1.e-2
#define STEP 1.e2

/*
int v_scf( float *V, float lambda)
*/
int v_scf( V,lambda)
float *V, lambda;
{
    int i;
    int no,o_number();
    int a_inc_d();
    int renormalize();
    int o_update_normals();
    float Just_Ritz_product();
    float x;

    no = o_number();
    if( no == 0 ) return 0;

    /* this is a terrible solution, but one we need to
    * use, the alternative is x+dx in all of the phixphi routines
    *  which is even worse
    */
    if( lambda != 0.) a_inc_d(lambda);
    o_update_normals();
    renormalize();
    x = Just_Ritz_product(NULL);
    if( lambda != 0) a_inc_d(-lambda);

    *V += x*HARTREE;
}

/*
int f_scf( float lambda)
*/
int f_scf( lambda)
float lambda;
{
    int i;
    int no,o_number();
    int na,a_number();
    ATOM *ap,*a_next();
    int a_inc_d();
    float INDO_Ritz_product();
    float baseV,x;
    float baseN,orth_justone_nuclear(),orth_nuclear();
    ORBITAL *o1,*o2,*o_next();
    int i1,i2;
    void orth_invert();
    int renormalize();
    int o_update_normals();
    float (*t)[],(*ov)[],(*ok)[];

    return 0;

    no = o_number();
    if( no == 0 ) return 0;
    t = malloc(no*no*sizeof(float));
    ov = malloc(no*no*sizeof(float));
    ok = malloc(no*no*sizeof(float));
    if( t == NULL || ov == NULL || ok == NULL){
        aaerror("cannot allocate memory in f_scf"); exit(0);}
    /* setup the overlap inverse */
    orth_invert(t,ov,ok,no);

    na = a_number();

    /* this is a terrible solution, but one we need to
    * use, the alternative is x+dx in all of the phixphi routines
    *  which is even worse
    */
    if( lambda != 0) a_inc_d(lambda);
    o_update_normals();
    renormalize();
    baseV = INDO_Ritz_product(NULL);
    ap = a_next(0);
    for( i=0; i< na; i++,ap=ap->next)
    {
        if( ap->active){
            if( ap->na > 0 ){
                ap->x += DSTEP;
                o_update_normals();
                x = INDO_Ritz_product(NULL);
                ap->fx -= STEP*(x-baseV)*HARTREE;
                ap->x -= DSTEP;
                ap->y += DSTEP;
                o_update_normals();
                x = INDO_Ritz_product(NULL);
                ap->fy -= STEP*(x-baseV)*HARTREE;
                ap->y -= DSTEP;
                ap->z += DSTEP;
                o_update_normals();
                x = INDO_Ritz_product(NULL);
                ap->fz -= STEP*(x-baseV)*HARTREE;
                ap->z -= DSTEP;
            } else {
                /* have to do the others here */
                baseN = 0.;
                o1 = o_next(-1);
                for( i1=0; i1<  no; i1++)
                {
                    for( i2=0; i2<  no; i2++)
                    {
                        baseN += orth_justone_nuclear(ap,o1,o2,1.)*(*t)[i1*no+i2];
                        o2 = o2->next;
                    } o1 = o1->next;
                }
                ap->x += DSTEP;
                x = 0.;
                for( i1=0; i1<  no; i1++)
                {
                    for( i2=0; i2<  no; i2++)
                    {
                        x += orth_justone_nuclear(ap,o1,o2,1.)*(*t)[i1*no+i2];
                        o2 = o2->next;
                    } o1 = o1->next;
                }
                ap->x -= DSTEP;
                ap->fx += STEP*(baseN-x)*HARTREE;
                ap->y += DSTEP;
                x = 0.;
                for( i1=0; i1<  no; i1++)
                {
                    for( i2=0; i2<  no; i2++)
                    {
                        x += orth_justone_nuclear(ap,o1,o2,1.)*(*t)[i1*no+i2];
                        o2 = o2->next;
                    } o1 = o1->next;
                }
                ap->y -= DSTEP;
                ap->fy += STEP*(baseN-x)*HARTREE;
                ap->z += DSTEP;
                x = 0.;
                for( i1=0; i1<  no; i1++)
                {
                    for( i2=0; i2<  no; i2++)
                    {
                        x += orth_justone_nuclear(ap,o1,o2,1.);
                        o2 = o2->next;
                    } o1 = o1->next;
                }
                ap->z -= DSTEP;
                ap->fz += STEP*(baseN-x)*HARTREE;
            } /* if na > 0 */
        }/* active */
    }/* for i */
    if( lambda != 0) a_inc_d(-lambda);
    free( ok); free(ov); free(t);
}
