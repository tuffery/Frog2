/*   Damp the maximum force
*
* force the maximum force to have a prescribed value
*  look at the value kdamp (setf kdamp xxx)
*
*  if the maximum force is less than the value it is left alone
*
*  if no value is defined then use   
*
* v_damp returns nothing
* 
*
*/
/*
*  copyright 2000 Robert W. Harrison
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

/* v_damp is a NULL function */
int v_damp( V,lambda)
float *V,lambda;
{
    return 0;
}


int f_damp( lambda)
float lambda;
{
    /* lambda is ignored and the dx == vx in pac et al are
    * also updated */
    float kdamp,get_f_variable();
    float fmax,a_max_f(void);
    ATOM *ap,*a_next(int);
    int natom,a_number(void);
    int i;

    /* check for silly calls */
    natom = a_number();
    if( natom < 1) return 0;
    kdamp = get_f_variable("kdamp");
    if( kdamp < 1.) kdamp = 10.;
    fmax = sqrt(a_max_f());
    if( fmax < kdamp) return 0;

    fmax = kdamp/fmax;
    for( i=0; i< natom; i++)
    {
        ap = a_next(i);
        ap->fx *= fmax;
        ap->fy *= fmax;
        ap->fz *= fmax;
    }

    return 0;
}/* endof f_damp */
