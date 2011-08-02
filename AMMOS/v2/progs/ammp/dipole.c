/* dipole.c
*
* calculate the dipole moment for a given range of atoms
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

void dipole( op,first,last )
FILE *op;
int first,last;
{
    ATOM *ap,*a_m_serial();
    float x,y,z,sumq,magnitude;
    int i,sum ;

    if( last < first )
    {
        i = last; last = first; first = i;
    }

    x = 0.; y = 0.; z = 0.;
    sumq = 0.;
    sum = 0;

    for( i=first; i<= last; i++)
    {
        ap = a_m_serial(i);
        if( ap != NULL )
        {
            x += ap->x * ap->q;
            y += ap->y * ap->q;
            z += ap->z * ap->q;
            sumq += ap->q;
            sum += 1;
        }
    }
    if(sum == 0)
    { fprintf(op,"Dipole warning, no atoms in %d %d\n",first,last);
        return;
    }
    /* check that a dipole exists */
    if( abs(sumq) > 1.e-5)
    {
        fprintf(op,"Dipole warning,%d %d sum charge is not zero %f\n"
                ,first,last,sumq);
        x = x/sum; y = y/sum; z = z/sum;
        fprintf(op,"Dipole %d %d center of charge is %f %f %f\n"
                ,first,last,x,y,z);
        return;
    }
    /* finally we have a dipole that is real */
    sumq = sqrt(x*x + y*y + z*z);
    fprintf(op,"Dipole %d %d moment %f vector %f %f %f\n",
            first,last,sumq,x,y,z);
    return;
}
