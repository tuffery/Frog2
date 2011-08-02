/* this is a quicky demonstrating the power of using a library
*
* we need to calculate radial distributions of water.
*
* cannot do this real well by dumping coordinates as we will 
* need O(10k) frames.  
*
*  this is a bitch.
*
* instead we call ammp (via eval) and sum the data as needed
*
* since this is a C program we can directly extract the data from the
* ATOM structures
*/

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "ammp.h"

#define NHIST 1001
#define MAXAT 1000
#define NSTEP 250
#define NSTEP 1000

void main()
{

    ATOM *mylist[MAXAT];
    ATOM *a_next();
    float histogram[NHIST];
    float rmax,dr,rms;
    float dx,dy,dz;

    int i, j,k,numatom,iter;
    FILE *ip,*op;


    ip = stdin;
    op = stdout;

    /* initialize */
    set_f_variable("mxdq", .05);
    set_f_variable("mxcut", 6.);
    set_i_variable("nostep",1 );
    eval( ip,op,"read radial.input");

    rmax = 10.; dr = (NHIST-1)/rmax;
    for( i=0; i< NHIST; i++)
        histogram[i] = 0.;

    for( i=0; i< MAXAT; i++)
    {
        mylist[i] = a_next(i);
        while( mylist[i] != NULL && strcmp( &mylist[i]->name ,"h2o.o") != 0 )
        {
            mylist[i] = a_next(1);
        }
        if( mylist[i] == NULL) { numatom = i; break;}
    }

    if( numatom < 2 )
    {
        aaerror( " 2 or fewer atoms found"); exit(0);
    }
    /* run */

    for( iter = 0; iter< NSTEP; iter++)
    {
        eval(ip,op,"tpac 10 .00001 300 ");

        for( i=0 ; i< numatom-1; i++)
            for( j = i+1 ; j< numatom ; j++)
            {
                dx = mylist[i]->x - mylist[j]->x;
                dy = mylist[i]->y - mylist[j]->y;
                dz = mylist[i]->z - mylist[j]->z;
                rms = dx*dx + dy*dy + dz*dz;

                k = (int) (sqrt(rms)*dr + .5);
                if( k < NHIST) histogram[k] += 1;
            }
    }
    /* output  */
    for( i=0; i< NHIST; i++)
    {
        histogram[i] /= NSTEP;
        printf(" %d %f\n", i, histogram[i]);

    }
    for( i=1; i< NHIST; i++)
    {
        rms = i/dr*4*3.14159265;
        histogram[i] /=rms;
        printf(" %d %f\n", i, histogram[i]);
    }


}/* end of main()*/

