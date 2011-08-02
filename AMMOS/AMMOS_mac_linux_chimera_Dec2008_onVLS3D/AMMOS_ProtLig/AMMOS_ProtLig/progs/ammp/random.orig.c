/* random number generators for ammp
*
*  randf() returns float in the range from 0 - 1.
*  randg() returns float gaussian distr sigma = 1
*  rand3( float *x1,*x2,*x3)  returns random unit 3vector
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

#include <math.h>
#include "ammp.h"
float randf()
{
    static float buff[55];
    static int ip = 0,jp = 0,kp;
    int i,seed,get_i_variable();
    static double xva;
    if( ip == 0 && jp == 0)
    {
        /* retrieve the seed from the storage */
        seed = get_i_variable("seed");
        /* initialize when the pointers are both zero   */
        for( ip=0; ip < 55; ip++)
        { seed = (seed*2349+14867)%32767;
            buff[ip] = (float)seed/32767.;
            if( buff[ip] > 1.) buff[ip] = buff[ip]-1.;
            if( buff[ip] < 0.) buff[ip] = buff[ip]+1.;
        }
        ip = 24; jp=55-ip; kp = 0;
    }
    i = kp;
    xva = buff[jp]+buff[ip];
    /* the following line can be non-deteriministic iff xva is
    the same precision as buff */
    if( xva > 1.) xva = xva -1.;
    buff[kp] = xva;
    kp = (kp+1)%55;
    ip = (ip+1)%55;
    jp = (jp+1)%55;
    return buff[i];
}
/* randg()
* return random guassian with unit sigma 
*
* use the polar method (see knuth) with a slight twist.
*
* rwh 5792
*/
float randg()
{
    float randf();
    float x1,x2,norm;

    norm = 2.;
    while( norm > 1.)
    {
        x1 = 2.*randf()-1; x2 = 2.*randf()-1;
        norm = x1*x1 + x2*x2;
    }

    if( norm < 1.e-9) norm = 1.e-9;
    return x1*sqrt( -2.*log(norm)/norm );
}
/* rand3( float *x,*y,*z)
*
* return a random unit three vector
* use knop's method (knuth v.2,p131)
*
* rwh 5/7/92
*/
void rand3( x,y,z)
float *x,*y,*z;
{
    float randf();
    float alpha,norm,x1,x2;
    norm = 2.;
    while( norm > 1.)
    {
        x1 = 2.*randf()-1; x2 = 2.*randf()-1;
        norm = x1*x1 + x2*x2;
    }
    /*	alpha = 2.*sqrt(1.-norm);
    	*x = x1*alpha;
    	*y = x2*alpha;
    	*z = 2.*norm-1.;
    */
    *x = x1; *y = x2;
    norm = sqrt(1.-norm);
    *z = norm;
    if( randf() < 0.5) *z = -norm;
}
