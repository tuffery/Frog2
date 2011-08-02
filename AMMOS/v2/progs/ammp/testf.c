#include <stdio.h>
#include <math.h>
#include "numeric.h"

main()
{ int i; float x, Fzero(),Fone();

    for( i=0; i< 100; i++)
    {
        x = i;

        printf("%f %f %f\n",x,Fzero(x*x),Fone(x*x));
    }

}

float Fzero( x)
float x;
{
#ifdef notanyoldSGI
    float accum,etox;
    accum = 0.;
    etox = exp(-x);
    /*
            accum = 1./(23.)*(2.*x*accum + etox);
            accum = 1./(21.)*(2.*x*accum + etox);
    */
    accum = 1./(19.)*(2.*x*accum + etox);
    accum = 1./(17.)*(2.*x*accum + etox);
    accum = 1./(15.)*(2.*x*accum + etox);
    accum = 1./(13.)*(2.*x*accum + etox);
    accum = 1./(11.)*(2.*x*accum + etox);
    accum = 1./( 9.)*(2.*x*accum + etox);
    accum = 1./( 7.)*(2.*x*accum + etox);
    accum = 1./( 5.)*(2.*x*accum + etox);
    accum = 1./( 3.)*(2.*x*accum + etox);
    accum = 2.*x*accum + etox;
    return accum;
#else
    if( x < 1.e-7) return 1.;
    x = sqrt(x);
    return .5*ROOTPI/x*erf(x);

#endif

}/* end of Fzero */

/* Fone(x) is -d/dx(Fzero(x)) */
float Fone(x)
float x;
{
    float y ;
    if( x < 1.e-7) return 0.;
    y = sqrt(x);
    y = .5*ROOTPI/y*erf(y);   /* Fzero */
    /* Fone is defined by the recursion
    *  Fzero = (2xFone + exp(-x)) 
    * or 
    *  Fone = 1./(2x)*(Fzero-exp(-x))
    */
    return (y-exp(-x))/(x+x);
} /* end of Fone */


