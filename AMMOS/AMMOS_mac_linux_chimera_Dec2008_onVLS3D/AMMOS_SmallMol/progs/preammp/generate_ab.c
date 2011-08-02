/* uberion.c
*
* generate super-radius ions in AMMP format
*
* 
*  use 6,12 potential for the 'hard shell'
*  this may need to be different for 'soft' ions
*
*  a = sqrt( 2*emin*r^6) where emin is the 
*             || of the depth
*  b = sqrt( emin * r^12 )
*
* then the coeffcient is just the product
* of a1*a2 and b1*b2.
*
*
*  these will depend on 'ionic strength'
* and are normalized to 
*  1molar = 1ion/55 waters
*   (actually ) 1000./18 waters.
*
*  the mean volume of a water is
*   1000ml/1000g*18g/M*6.02e-23 * 10^24A^3/ml
*
* = 10.84 A^3
*/
#define WATER 10.84
/*
*  the volume of the uberion is then
*
*   nwater/nion *10.84  596 A^3
*  when converted to radius this is 5.22 A
*  which is not so crazy.
*
*  the depth of the potential is chosen to be the average
*  water vdw term of 0.06 kcal/M
*
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
	char line[BUFSIZ];
	float q,a,b,r,emin;
	float mass;
	float ion;

	printf("the charge\n");
	fgets( line,80,stdin);
	sscanf(line,"%f",&q);
	
	printf(" the radius\n");
	fgets( line,80,stdin);
	sscanf(line,"%f",&r);

	printf(" the potential depth\n");
	fgets( line,80,stdin);
	sscanf(line,"%f",&emin);
	emin = fabs(emin);

	printf(" the ion mass\n");
	fgets( line,80,stdin);
	sscanf(line,"%f",&mass);

/* volume  */
	r = r*r*r;
	r *= r; /* r^6 */
	a = sqrt( 2.*emin*r);
	b = sqrt(emin*r*r);
	printf("atom 0 0 0 myserial uber.ion %f %f %f %f;\n",
		q,a,b,mass); 

}/* end of main */
