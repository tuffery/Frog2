/* richard.c
*
*  thermodynamic path integrals using the
*  Feynmann Weiner Kac approach (actually the Wiener version)
*
*  trying the complex Feynman form
*  the occupancy it the modulus
*
*
*  basically:
*
* follow the integral of the exp(lagrangian)
* over the path
*
*  exp( - ( 0.5 m Delta((v )^2)dt - Delta(V) dt) / (kT dt) )
*  where v is velocity, V potential energy and everything
* else is obvious
*
*  Fixed version
*
*  normalize by the mean kinetic energy
*
*   exp( - INtegral(V dt)/ (1.5*NRT *delta_time) )
*   this is correct and cancels out the energy units
*  which should be exp( -Integral/(planck' constant))
*
*  we're just taking the probability relative to a normal kinetic
*  energy distribution
*
*
*  this applies to a global occupancy as we can't
*
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


int richard( op,echo, vfs,ffs,nfs,who,niter,p1,p2,p3)
FILE *op;
int echo;
int (*vfs[])(),(*ffs[])(),nfs;
int niter;
char *who;
float p1,p2,p3;
{
	int pac(),pacpac(),tpac(),hpac(),ptpac(),verlet();
	int (*touse)();
	int i,j;
	float V,Vold,T,K,Kold;
	int set_f_variable();
	float get_f_variable(),o;
	float oreal,oimag;
	float dL,dreal,dimag;
	ATOM *ap,*a_next();
	int numatm, a_number();
/* figure out which  function to use */
	touse = NULL;
	if( strcmp( who,"pac") == 0){  touse = pac;}
	else if( strcmp(who,"pacpac") == 0){touse = pacpac;}
	else if( strcmp(who,"tpac") == 0){touse = tpac;}
	else if( strcmp(who,"hpac") == 0){touse = hpac;}
	else if( strcmp(who,"ptpac") == 0){touse = ptpac;}
	else if( strcmp(who,"verlet") == 0){touse = verlet;}
	if( touse == NULL ){
	aaerror("usage richard <pac,pacpac,tpac,hpac,ptpac> nstep \n");
	return 1;
	}

	T = p2; if( touse == ptpac) T = p3; if( T <= 0.00001) T = 300.;
	V= 0.;
	for( i=0; i<nfs; i++)
		(*vfs[i])(&V,0.);
	numatm = a_number();
	o = get_f_variable("occup");
	if( o <= 1.e-10) o = 1.0;
	V= 0.;
	for( i=0; i<nfs; i++)
		(*vfs[i])(&V,0.);

	for( i=0; i< niter; i++)
	{
		if( touse == pac || touse == pacpac || touse == verlet)
		{ (*touse)(ffs,nfs,1,p1);}
			else if( touse == tpac)
			{(*touse)(ffs,nfs,1,p1,p2);}
			else if( touse == hpac)
			{ (*touse)(ffs,vfs,nfs,1,p1,p2);}
			else 
			{ (*touse)(ffs,nfs,1,p1,p2,p3);}
	Vold = V;
	V= 0.;
	for( j=0; j<nfs; j++)
		(*vfs[j])(&V,0.);

		dL = exp( -(V-Vold)*p1/ ( 1.5*1.987*T*0.001*p1*numatm) );
		o *= dL;
		if( echo) fprintf(op,"%f %f %f %f\n", i*p1, V, dL, o);
		
	}/* end of for(i) */


	set_f_variable( "occup",o);
	return 0;
}/* end of richard */
