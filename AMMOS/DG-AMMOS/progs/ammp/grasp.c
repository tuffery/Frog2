/* grasp.c
*  search a space of torsion's using the
*  current potentials and greedy algorithms
*
* SINCE AMMP ONLY HAS ONE COORDINATE SET
* ONLY THE BEST IS KEPT
*
* uses the routines in tset.c and random.c
* these routines require atom dot notation (x.ca)
* and use the pdb mode 
*
*
*/


/*
*  copyright 1993,1995,1996 Robert W. Harrison
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
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"



int grasp( op,echo,vfs,ffs,nfs, nstep,nkeep,nsearch, imin,imax,whom)
FILE *op;
int echo;
int nfs,(*vfs[])(),(*ffs[])();
int nstep,nkeep,nsearch,imin,imax;
char *whom;
{
    int natom,i,j,k,l,lmax,a_number();
    float x,delta,randf();
    int run;
    ATOM *(*apa)[],*a_next(),*ap;
    char *cp;
    float the_best_energy;
    float (*the_best_xyz)[];
    /* defined in tset.c */
    float get_torsion_value();
    int tset_bond_build();
    int set_torsion();
    int grasp_valid_coordinate();

    if( imax < imin ) { j = imax; imax = imin;imin = j;}

    the_best_energy = 1.e10;
    j = imax-imin;
    if( j < 4 || a_number() < 4) {
        aaerror(" you need at least four atoms for grasp");
        return 1==1;
    }
    apa = malloc( j*sizeof( ATOM * ));
    if( apa == NULL){ goto ERROR ; }
    the_best_xyz = malloc( 3*a_number()*sizeof( float ) );
    if( the_best_xyz == NULL ) {goto ERROR ;}
    /* loop through the atoms and find all in the range
    * from imin to imax which are of type *whom
    */
    for( i=0; i< j; i++)
        (*apa)[i] = NULL;
    ap = a_next(-1);
    natom = 0;
    lmax = a_number();
    for( i=0; i< j; i++)
    {
        (*apa)[natom] = NULL;
        for( l=0; l< lmax; l++){
            k = ap->serial/100;
            if( k == i+imin){
                cp = &ap->name[0];
                while(*cp != '.' && *cp != '\0') cp++;
                if(*cp == '.') cp++;
                if( strcmp(cp,whom) == 0 )
                {
                    (*apa)[natom] = ap;
                    natom++;
                    break;
                }
            }/* the residue is the right one */
            if( ap->next == ap || ap == NULL || ap->next == NULL)
                ap = a_next(-1);
            else
                ap = ap->next;
        }/* l */
    }/* i */
    /* now the array apa has the atom pointers, in order and
    * natom is the number of atoms 
    */
    for( i=0; i< nstep; i++)
    {
        for( j=0; j< natom-3; )
        {
            run = 10*randf() + 1;

            /*
            	 delta = 2*3.141592653589793*randf()  ; 
            */
            delta = 3.141592653589793*(2*randf()-1.)/10.  ;
            /*
            	for( l=0; l< run; l++)
            	{
            */
            if( (*apa)[j] == NULL ) break;
            if( (*apa)[j+1] == NULL ) break;
            if( (*apa)[j+2] == NULL ) break;
            if( (*apa)[j+3] == NULL ) break;
            /*
            	x = get_torsion_value( (*apa)[j],(*apa)[j+1],(*apa)[j+2],(*apa)[j+3]);
            	tset_bond_build( (*apa)[j],(*apa)[j+1],(*apa)[j+2],(*apa)[j+3]);
            	set_torsion( (*apa)[j],(*apa)[j+1],(*apa)[j+2],(*apa)[j+3],delta-x);
            */
            tset_bond_build( (*apa)[j],(*apa)[j+1],(*apa)[j+2],(*apa)[j+3]);
            set_torsion( (*apa)[j],(*apa)[j+1],(*apa)[j+2],(*apa)[j+3],delta-x);
            j += run;
            /*
            	j += 1; 
            	if( j > natom-4) break;

            	}*//*l*/
        }/*j*/

        if( !grasp_valid_coordinate()) continue ;
        cngdel( vfs,ffs, nfs, nsearch,nsearch, 0.,0);
        delta = 0.;
        for( j=0; j< nfs; j++)
            (*vfs[j])( &delta,0.);
        if(echo)fprintf(op,"step %d the best energy is %f %f\n",i,
                            the_best_energy,delta);
        if( delta < the_best_energy)
        {
            if(echo)fprintf(op,"step %d the best energy is %f\n",i,delta);
            the_best_energy = delta;
            for( l=0; l < lmax; l++)
            {
                ap = a_next(l);
                (*the_best_xyz)[3*l] = ap->x;
                (*the_best_xyz)[3*l+1] = ap->y;
                (*the_best_xyz)[3*l+2] = ap->z;
            }

        }

    }
    for( l=0; l < lmax; l++)
    {
        ap = a_next(l);
        ap->x = (*the_best_xyz)[3*l];
        ap->y = (*the_best_xyz)[3*l+1];
        ap->z = (*the_best_xyz)[3*l+2];
    }

    free( the_best_xyz);
    free(apa);
    return 1==1;
ERROR: ;
    aaerror("cannot allocate memory in grasp");
    if( the_best_xyz != NULL ) free(the_best_xyz);
    if( apa != NULL ) free(apa);
    return 1==0;
}/* end of grasp */



int grasp_valid_coordinate()
{

    int na,a_number();
    ATOM *ap1,*ap2,*a_next();
    int i,j;
    float x,y,z;

    na = a_number();
    if( na < 1) return ;

    ap1 = a_next(-1);
    ap1 = ap1->next;
    for( i=1; i< na; i++)
    {
        for( j=0; j< i; j++)
        {
            ap2 = a_next(j);

            x = fabs(ap1->x -ap2->x);
            if( x < 1.e-1)
            {
                y = fabs(ap1->y -ap2->y);
                if( y < 1.e-1)
                {
                    z = fabs(ap1->z -ap2->z);
                    if( z< 1.e-1)
                    {
                        return 1==0;
                    }}}
        }
        ap1 = ap1->next;
    }

    return 1==1;
}
