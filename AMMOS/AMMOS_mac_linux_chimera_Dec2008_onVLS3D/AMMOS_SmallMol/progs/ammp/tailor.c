/* tailor.c
*
* collection of routines to service atom memory storage
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of atoms, with the associated forces, and misc consts
*
* routines for tailoring nonbonded parameters during a run
*  tailor_qab
*  tailor_include
*  tailor_exclude
*
*/
/*
*  copyright 1992 1993 Robert W. Harrison
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
* but otherwise is self-contained. Note the hooks for Non-bonded potentials
*/
void tailor_qab( who, q,a,b)
int who;
float q, a, b;
{
    ATOM *ap,*a_m_serial();

    ap = a_m_serial( who );
    if( ap == NULL ) {aaerror(" undefined atom in tailor_qab"); return;}
    ap->q = q; ap->a = a; ap->b = b;
    return;
}
void tailor_include( from, to )
int from , to;
{
    ATOM *ap1,*ap2,*a_m_serial();
    int i;
    ap1 = a_m_serial( from );
    if( ap1 == NULL ) {/*aaerror(" undefined atom in tailor_include"); */return;}
    ap2 = a_m_serial( to );
    if( ap2 == NULL ) {/*aaerror(" undefined atom in tailor_include");*/ return;}

    for( i=0; i< ap1->dontuse; i++)
    { if( ap1->excluded[i] == ap2) break; }
    if( i < ap1->dontuse-1 )
    {
        for(; i< ap1->dontuse-1; i++)
        {
            ap1->excluded[i] = ap1->excluded[i+1];
            ap1->exkind[i] = ap1->exkind[i+1];
        }
    }
    ap1->exkind[(ap1->dontuse)] = 0;
    ap1->dontuse -= 1;
    for( i=0; i< ap2->dontuse; i++)
    { if( ap2->excluded[i] == ap1) break; }
    if( i < ap2->dontuse-1 )
    {
        for(; i< ap2->dontuse-1; i++)
        {
            ap2->excluded[i] = ap2->excluded[i+1];
            ap2->exkind[i] = ap2->exkind[i+1];
        }
    }
    ap2->exkind[(ap2->dontuse)] = 0;
    ap2->dontuse -= 1;

}
void tailor_exclude( from, to )
int from , to;
{
    ATOM *ap1,*ap2,*a_m_serial();
    int set_i_variable(),get_i_variable();
    int i;
    ap1 = a_m_serial( from );
    if( ap1 == NULL ) {/*aaerror(" undefined atom in tailor_exclude"); */return;}
    ap2 = a_m_serial( to );
    if( ap2 == NULL ) {/*aaerror(" undefined atom in tailor_exclude");*/ return;}

    if( ap1->dontuse == NEXCLUDE)
    {aaerror(" cannot add an atom to the exclude list "); return; }
    for( i=0; i< ap1->dontuse; i++)
    {  if(ap1->excluded[i] == ap2) goto FOUND; }
    ap1->exkind[(ap1->dontuse)] = 1;
    ap1->excluded[(ap1->dontuse)++] = ap2;
FOUND:
    if( ap2->dontuse == NEXCLUDE)
    {aaerror(" cannot add an atom to the exclude list "); return; }
    for( i=0; i< ap2->dontuse; i++)
    {  if(ap2->excluded[i] == ap1) return; }
    ap2->exkind[(ap2->dontuse)] = 1;
    ap2->excluded[(ap2->dontuse)++] = ap1;

    i = 0;
    i = get_i_variable("numtail");
    i = i + 1;
    set_i_variable("numtail", i );
}

