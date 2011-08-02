/* analyze.c
*
*  routine to analyze energy and force for AMMP.
*
*  analyzes the potential due to each kind of potential used
*
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
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-bonded potentials
*/

void analyze( vfs,nfs,ilow,ihigh,op )
int  (*vfs[])();
int nfs;
FILE *op;
int ilow,ihigh;
{
    /* block of function used in eval()
    *   only the v_stuff are needed
    */
    int v_bond(),f_bond(),v_angle(),f_angle();
    int v_mmbond(),f_mmbond(),v_mmangle(),f_mmangle();
    int v_periodic(),f_periodic();
    int v_nonbon(),f_nonbon(),v_torsion(),f_torsion();
    int v_ttarget(),f_ttarget();
    int v_screen(),f_screen();
    int atom(),bond(),angle(),torsion();
    int v_hybrid(),f_hybrid();
    int v_step(),f_step();
    int v_swarm(),f_swarm(),a_swarm();
    int morse(),v_morse(),f_morse();
    int restrain(),v_restrain(),f_restrain();
    int tether(),v_tether(),f_tether();
    int u_v_nonbon(), u_f_nonbon();

    int v_abc();

    int v_noel(),a_noel();
    int v_ho_noel();
    int a_bond(),a_mmbond(),a_angle(),a_mmangle();
    int a_nonbon(),a_torsion(),a_hybrid(),a_restrain();
    int a_react();
    int a_ttarget();
    int a_tether();
    int a_screen();
    int a_abc(),a_step();
    int v_image(),a_image();

    float V,vt;
    int ifs;
    int i,j;
    i = ilow;
    j = ihigh;
    if( ihigh < ilow ) j = ilow;
    V = 0.;
    for( ifs = 0; ifs < nfs; ifs++ )
    {
        vt = 0.;
        /*	(*vfs[ifs])(&vt,0.); */
        if( vfs[ifs] == v_bond)
        { a_bond(&vt,0.,i,j,op);fprintf( op," %f bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_mmbond)
        {a_mmbond(&vt,0.,i,j,op); fprintf( op," %f mm bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_mmangle)
        {a_mmangle(&vt,0.,i,j,op); fprintf( op," %f mm angle energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_angle)
        {a_angle(&vt,0.,i,j,op); fprintf( op," %f angle energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_abc)
        {a_abc(&vt,0.,i,j,op); fprintf( op," %f abc energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_noel)
        {a_noel(&vt,0.,i,j,op); fprintf( op," %f noel energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_ho_noel)
        {a_noel(&vt,0.,i,j,op); fprintf( op," %f noel energy\n",vt); goto DONE;}
        if( vfs[ifs] == u_v_nonbon)
        {a_nonbon(&vt,0.,i,j,op); fprintf( op," %f non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_nonbon)
        {a_nonbon(&vt,0.,i,j,op); fprintf( op," %f non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_screen)
        {a_screen(&vt,0.,i,j,op); fprintf( op," %f screened non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_torsion)
        {a_torsion(&vt,0.,i,j,op); fprintf( op," %f torsion energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_hybrid)
        {a_hybrid(&vt,0.,i,j,op); fprintf( op," %f hybrid energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_periodic)
        {  goto DONE;}
        if( vfs[ifs] == v_tether)
        {a_tether(&vt,0.,i,j,op); fprintf( op," %f tether restraint energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_restrain)
        {a_restrain(&vt,0.,i,j,op); fprintf( op," %f restraint bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_morse)
        {  goto DONE;}
        if( vfs[ifs] == v_image)
        { a_image(&vt,0., i,j,op); fprintf(op," %f image energy\n",vt);
            goto DONE; }
        if( vfs[ifs] == v_step)
        { a_step(&vt,0.,i,j,op); fprintf(op,"%f step energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_swarm)
        {a_swarm(&vt,0.,i,j,op); fprintf(op,"%f swarm energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_ttarget)
        {a_ttarget(&vt,0.,i,j,op); fprintf( op," %f torsion target energy\n",vt); goto DONE;}
DONE:
        /* next statement is needed because cannot have a label at an end loop */
        V += vt;
        vt = 0.;
    }
    fprintf( op," %f total potential energy\n",V);
    /* end of routine */
}

