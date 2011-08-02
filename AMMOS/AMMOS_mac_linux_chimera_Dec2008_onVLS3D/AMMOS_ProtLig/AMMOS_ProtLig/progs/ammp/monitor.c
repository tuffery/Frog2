/* monitor.c
*
*  routine to monitor energy and force for AMMP.
*
*  monitors the potential due to each kind of potential used
*
*  reports kinetic energy and total potential and action (T-V)
*
*  reports maximum l_infinity force
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

void monitor( vfs,ffs,nfs,op )
int  (*vfs[])(),(*ffs[])();
int nfs;
FILE *op;
{
    /* block of function used in eval()
    *   only the v_stuff are needed
    */
    int v_bond(),f_bond(),v_angle(),f_angle();
    int v_mmbond(),f_mmbond(),v_mmangle(),f_mmangle();
    int f_c_angle(),v_c_angle();
    int v_periodic(),f_periodic();
    int v_nonbon(),f_nonbon(),v_torsion(),f_torsion();
    int v_screen(),f_screen();
    int v_hard(),f_hard();
    int atom(),bond(),angle(),torsion();
    int v_hybrid(),f_hybrid();
    int morse(),v_morse(),f_morse();
    int restrain(),v_restrain(),f_restrain();
    int tether(),v_tether(),f_tether();
    int u_v_nonbon(), u_f_nonbon();
    int v_abc(),f_abc();
    int v_noel(),f_noel();
    int v_debye(),f_debye();
    int v_shadow(),f_shadow();
    int v_fourd(),f_fourd();
    int v_step(),f_step();
    int v_swarm(),f_swarm();

    int a_number();
    float mxdq,get_f_variable();

    float V,T,vt;
    ATOM *ap,*a_next();
    int ifs,a_f_zero();

    if( a_number() < 1)
    { aaerror(" no atoms defined - nothing to calculate \n"); return;}

    V = 0.; T = 0.;
    a_f_zero();
    for( ifs = 0; ifs < nfs; ifs++ )
    {
        vt = 0.;
        (*vfs[ifs])(&vt,0.);
        mxdq = get_f_variable("mxdq");
        set_f_variable("mxdq",100.);
        (*ffs[ifs])(0.);
        set_f_variable("mxdq",mxdq);

        V += vt;
        if( vfs[ifs] == v_bond)
        { fprintf( op," %f bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_mmbond)
        { fprintf( op," %f mm bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_c_angle)
        { fprintf( op," %f cangle energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_mmangle)
        { fprintf( op," %f mm angle energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_angle)
        { fprintf( op," %f angle energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_abc)
        { fprintf( op," %f abc energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_noel)
        { fprintf( op," %f noel energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_screen)
        { fprintf( op," %f screened non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == u_v_nonbon)
        { fprintf( op," %f non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_nonbon)
        { fprintf( op," %f non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_hard)
        { fprintf( op," %f hard shell non-bonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_torsion)
        { fprintf( op," %f torsion energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_hybrid)
        { fprintf( op," %f hybrid energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_periodic)
        { fprintf( op," %f periodic nonbonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_tether)
        { fprintf( op," %f tether restraint energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_restrain)
        { fprintf( op," %f restraint bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_fourd)
        { fprintf( op," %f 4-d to 3-d projection energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_debye)
        { fprintf( op," %f debye screened nonbonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_shadow)
        { fprintf( op," %f shadow nonbonded energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_morse)
        { fprintf( op," %f morse bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_step)
        { fprintf( op," %f step pseudo-bond energy\n",vt); goto DONE;}
        if( vfs[ifs] == v_swarm)
        { fprintf( op," %f swarm energy\n",vt); goto DONE;}
        fprintf( op," %f unknown potential type\n",vt);
DONE:
        /* next statement is needed because cannot have a label at an end loop */
        vt = 0.;
    }
    fprintf( op," %f total potential energy\n",V);
    ifs = -1;
    while( (ap = a_next(ifs)) != NULL)
    {
        ifs = 1;
        T += ap->vx*ap->vx*ap->mass;
        T += ap->vy*ap->vy*ap->mass;
        T += ap->vz*ap->vz*ap->mass;
    }
    T = T*.5/4.184/1000/1000;
    fprintf( op," %f total kinetic energy\n",T);
    fprintf( op," %f total energy\n",T+V);
    fprintf( op," %f total action\n",T-V);
    set_f_variable("etotal", T+V);
    set_f_variable("epot", V);
    /* end of routine */
}

