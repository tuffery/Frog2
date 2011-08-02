/* atoms.c
*
* collection of routines to service atom memory storage
*
* POOP (Poor-mans Object Oriented Programming) using scope rules
*
* these routines hold a data base (in terms of array indeces)
* of atoms, with the associated forces, and misc consts
*
* routines
*  atom - adds an atom to the table
*  a_m_serial returns pointer to ATOM structure for matching serial
*  a_next  gets next atom in line
*  a_f_zero  zeros force (fx..) entries
*  a_d_zero  zeros dx entries
*  a_v_zero zeros velocity entries 
*  a_number  returns number of atoms
* (this could be table driven but what the hell memories cheap)
*
*  a_readvelocity( serial,vx,vy,vz) sets the velocities
*  dump_atom,dump_velocity,dump_force dump the information
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
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
/* ATOM structure contains a serial number for indexing into
* arrays and the like (a Hessian)
* but otherwise is self-contained. Note the hooks for Non-bonded potentials
*/
#define ALONG sizeof(ATOM)

ATOM *first = NULL;
ATOM *last = NULL;
static int atomNUMBER=0,atomUPDATE=0;
/* function atom adds an atom to the atom list
* returns 1 if ok
* returns 0 if not
* passed the atomic data (at least the initial values )
* allocates the new memory, initializes it and
* returns
*/
int atom(x,y,z,serial,q,a,b,mass,name )
float x,y,z,q,a,b,mass;
int serial;
char *name;
{
    int i;
    static int highest=-1,lowest=-1;
    ATOM *new, *a_m_serial();
    char *ctemp;
    new = NULL;
    if( highest >= serial && lowest <= serial) new = a_m_serial( serial);
    if( new == NULL)
    {
        if( ( new = malloc( ALONG ) ) == NULL)
        {
            return 0;
        }
        new ->dontuse = 0;
        new ->dontuse = 1;
        new ->excluded[0] = new;
        new -> active = 1;
        for( i=0; i< NEXCLUDE; i++)
            new->exkind[i] = 0;
        new->next = NULL;
    }
    /* initialize the pointers */
    if( first == NULL){ first = new;
        highest = serial; lowest = serial; }
    if( last == NULL) last = new;
    new ->x = x;
    new ->y = y;
    new ->z = z;
    new->w = 0.;
    new ->fx = 0.;
    new ->fy = 0.;
    new ->fz = 0.;
    new->fw = 0.;
    new ->dx = 0.;
    new ->dy = 0.;
    new ->dz = 0.;
    new->dw = 0.;
    new ->vx = 0.;
    new ->vy = 0.;
    new ->vz = 0.;
    new->vw = 0.;
    new ->jaa = -1;	new->chi = -1;
    new->rdebye = -1.;
    /* for reinterpolates */
    new ->px = 10e10;
    new ->py = 10e10;
    new ->pz = 10e10;
    new ->pw = 10e10;
    /*
    	new ->qx = 10e10;
    	new ->qy = 10e10;
    	new ->qz = 10e10;
    */
    new ->q = q;
    new ->a = a;
    new ->b = b;
    new ->serial = serial;
    new ->mass = mass;
    for(i=0;i<8; i++)
    {
        new->name[i] = *name;
        new->name[i+1] = '\0';
        if( *name == '\0') break;
        name++;
    }
    if( new->next == NULL)
    {
        new -> next = new;
        last -> next = new;
        last = new;
    }
    atomUPDATE = 1;
    if( highest < serial ) highest = serial;
    if( lowest > serial ) lowest = serial;
    return 1;
}
/* function a_number()
* returns number of atoms defined
*  this is just atomNUMBER if atomUPDATE == 0 
*  other wise just figure it out 
*/
int a_number()
{
    ATOM *ap;
    if( atomUPDATE )
    {
        atomUPDATE = 0;
        atomNUMBER = 0;
        if( first == NULL ) return 0 ;
        ap = first;
        while(1)
        {
            if( ap->next == NULL) break;
            atomNUMBER++;
            if( ap->next == ap ) break;
            ap = ap->next;
        }
    }
    return atomNUMBER;
}
/* function a_m_serial( serial )
* returns NULL on error or returns the address of the ATOM
* which matches serial
* cute?
*/
ATOM *a_m_serial( serial )
int serial;
{
    static ATOM *ap = NULL;
    static ATOM *lastmatched = NULL;
    int i , n, a_number();
    if( atomUPDATE) n= a_number();
    else n = atomNUMBER;

    ap = first; /* static pointer is hook for more efficient search */
    if( ap == NULL) return NULL;
    if( lastmatched == NULL ) lastmatched = first;

    if( serial == lastmatched->serial) return lastmatched;
    if( serial > lastmatched->serial) ap = lastmatched;
    for( i=0; i< n; i++ )
    {
        if( ap-> serial == serial) {lastmatched = ap;return ap;}
        if( ap == ap->next)ap = first ;
        else ap = ap->next;
    }
    return NULL;
}
/* function a_next( flag )
* returns NULL on error or last atom
* then steps to the next
* cute?
* flag <= 0 starts it off
*/
ATOM *a_next( flag )
int flag;
{
    static ATOM *ap = NULL;
    if( ap == NULL) ap = first ;
    if( ap == NULL) return NULL;
if( flag <= 0){ ap = first; return ap;}
    if( ap == ap->next) return NULL;
    ap = ap->next;
    return ap;
}
/* function a_f_zero()
* zeros the forces in each atom element
*/
/* return is 0 on error 1 iff OK */
int a_f_zero()
{
    ATOM *ap;
    ap = first;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> fx = 0.; ap -> fy = 0.; ap -> fz = 0.;
        ap->fw = 0.;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_d_zero()
* zeros the dx,dy,dz storage for each atom element 
*/
/* return is 0 on error 1 iff OK */
int a_d_zero()
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> dx = 0.; ap -> dy = 0.; ap -> dz = 0.;
        ap->dw = 0.;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_g_zero()
* zeros the velocities in each atom element
*/
/* return is 0 on error 1 iff OK */
int a_g_zero()
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> gx = 0.; ap -> gy = 0.; ap -> gz = 0.;
        ap->gw = 0.;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_v_zero()
* zeros the velocities in each atom element
*/
/* return is 0 on error 1 iff OK */
int a_v_zero()
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> vx = 0.; ap -> vy = 0.; ap -> vz = 0.;
        ap->vw = 0.;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_inc_f( lambda )
*  moves the atoms lambda times the forces
*/
int a_inc_f( lambda )
float lambda;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> x += ap->fx*lambda;
        ap->y += ap->fy*lambda; ap->z += ap->fz*lambda;
        ap->w += ap->fw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_inc_d( lambda )
*  moves the atoms lambda times the dx
*/
int a_inc_d( lambda )
float lambda;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> x += ap->dx*lambda;
        ap->y += ap->dy*lambda; ap->z += ap->dz*lambda;
        ap -> w += ap->dw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_inc_v( lambda )
*  moves the atoms lambda times the velocities
*/
int a_inc_v( lambda )
float lambda;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> x += ap->vx*lambda;
        ap->y += ap->vy*lambda; ap->z += ap->vz*lambda;
        ap -> w += ap->vw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_ftodx( lambda )
*  moves the atom force components to the dx slots
* with a scale factor.
*/
int a_ftodx( lambda,lamold )
float lambda,lamold;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap -> dx =ap->dx*lamold+ ap->fx*lambda;
        ap->dy =ap->dy*lamold+ ap->fy*lambda;
        ap->dz =ap->dz*lamold+ ap->fz*lambda;
        ap->dw =ap->dw*lamold+ ap->fw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_ftogx( lambda,lamold )
*  moves the atom force components to the gx slots
* with a scale factor.
*/
int a_ftogx( lambda ,lamold)
float lambda,lamold;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap->gx = ap->gx*lamold +  ap->fx*lambda;
        ap->gy = ap->gy*lamold +  ap->fy*lambda;
        ap->gz = ap->gz*lamold +  ap->fz*lambda;
        ap->gw = ap->gw*lamold +  ap->fw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_ftovx( lambda,lamold )
*  moves the atom force components to the vx slots
* with a scale factor.
*/
int a_ftovx( lambda ,lamold)
float lambda,lamold;
{
    ATOM *ap;
    ap = first;
    if( ap == NULL) return 0;
    while(1)
    {
        if( ap->next == NULL) return 0;
        ap->vx = ap->vx*lamold +  ap->fx*lambda;
        ap->vy = ap->vy*lamold +  ap->fy*lambda;
        ap->vz = ap->vz*lamold +  ap->fz*lambda;
        ap->vw = ap->vw*lamold +  ap->fw*lambda;
        if( ap == ap->next) return 1;
        ap = ap->next;
    }
}
/* function a_max_f()
* returns the maximum l2 metric of a force on an atom
*/
float a_max_f()
{
    float l2norm ,l2max;
    ATOM *ap;
    ap = first;
    l2max = -1.;
    if( ap == NULL) return l2max;
    while(1)
    {
        if( ap->next == NULL) return l2max;
        l2norm = ap->fx*ap->fx;
        l2norm += ap->fy*ap->fy;
        l2norm += ap->fz*ap->fz;
        l2norm += ap->fw*ap->fw;
        if( l2norm > l2max )
            l2max = l2norm;
        if( ap == ap->next) return l2max;
        ap = ap->next;
    }
}
/* function a_max_d()
* returns the maximum l2 metric of a displacement of an atom
*/
float a_max_d()
{
    float l2norm ,l2max;
    ATOM *ap;
    ap = first;
    l2max = -1.;
    if( ap == NULL) return l2max;
    while(1)
    {
        if( ap->next == NULL) return l2max;
        l2norm = ap->dx*ap->dx;
        l2norm += ap->dy*ap->dy;
        l2norm += ap->dz*ap->dz;
        l2norm += ap->dw*ap->dw;
        if( l2norm > l2max )
            l2max = l2norm;
        if( ap == ap->next) return l2max;
        ap = ap->next;
    }
}
/* function a_l2_f(  )
*  return l2 norm of the forces
*/
float a_l2_f(  )
{
    ATOM *ap;
    float l2;
    ap = first;
    if( ap == NULL) return 0.;
    l2 = 0.;
    while(1)
    {
        if( ap->next == NULL) return -l2;
        l2 += ap->fx*ap->fx ;
        l2 += ap->fy*ap->fy ;
        l2 += ap->fz*ap->fz ;
        l2 += ap->fw*ap->fw ;
        if( ap == ap->next) return l2;
        ap = ap->next;
    }
}
/* function a_l2_g(  )
*  return l2 norm of the velocities
*/
float a_l2_g(  )
{
    ATOM *ap;
    float l2;
    ap = first;
    if( ap == NULL) return 0.;
    l2 = 0.;
    while(1)
    {
        if( ap->next == NULL) return -l2;
        l2 += ap->gx*ap->gx ;
        l2 += ap->gy*ap->gy ;
        l2 += ap->gz*ap->gz ;
        l2 += ap->gw*ap->gw ;
        if( ap == ap->next) return l2;
        ap = ap->next;
    }
}
/* function a_l2_v(  )
*  return l2 norm of the velocities
*/
float a_l2_v(  )
{
    ATOM *ap;
    float l2;
    ap = first;
    if( ap == NULL) return 0.;
    l2 = 0.;
    while(1)
    {
        if( ap->next == NULL) return -l2;
        l2 += ap->vx*ap->vx ;
        l2 += ap->vy*ap->vy ;
        l2 += ap->vz*ap->vz ;
        l2 += ap->vw*ap->vw ;
        if( ap == ap->next) return l2;
        ap = ap->next;
    }
}
/* function a_l2_d(  )
*  return l2 norm of the dx terms
*/
float a_l2_d(  )
{
    ATOM *ap;
    float l2;
    ap = first;
    if( ap == NULL) return 0.;
    l2 = 0.;
    while(1)
    {
        if( ap->next == NULL) return -l2;
        l2 += ap->dx*ap->dx ;
        l2 += ap->dy*ap->dy ;
        l2 += ap->dz*ap->dz ;
        l2 += ap->dw*ap->dw ;
        if( ap == ap->next) return l2;
        ap = ap->next;
    }
}
/* routine dump_atoms
* this function outputs the atomic parameters
* and does it in a simple form
*  atom x,y,z,serial,name,q,a,b,mass
* where atom is the string atom
* the rest is just free format
*/
void dump_atoms( where )
FILE *where;
{
    ATOM *a,*ap;
    ATOM *bonded[20];
    int i ,j ;
    void dump_excludes();
    a = first;

    if( a == NULL) return;
    while( (a->next != a)  )
    {
        if( a->next == NULL) return;
        fprintf( where,"atom %f %f %f %d %s %f %f %f %f \;\n",
                 a->x,a->y,a->z,a->serial,a->name,a->q,a->a,a->b,
                 a->mass );
        if( a->chi > 0 && a->jaa > 0)
            fprintf( where,"mompar %d %f %f;\n",
                     a->serial,a->chi,a->jaa);
        if( !a->active) fprintf(where," inactive %d ;\n", a->serial);
        if( a->rdebye > 0.)fprintf(where,"rdebye %d %f;\n",a->serial,a->rdebye);
        a = a->next;
    }
    if( a->next == NULL) return;
    fprintf( where,"atom %f %f %f %d %s %f %f %f %f \;\n",
             a->x,a->y,a->z,a->serial,a->name,a->q,a->a,a->b,
             a->mass );
    if( a->chi > 0 && a->jaa > 0)
        fprintf( where,"mompar %d %f %f;\n",
                 a->serial,a->chi,a->jaa);
    if( !a->active) fprintf(where," inactive %d ;\n", a->serial);
    if( a->rdebye > 0.)fprintf(where,"rdebye %d %f;\n",a->serial,a->rdebye);
    /* all of the atoms have been dumped so now dump the excludes */
    dump_excludes(where);
} /* end of dump atoms */
/* routine dump_excludes( FILE *where )
*  write out all of the excludes */
void dump_excludes( where )
FILE *where ;
{
    ATOM *a,*ap,*a_next();
    int istailored;
    int get_i_variable();
    int numatm,a_number();
    int i,j;

    /*
    	istailored = 0;
    	istailored = get_i_variable("numtail");
    	if( istailored <= 0 ) return;
    */

    numatm = a_number();
    if( numatm <= 0 ) return;

    for( i=0; i< numatm; i++)
    {
        a = a_next(i);
        for( j=0; j< a->dontuse; j++)
        {
            if( a->exkind[j] > 0) {
                ap = a->excluded[j];
                if( ap->serial > a->serial)
                    fprintf( where," tailor exclude %d %d ;\n",
                             a->serial, ap->serial);
                /*
                		if( ap->serial < a->serial)
                		fprintf( where," tailor exclude %d %d ;\n",
                		a->serial, ap->serial);
                */
            }
        }
    }
}/* end of dump_excludes */
/* routine dump_velocity
* this function outputs the atomic parameters
* and does it in a simple form
*  velocity serial vx,vy,vz
* where atom is the string atom
* the rest is just free format
*/
void dump_velocity( where )
FILE *where;
{
    ATOM *a;
    a = first;
    if( a == NULL) return;
    while( (a->next != a)  )
    {
        if( a->next == NULL) return;
        fprintf( where,"velocity %d %f %f %f \;\n",
                 a->serial,a->vx,a->vy,a->vz );
        a = a->next;
    }
    if( a->next == NULL) return;
    fprintf( where,"velocity %d %f %f %f \;\n",
             a->serial,a->vx,a->vy,a->vz );
}
/* int a_readvelocity( serial,vx,vy,vz)
*  update the velocity field of the atom structure
*/
int a_readvelocity( serial,vx,vy,vz)
int serial;
float vx,vy,vz;
{
    ATOM *ap,*a_m_serial();
    ap = a_m_serial( serial);
    if( ap == NULL) return 0;
    ap ->vx = vx;
    ap ->vy = vy;
    ap ->vz = vz;
    return 1;
}
/* routine dump_force
* this function outputs the atomic parameters
* and does it in a simple form
*  force serial x,y,z,fx,fy,fz
* where atom is the string atom
* the rest is just free format
*/
void dump_force( where )
FILE *where;
{
    ATOM *a;
    a = first;
    if( a == NULL) return;
    while( (a->next != a)  )
    {
        if( a->next == NULL) return;
        fprintf( where,"force %d %f %f %f %f %f %f \;\n",
                 a->serial,a->x,a->y,a->z,a->fx,a->fy,a->fz );
        a = a->next;
    }
    if( a->next == NULL) return;
    fprintf( where,"force %d %f %f %f %f %f %f \;\n",
             a->serial,a->x,a->y,a->z,a->fx,a->fy,a->fz );
}
/* routine dump_pdb
* this function outputs the atomic parameters
* in pdb format
*
* res_mod is used when encoding atoms to serial
*  code is  (residue_number-1) * res_mod + (atom_number-1)
*   the -1 is FORTRAN style uugh!!!  (i.e.  zero is not allowed number)
*
*  residue name and atom name are dot coded as
*   wat.oh
*
* other codes are possible, but this is what we chose for the moment
*/
void dump_pdb( where ,res_mod)
FILE *where;
int res_mod;
{
    ATOM *a;
    char *np,resid[5],atid[5];
    float o,get_f_variable();
    int i,ires,iatom;
    a = first;
    iatom = 0;
    if( a == NULL) return;
    if( res_mod == 0 )
    {
        aaerror( " need a non-zero residue modulus in dump_pdb\n");
        return ;
    }
    o = get_f_variable("occup");
    if( o >= 0)
        fprintf(where,"REMARK WKF occupancy is %f\n",o);
    while( (a->next != a)  )
    {
        if( a->next == NULL) return;
        iatom++;
        /*		ires = a->serial/res_mod +1 ; */
        ires = a->serial/res_mod  ;
        np = a->name;
        while( strcmp(np,"sna.rkq") == 0)
        { a= a->next;
            if( a->next == NULL ) return; ires = a->serial/res_mod;
            np = a->name; }
        for( i=0; i<5;i++)
            { if(*np != '.')
            { if( islower(*np)) {resid[i] = toupper(*np);}
                else{ resid[i] = *np;}}
            else{ resid[i] = '\0'; break; }
            if( *np == '\0') break;
            np++;
        }
        if( *np == '.') np++;
        for( i=0; i<5;i++)
            { if(*np != '.')
            { if( islower(*np)){ atid[i] = toupper(*np);}
                else{ atid[i] = *np;}}
            else{ atid[i] = '\0'; break; }
            if( *np == '\0' ) break;
            np++;
        }
        /* brookhaven format ,(sort of) */
        /*
        	fprintf(where,
        	"ATOM  %5d %-4s%c%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        	iatom,atid,' ',resid,ires,a->x,a->y,a->z,1.,10.);
        */
        if( atid[0] == 'H')
            fprintf(where,
                    "ATOM  %5d %-4s%c%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    iatom,atid,' ',resid,ires,a->x,a->y,a->z,1.,10.);
        else
            fprintf(where,
                    "ATOM  %5d  %-4s%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    iatom,atid,resid,ires,a->x,a->y,a->z,1.,10.);

        a = a->next;
    }
    if( a->next == NULL) return;
    iatom++;
    ires = a->serial/res_mod ;
    np = a->name;
    if( strcmp(np,"sna.rkq") != 0)
    {
        for( i=0; i<5;i++)
            { if(*np != '.')
            { if( islower(*np)){ resid[i] = toupper(*np);}
                else{ resid[i] = *np;}}
            else{ resid[i] = '\0'; break; }
            if( *np == '\0') break;
            np++;
        }
        if( *np == '.') np++;
        for( i=0; i<5;i++)
            { if(*np != '.')
            {if( islower(*np)){ atid[i] = toupper(*np);}
                else{ atid[i] = *np;}}
            else{ atid[i] = '\0'; break; }
            if( *np == '\0' ) break;
            np++;
        }
        /* brookhaven format ,(sort of) */
        /*
        	fprintf(where,
        	"ATOM  %5d %-4s%c%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
        	iatom,atid,' ',resid,ires,a->x,a->y,a->z,1.,10.);
        */
        if( atid[0] == 'H')
            fprintf(where,
                    "ATOM  %5d %-4s%c%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    iatom,atid,' ',resid,ires,a->x,a->y,a->z,1.,10.);
        else
            fprintf(where,
                    "ATOM  %5d  %-4s%-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    iatom,atid,resid,ires,a->x,a->y,a->z,1.,10.);
    }/* end of if for SNARK.Q */
    fprintf(where,"END   \n");
}
/* function a_pr_beta()
* a_pr_beta() returns the Poliak Ribeire beta
* for conjugate gradients
*
*/
float a_pr_beta()
{
    float a,b;
    ATOM *ap;
    ap = first;
    a = 0.; b = 0.;
    if( ap == NULL) return a;
    while(1)
    {
        if( ap->next == NULL) return 0.;
        a += ap->fx*ap->fx;
        a += ap->fy*ap->fy;
        a += ap->fz*ap->fz;
        a += ap->fw*ap->fw;
        a -= ap->gx*ap->fx;
        a -= ap->gy*ap->fy;
        a -= ap->gz*ap->fz;
        a -= ap->gw*ap->fw;
        b += ap->gx*ap->gx;
        b += ap->gy*ap->gy;
        b += ap->gz*ap->gz;
        b += ap->gw*ap->gw;
        if( ap == ap->next)
        {
            if( b  <= 1.e-5) { a = 0.; b = 1.;}
            return a/b;
        }
        ap = ap->next;
    }
}
/* inactivate ( i1,i2 )
*  i2  == 0 just do i 1 
*  else do the range from i1 to i2 
*
* set the active flag in the atoms to 0 and
* turn them off
*/
inactivate_non_zero (i1,i2)
int i1,i2 ;
{
    int upper, lower;
    ATOM *ap,*a_m_serial(),*a_next();
    int i ,numatm,a_number();

    if( i2 == 0 )
    {
        ap = a_m_serial( i1) ;
        if( ap != NULL)
            ap -> active = 0;
        return ;
    }

    upper = i2; lower = i1;
if( i2 < i1 ) { lower = i2; upper = i1;}

    numatm = a_number();
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->serial >= lower && ap->serial <= upper)
            if( ap->x != 0. && ap->y != 0. && ap->z != 0.)
                ap->active = 0;
    }
}
/* inactivate ( i1,i2 )
*  i2  == 0 just do i 1 
*  else do the range from i1 to i2 
*
* set the active flag in the atoms to 0 and
* turn them off
*/
inactivate (i1,i2)
int i1,i2 ;
{
    int upper, lower;
    ATOM *ap,*a_m_serial(),*a_next();
    int i ,numatm,a_number();

    if( i2 == 0 )
    {
        ap = a_m_serial( i1) ;
        if( ap != NULL)
            ap -> active = 0;
        return ;
    }

    upper = i2; lower = i1;
if( i2 < i1 ) { lower = i2; upper = i1;}

    numatm = a_number();
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->serial >= lower && ap->serial <= upper)
            ap->active = 0;
    }
}
/* activate ( i1,i2 )
*  i2  == 0 just do i 1 
*  else do the range from i1 to i2 
*
* set the active flag in the atoms to 1 and
* turn them on
*/
activate (i1,i2)
int i1,i2 ;
{
    int upper, lower;
    ATOM *ap,*a_m_serial(),*a_next();
    int i ,numatm,a_number();

    if( i2 == 0 )
    {
        ap = a_m_serial( i1) ;
        if( ap != NULL)
            ap -> active = 1;
        return ;
    }

    upper = i2; lower = i1;
if( i2 < i1 ) { lower = i2; upper = i1;}

    numatm = a_number();
    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->serial >= lower && ap->serial <= upper)
            ap->active = 1;
    }
}
/* a_inactive_f_zero()
* loop through the atoms and zero the inactive forces 
*/
a_inactive_f_zero()
{
    int i ,numatom,a_number();
    ATOM *ap,*a_next();
    numatom = a_number();
    for( i=0; i< numatom ; i++)
    {
        ap = a_next(i);
        if( ap->active == 0)
        {  ap->fx = 0.; ap->fy = 0.; ap->fz = 0.; ap->fw = 0.; }
    }
}
