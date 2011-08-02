/* just the main routine of ammp
*  read_eval_do and eval are cut out into eval.c and
* used there.
*
* this allows ammp to be called directly from a program
* by making calls to eval or read_eval_do 
*
*  RWH 8/13/93  
*/
/* ammp.c
* Another Molecular Mechanics Program
*
*  this essentially runs the intermediate code for 
*  a molecular mechanics program
*
* instructions are of the form
*  ident <parameters> ;
*  # <stuff> ; is a comment
*  most instructions can be nested, but NOT loop<if> and labels 
*
*  allowed idents
*
*	atom   - atom record
*	bond   - bond record
*       morse  - morse record
*	angle  - angle record
* 	torsion - torsion record
*       hybrid  - hybrid (pyramid height) record
*	abc    - angle bond correlation record
*       velocity  - velocity record
*	read <file>  open and read from file untill done
*	output <file> <vers>  open and use for output file
*	dump <atom,bond,angle,abc,torsion,hybrid,morse,pdb,variable,velocity,force> 
*                         write out the results
*	analyze ilow,ihigh  write out the errors in the current potential for atoms 
*				ilow to ihigh. if ilow > ihigh ilow to ilow
*	close  		close the current output file if not stdout
*	steep  niter,toler   steepest descents
*	bfgs  niter,toler  bfgs quasi newton 
*	cngdel  niter,ncut,toler  conjugate del  
*	trust   niter,dtoler,toler   trust optimizer
*	echo <off>   echo to the user (turn off when dumping !!)
*       use  < none,bond,angle,abc,torsion,nonbon,morse,restrain,tether,periodic,
*              mmbond,mmangle,cangle>
*               flag on potentials
*       restrain    - restrain a distance 
*	tether      - tether an atom to a positon
*          tether serial fk x y z
*          tether all fk x y z  do all of them
*	tailor  qab   number q a b  - set the qab parameters of an atom
*       tailor  exclude  number number  - add an interaction to the nonbon exclude list
*       tailor  include number number  - delete an interaction from the nonbon exclude list
*	setf name value  set a float into the variable store
*       seti name value   set an int into the variable store
*       loopi label init max delta  loop to label while init < max integer vers.
*       loopf label init max delta  loop to label while init < max float vers.
*       label:    
*	monitor    find potential energy and kinetic energy, and calculate the forces
*       v_maxwell  temperature,dx,dy,dz
*	v_rescale   temperature
*       verlet       nstep,dtime (dtime is in m/s = .01A/ps)
*       pac          nstep,dtime (dtime is in m/s = .01A/ps)
*       tpac          nstep,dtime,Temp (dtime in m/s = .01A/ps,1fs = .00001)
*       hpac          nstep,dtime,Htarget (dtime in m/s = .01A/ps,1fs = .00001)
*       pacpac       nstep,dtime (dtime is in m/s = .01A/ps)
*	dipole first,last  calculate the dipole moment for atoms first to last
*                          assumes sequential atom numbers...
*	tgroup id serial1 serial2 serial3 serial4 base number
*            define a tgroup( torsion by serial numbers) base = zeropoint
*	     number == number of steps.  The group of atoms is everything bonded to 
*	      serial3 that isn't serial 2.
*	tsearch id id id id (up to 8  - terminated by 0 or ; ) 
*             search the tgroups defined
*
* 	mompar  serial,chi,jaa  add electronegativity and self colomb to atom serial
*	momadd  serial serial  adds atoms to the MOM stack( can just be called with one)
*       mom   tq, niter   solves current mom stack for charges  
*			tq = total charge, niter = number of iterations (20 default)
*
*
*	math routines  see math.c
*		add a b ;
*		sub a b ;
*		mul a b;
*		div a b;
*		nop a;  these routines can work with atomic parameters 
*		mov a b;  variables, and imeadiate values.
*		max a b;
*		min a b;
*		randf a ;
*
*	  serial a i atomid;  put the serial number or residue i, atom atomid
*                   into a
*	index a i;  put the serial number of the ith atom into a;
*
*        je a b label: ;   jump a == b
*        jl a b label: ;   jump a < b
*        jg a b label: ;   jump a > b
*	jes a string label: ; dump to label if a->name == string
*	jnes a string label: ; dump to label if a->name != string
*           jumps are restricted to the current file
*
*	exit         - exit the routine - in case EOF is not defined
*
*
*	others like fix,and... TBD
*   first nonblank == '#' is a comment and the line is skipped 
*/
/*
*  copyright 1992,1993 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/
#define ANSI 1
#define MAXTOKEN 10 
#define TOKENLENGTH 80 
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
void main()
{
    void read_eval_do();
    /* set some defaults */
    set_f_variable( "mxdq",  .05);
    /* mxcut of 6 requires NCLOSE == 100 */
    /* mxcut of 8 requires NCLOSE == 200 */
    set_f_variable( "mxcut",  6.);
    set_i_variable( "nostep", 1);
    /* read_eval_do is called this way so it can recurse */
    read_eval_do( stdin,stdout);
}

