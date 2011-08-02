/* ammp.c
* Another Molecular Mechanics Program
*
*  this essentially runs the intermediate code for 
*  a molecular mechanics program
*
* instructions are of the form
*  ident <parameters> ;
*  # <stuff> ; is a comment
*   "<stuff>" is  a literal string
*  most instructions can be nested, but NOT loop<if> and labels 
*
*  allowed idents
*
*	atom   - atom record
*	bond   - bond record
*       morse  - morse record
*	angle  - angle record
* 	torsion - torsion record
*       ttarget - torsion target i j k l angle fk;
*       swarm - swarm k n (ori end),... n/2 times restrain a bunch of pairs
*                                           to the same distance 
*       clone i ; clone to i
*       statclone; status of clones
*       subclone <list>; average clones in list  subclone ; average all clones
*       nnclone <list>; average clones in list  generate noel terms
*
*       hybrid  - hybrid (pyramid height) record
*	abc    - angle bond correlation record 
*	dill - dill style nonbonded term
*                  i1 i2 i3 angle zero_angle  dr/da dk12/da dk23/da
*	av5  - tetrahedral 'volume' or centroid
*		i1,i2,i3,i4,i5, k, value
*	noel -noe distance restraint
* 	noegen dm,dp;  set the noels to the current geometry
*       central atom  weight exponent;  central force terms for folding
*	rdebye atom value;  context-variable debye radius
*       velocity  - velocity record
*	read <file>  open and read from file untill done
*	output <file> <vers>  open and use for output file
*	dump <atom,bond,angle,abc,ttarget,torsion,hybrid,av5,morse,pdb,variable,velocity,force> 
*                         write out the results
*	analyze ilow,ihigh  write out the errors in the current potential for atoms 
*				ilow to ihigh. if ilow > ihigh ilow to ilow
*	close  		close the current output file if not stdout
*	steep  niter,toler   steepest descents
*	bfgs  niter,toler  bfgs quasi newton 
*	cngdel  niter,ncut,toler  conjugate del  
*	trust   niter,dtoler,toler   trust optimizer
*       polytope imin,imax,niter,vstart,vfinal  polytope of a range
*	rigid imin,imax,niter, vstart,vfinal  polytope rigid body solver
*  ghost id low high;  set up a ghost for energy ghosts are invisible to each other
*  buster id;  use a ghost
*
*	gdock toler,ngene,niter, var_angle,var_xyz,imin,imax;  genetic dock
*                                           routine
*	echo <off>   echo to the user (turn off when dumping !!)
*       use  < none,bond,angle,abc,ttarget,torsion,nonbon,morse,restrain,tether
*		,periodic,mmbond,mmangle,cangle,gauss,screen,debye,shadow,fourd
*		hobond hoangle trace honoel hotether image >  
*               flag on potentials
*       restrain    - restrain a distance 
*	tether      - tether an atom to a positon
*          tether serial fk x y z
*          tether all fk x y z  do all of them
*         
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
*       ppac          nstep,dtime,pressure (dtime in m/s = .01A/ps,1fs = .00001)
*       ptpac          nstep,dtime,pressure,Temp (dtime in m/s = .01A/ps,1fs = .00001)
*       hpac          nstep,dtime,Htarget (dtime in m/s = .01A/ps,1fs = .00001)
*       pacpac       nstep,dtime (dtime is in m/s = .01A/ps)
* 	richard <pac ... verlet> nstep <md parameters>;  WFK integral 
*                                     sets the parameter "occup"
*
*	doubletime   nstep,dlong,dshort,temper  double time scale dynamics
*	dipole first,last  calculate the dipole moment for atoms first to last
*                          assumes sequential atom numbers...
*	tgroup id serial1 serial2 serial3 serial4 base number
*            define a tgroup( torsion by serial numbers) base = zeropoint
*	     number == number of steps.  The group of atoms is everything bonded to 
*	      serial3 that isn't serial 2.
*	tsearch id id id id (up to 8  - terminated by 0 or ; ) 
*             search the tgroups defined
*
*       tset i1 i2 i3 i4 where
*            set the torsion angle defined by i1...i4 to where
*            unlike tgroup,tsearch only one angle at a time, and
*            no limit to the number of atoms rotated
*	tmin i1 i2 i3 i4 nstep
*            search the torsion angle  i1...i4 for the minimum
*            energy in nsteps
*	tmap i1 i2 i3 i4 j1 j2 j3 j4 ni nj;  map in ni nj steps
*              the i j atoms over all 360 degrees;
* 
* 	mompar  serial,chi,jaa  add electronegativity and self colomb to atom serial
*	momadd  serial serial  adds atoms to the MOM stack( can just be called with one)
*       mom   tq, niter   solves current mom stack for charges  
*			tq = total charge, niter = number of iterations (20 default)
*
*       time  return time of day and elapsed time (not on all machines)
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
*  	active i1 i2; <i2 optional> active atoms i1 to i2 (default is active)
*       inactive i1 i2; < i2 optional> inactivate atoms i1 to i2 
*	inactive tether;  loop throught the tether list and inactivate those atoms
*       nzinactive i1 i2; < i2 optional> inactivate atoms i1 to i2 that
*                               are not 0 0 0  
*
*
* 	grasp nstep nopt imin imax atom;  GRASP in torsion space
*	genetic nstep ndeep sigma target n_opt_steps ; genetic optimizer
*	gsdg  niter min_atom max_atom; iterative distance geometry bounded by
*                                       serial numbers
*       abuild niter min_atom max_atom;  analytic builder using bond,angle,noel
*                                       bounded by serial numbers
*	bell  niter min_atom max_atom; iterative distance geometry bounded by
*                                       serial numbers
*       Kohonen niter radius <1 initialize -1 continue> < r | r z | rx ry rz>
*                         kohonen search
*       Kdock niter radius <1 initialize -1 continue>  rx ry rz
*                         kohonen search on limited region around rx ry rz
*                          niter is the number of steps, not numatom*steps
*
*	dgeom niter origin shift;  standard distance geometry
*                             implemented with the power method
*                              origin is the atom to use as the key
*                              shift is the amount of eigenvalue shift
*
*	normal damp    ;    calculate the normal modes  if damp > 0 output them
*
*	table id n ; create empty sorted table
*       tableent id who r v ; add the who'th element to the table it
*       access with use tbond
*
* direct SCF terms
*       orbit <o1,o1o,o2,o2p,o3,o4s,o4p,om> i1,<i2-i5>,osn, parameters, ipair ;
*               ipair == 2 (doublet) ipair == 1 (singlet)
*       expand osn,n,a,r,a,r (up to 6)  ;
*                          these define an orbital
*
*       dscf <coef,expo,xyz,geom,anal> n toler;  optimize the orbitals
*            <coefficients, exponents, atom center, orbital geometry>;
*
*
*    fdfield Dielectric, Z, T, gaurd, dgrid, radius;
*         finite element field solver  
*
*
*	others like fix,and... TBD
*   first nonblank == '#' is a comment and the line is skipped 
*
*
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
#define MAXTOKEN 20 
#define TOKENLENGTH 80 
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <time.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#ifdef ESV
#define NOTIME
#endif
#ifndef NOTIME
#define TIME
#endif

#include "ammp.h"
#include "orbit.h"
#ifdef GRACELESS
#include "graceless.h"
#endif

void read_eval_do( ip,op )
FILE *ip,*op;
{
    char line[4096], *ap, *ap1; /* buffer and pointers for sscanf */
    int  itemp[MAXTOKEN];  /* integers to read */
    float ftemp[MAXTOKEN]; /* floats to read */
    int inbuff;     /* where am i in the buffer */
    int inliteral ; /* inside of a literal string */
    int eval();
    /* always start outside of a literal */

    inliteral = (1 == 0); /* portable lexpresion */
    /* error checking */
    if( ip == NULL)
    {
        aaerror(" cannot use input file \n");
        return ;
    }
    if( op == NULL)
    {
        aaerror(" cannot use output file \n");
        return;
    }
    /* lex out a line  (i.e.  <stuff> ; )
    * filter out \n and the \; */
    inbuff = 0;
    while((line[inbuff]=fgetc(ip))!= (char)EOF)
    {/* start of the lex while */
        if( !inliteral && line[inbuff] == '"') inliteral = (1==1);
        if(  inliteral && line[inbuff] == '"') inliteral = (1==0);
        if( !inliteral ){
            if( line[inbuff] == ',') line[inbuff] = ' ';
            if( line[inbuff] == '\t') line[inbuff] = ' ';
            if( line[inbuff] == '\n') line[inbuff] = ' ';
            /* ascii assumption */
            if( line[inbuff] <  ' ') line[inbuff] = ' ';
            if( line[inbuff] == ';')
            {
                line[inbuff] = '\0';
                if(eval(ip,op,line )<0) return;
                fflush(op);
                inbuff = 0;
            } else if( line[inbuff] != '\n') inbuff++;
        }/* end of inliteral */
    }/* end of the lex while */
}/* end of routine */

/* eval actually parses the line */
/* original version used sscanf *
*  current version lexes tokens and if numeric
*  converts them to integer and floating point versions
*/ 
int  (*potentials[10])(),(*forces[10])(),nused=(-1);
int eval( ip,op,line )
FILE *ip,*op;
char *line;
{
    FILE *newfile,*fopen(),*tmpfile();
    char  token[MAXTOKEN][TOKENLENGTH],*ap, *ap1;
    /* buffer and pointers for sscanf */
    char   errmes[BUFSIZ];
    int  itemp[MAXTOKEN],itoken,tisvariable(),tisint();  /* integers to read */
    float ftemp[MAXTOKEN]; /* floats to read */
    /*static  int  (*potentials[10])(),(*forces[10])(),nused=(-1); */
    static  int  echo=1;
    static int inloop = 1;
    int get_i_variable(),set_i_variable(),set_f_variable();
    int clone(),statclone(),subclone(),nnclone();
    float get_f_variable();
    int v_bond(),f_bond(),v_angle(),f_angle();
    int f_tbond(),v_tbond();
    int v_mmbond(),f_mmbond(),v_mmangle(),f_mmangle();
    int v_c_angle(), f_c_angle();
    int v_periodic(),f_periodic();
    int v_nonbon(),f_nonbon(),v_torsion(),f_torsion();
    int ttarget(),v_ttarget(),f_ttarget();
    int swarm(),v_swarm(),f_swarm();
    int v_box(),f_box();
    int v_damp(),f_damp();
    int v_cube(),f_cube();
    int v_react(),f_react();
    int v_screen(),f_screen();
    int v_gauss(),f_gauss();
    int u_v_nonbon(),u_f_nonbon();
    int v_hard(),f_hard();
    int v_ho_bond(),f_ho_bond();
    int v_ho_angle(),f_ho_angle();
    int v_trace(),f_trace();
    int v_image(),f_image();
    int v_scf(),f_scf();
    int v_central(),f_central();
    int atom(),bond(),angle(),torsion(),a_readvelocity();
    int morse(),v_morse(),f_morse();
    int restrain(),v_restrain(),f_restrain();
    int tether(),v_tether(),f_tether(),alltether(), inactive_tether();
    int v_ho_tether(),f_ho_tether();
    int hybrid(),v_hybrid(),f_hybrid();
    int av5(),v_av5(),f_av5();
    int v_ho_av5(),f_ho_av5();
    int abc(),v_abc(),f_abc();
    int noel(),v_noel(),f_noel();
    int v_noe_lin(), f_noe_lin();
    int dill(),v_dill(),f_dill();
    int step(),v_step(),f_step();
    int v_ho_noel(),f_ho_noel();
    int v_debye(),f_debye();
    int v_shadow(),f_shadow();
    int init_fourd(),v_fourd(),f_fourd();
    int math();
    int tbond(),create_table(),add_pair_to_table();
    int gtor_define(), gtor_gene();
#ifdef TIME
    clock_t clock();
#endif
    void gsdg();
    int a_build();
    void bellman();
    void kohonen();
    void kdock();
    void dgeom();
    int grasp();
    void gene();
    void analyze();
    void monitor();
    void mom_add(),mom();
    void dipole();
    void tailor_qab();
    void tailor_include();
    void tailor_exclude();
    void gdock();
    int fd_field();
    int verlet(),v_maxwell(),v_rescale();
    int pac(),pacpac(),tpac(),hpac(),ppac();
    int ptpac();
    int doubletime();
    int FDnormal();
    /* default setup of potentials and forces */
    if( nused == -1) {
        potentials[0] = v_bond;
        potentials[1] = v_angle;
        potentials[2] = u_v_nonbon;
        potentials[3] = v_torsion;
        potentials[4] = v_hybrid;
        forces[0] = f_bond;
        forces[1] = f_angle;
        forces[2] = u_f_nonbon;
        forces[3] = f_torsion;
        forces[4] = f_hybrid;
        nused = 5;
    }
    /* for safety and to avoid side effects the token arrays are zero'd */
    for( itoken=0; itoken<MAXTOKEN; itoken++)
    {
        token[itoken][0] = '\0';
        itemp[itoken] = 0; ftemp[itoken] = 0.;
    }
    /* now extract tokens and prepare to match it */

    if( echo ) fprintf(op,"%s;\n",line);
    ap = line;
    for( itoken=0; itoken< MAXTOKEN; itoken++)
    {
        ap1 = &token[itoken][0];
        *ap1 = '\0';
        while(*ap == ' ') ap++;
        if( *ap == '"') { /* its a literal copy until '"' is seen */
            ap++;
            while( *ap != '"' && *ap != '\0')
            {
                *(ap1++) = *(ap++);
            }
            if( *ap == '"' ) ap++;
        }else {/* not literal */
            if( itoken== 0 && *ap == '#') return 1;
            while(*ap != ' ' && *ap != '\0')
            {
                if( itoken == 0 ||  ( strcmp(&token[0][0],"read") != 0 &&
                                      strcmp(&token[0][0],"output") != 0)  )
                {
                    if( isupper(*ap))
                    {*ap1 = tolower(*ap); }else{ *ap1 = *ap;}
                    ap1++; ap++;
                } else {
                    *ap1 = *ap; ap1++; ap++;
                }
            }
        }/* end of if not literal */
        *ap1 = '\0';
        /*  now have a list of lexed tokens */
        /*	printf(" %d %s \n",itoken,&token[itoken][0]);  */
        /* if the token is a number atof or atoi it */
        ap1 = &token[itoken][0];
        if( tisvariable(ap1) )
        {
            /* printf(" %s is a variable\n",ap1);  */
            /*  here is where the variable fetch command will go */
            ftemp[itoken] = get_f_variable( ap1);
            itemp[itoken] = get_i_variable( ap1);
        } else{
            if( tisint(ap1) == 1  )
            {
                /*		printf(" %s is an integer\n",ap1); */
                itemp[itoken] = atoi(ap1);
                ftemp[itoken] = itemp[itoken];
            }else{
                /*		printf(" %s is a float\n",ap1);   */
                ftemp[itoken] = atof(ap1);
                /*		sscanf( ap1,"%g",&ftemp[itoken]); */
                itemp[itoken] = (int)ftemp[itoken];
            }
        }
        if(*ap == '\0') break;
    }
    if( token[0][0]  == '\0') return 1 ;
    /* blank lines are not an error */

    /* the block ifs are used rather than a switch to manage
    * potential complexity of  man commands 
    * each if done in the general pattern will have no  side
    * effects and is therfor well phrased.
    * needless to say common commands should be first 
    */ 
    if( strcmp( &token[0][0], "atom" ) == 0 )
    {
        if(atom( ftemp[1],ftemp[2],ftemp[3],itemp[4],ftemp[6],
                 ftemp[7],ftemp[8],ftemp[9],&token[5][0]) )
        { } else {
            aaerror(" cannot add to atom structure -data structure error");
            exit(0); }
        goto DONE;
    }

    if( strcmp( &token[0][0], "bond" ) == 0 )
    {
        if( bond(itemp[1],itemp[2],ftemp[3],ftemp[4]))
        { }else
        {
            aaerror(" cannot add to bond structure -data structure error");
            exit(0); }
        goto DONE;
    }

    if( strcmp( &token[0][0], "restrain" ) == 0 )
    {
        if( restrain(itemp[1],itemp[2],ftemp[3],ftemp[4]))
        { }else
        { aaerror(" cannot add to restrain structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "morse" ) == 0 )
    {
        if( morse(itemp[1],itemp[2],ftemp[3],ftemp[4],ftemp[5]))
        { }else
        { aaerror(" cannot add to morse structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "angle" ) == 0 )
    {
        ftemp[2] = 3.141592653589793/180.;
        ftemp[5] = ftemp[5]*ftemp[2];
        if( angle( itemp[1],itemp[2],itemp[3],ftemp[4],ftemp[5]) )
        { } else
        { aaerror(" cannot add to angle structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0],"central") == 0)
    {
        central(itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "step" ) == 0)
    {
        step( itemp[1],itemp[2],ftemp[3],ftemp[4],ftemp[5],ftemp[6],ftemp[7],ftemp[8]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "dill" ) == 0)
    {
        dill( itemp[1],itemp[2],ftemp[3],ftemp[4],ftemp[5]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "noel" ) == 0 )
    {
        if( noel(itemp[1],itemp[2],ftemp[3],ftemp[4],ftemp[5],
                 ftemp[6],ftemp[7]))
        { }else
        {
            aaerror(" cannot add to noel structure -data structure error");
            exit(0); }
        goto DONE;
    }

    if( strcmp (&token[0][0],"abc") == 0 )
    {
        ftemp[4] *= 3.14159265/180. ;
        ftemp[5] *= 3.14159265/180. ;
        if(
            abc(itemp[1],itemp[2],itemp[3],ftemp[4],ftemp[5],ftemp[6],ftemp[7],ftemp[8]) )
        { } else
        { aaerror(" cannot add to abc structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp (&token[0][0],"av5") == 0 )
    {
        if(
            av5(itemp[1],itemp[2],itemp[3],itemp[4],itemp[5],ftemp[6],ftemp[7] ))
        { } else
        { aaerror(" cannot add to av5 structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "ttarget" ) == 0 )
    {
        if(
            ttarget(itemp[1],itemp[2],itemp[3],itemp[4],3.14159265*ftemp[5]/180.,ftemp[6]))
        { } else{ aaerror(" cannot add to a ttarget structure"); }
        goto DONE;
    }

    if( strcmp( &token[0][0], "torsion" ) == 0 )
    {
        ftemp[2] = acos(-1.)/180.; ftemp[7] = ftemp[7]*ftemp[2];
        if( torsion(itemp[1],itemp[2],itemp[3],itemp[4],ftemp[5],itemp[6],
                    ftemp[7]) )
        { } else
        { aaerror(" cannot add to torsion structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "hybrid" ) == 0 )
    {
        if(hybrid(itemp[1],itemp[2],itemp[3],itemp[4],ftemp[5],ftemp[6]))
        { } else
        { aaerror(" cannot add to hybrid structure -data structure error");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "swarm") == 0)
    {
        if( itemp[2] > 8) itemp[2] = 8;
        swarm( ftemp[1],itemp[2],itemp[3],itemp[4],itemp[5],itemp[6],
               itemp[7],itemp[8],itemp[9],itemp[10]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "table") == 0)
    {
        create_table( itemp[1],itemp[2]);
        /*
        old format requires huge token array
        		for( itemp[0]=0; itemp[0]< itemp[2]; itemp[0]++)
        		{
        		add_pair_to_table( itemp[1], itemp[0], ftemp[3+itemp[0]*2],
        				ftemp[4+itemp[0]*2]);
        		}
        */
        goto DONE;
    }
    if( strcmp( &token[0][0],"tableent") == 0)
    {
        add_pair_to_table( itemp[1],itemp[2],ftemp[3],ftemp[4]);
        goto DONE;
    }
    if( strcmp( &token[0][0] ,"tbond") == 0)
    {
        if( tbond( itemp[1],itemp[2],itemp[3],ftemp[4]))
        {}else{ aaerror(" error in TBOND");}
        goto DONE;
    }

    if( strcmp(&token[0][0],"clone") == 0)
    {
        clone( itemp[1]);
        goto DONE;
    }
    if( strcmp(&token[0][0],"statclone") == 0)
    {
        statclone(op);
        goto DONE;
    }
    if( strcmp(&token[0][0],"subclone") == 0)
    {
        for( itoken=1; itoken< MAXTOKEN; itoken++)
        {
            if( token[itoken][0] == '\0') break;
        }
        subclone(op, &itemp[1], itoken-1);

        goto DONE;
    }
    if( strcmp(&token[0][0],"nnclone") == 0)
    {
        for( itoken=1; itoken< MAXTOKEN; itoken++)
        {
            if( token[itoken][0] == '\0') break;
        }
        nnclone(op, &itemp[1], itoken-1,4.,8.);

        goto DONE;
    }




    if( strcmp( &token[0][0], "orbit" ) == 0 )
    {
        /* the orbital types (Or1 ... ) are defined in orbit.h */
        if( strcmp(&token[1][0], "1")== 0){
            orbital(Or1,itemp[2],-1,-1,-1,-1,-1,itemp[3],0.,0.,0.,0., itemp[4],itemp[5]);        goto DONE;
        }
        if( strcmp(&token[1][0], "1o")== 0){
            orbital(Or1o,itemp[2],-1,-1,-1,-1,-1,
                    itemp[3],ftemp[4],ftemp[5],ftemp[6],ftemp[7], itemp[8],itemp[9]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "2")== 0){
            orbital(Or2,itemp[2],itemp[3],-1,-1,-1,-1,
                    itemp[4],ftemp[5],0.,0.,0., itemp[6],itemp[7]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "2p")== 0){
            orbital(Or2p,itemp[2],itemp[3],-1,-1,-1,-1,
                    itemp[4],ftemp[5],0.,0.,0., itemp[6],itemp[7]);
            goto DONE;
        }

        if( strcmp(&token[1][0], "3")== 0){
            orbital(Or3,itemp[2],itemp[3],itemp[4],-1,-1,-1,
                    itemp[5],ftemp[6],ftemp[7],0.,0., itemp[8],itemp[9]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "4s")== 0){
            orbital(Or4s,itemp[2],itemp[3],itemp[4],itemp[5],-1,
                    -1,itemp[6],ftemp[7],0.,0.,0., itemp[8],itemp[9]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "4p")== 0){
            orbital(Or4p,itemp[2],itemp[3],itemp[4],itemp[5],-1,
                    -1,itemp[6],ftemp[7],0.,0.,0., itemp[8],itemp[9]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "m")== 0){
            orbital(Orm,itemp[2],itemp[3],itemp[4],itemp[5],
                    itemp[6],itemp[7],itemp[8],0.,0.,0.,0., itemp[9],itemp[10]);
            goto DONE;
        }
        goto DONE;
    }


    if( strcmp( &token[0][0], "expand" ) == 0 )
    {
        expand(itemp[1],itemp[2],
               ftemp[3],ftemp[4], ftemp[5],ftemp[6], ftemp[7],ftemp[8],
               ftemp[9],ftemp[10], ftemp[11],ftemp[12], ftemp[13],ftemp[14]);
        goto DONE;
    }



    if( strcmp( &token[0][0], "velocity" ) == 0 )
    {
        if( a_readvelocity(itemp[1],ftemp[2],ftemp[3],ftemp[4]) )
        { } else
        { aaerror(" cannot update velocity -is this atom defined? ");
            exit(0);}
        goto DONE;
    }
    if( strcmp( &token[0][0], "tether" ) == 0 )
    {
        if( strcmp( &token[1][0], "all") == 0 )
        {
            if( alltether( ftemp[2] ) )
            {} else
            { aaerror(" cannot add to tether structure -data structure error");
                exit(0);}
        }else
        {
            if( tether( itemp[1],ftemp[2],ftemp[3],ftemp[4],ftemp[5]) )
            { } else
            { aaerror(" cannot add to tether structure -data structure error");
                exit(0);}
        }
        goto DONE;
    }

    if( strcmp( &token[0][0],"gtgroup") == 0)
    {
        gtor_define( itemp[1],itemp[2],itemp[3],itemp[4],ftemp[5],itemp[6]);
        goto DONE ;
    }
    if( strcmp( &token[0][0],"tgroup") == 0)
    {
        tgroup( itemp[1],itemp[2],itemp[3],itemp[4],itemp[5],ftemp[6],itemp[7]);
        goto DONE ;
    }

    if( strcmp( &token[0][0],"tsearch") == 0)
    {
        tsearch( itemp[1],itemp[2],itemp[3],itemp[4],itemp[5],itemp[6],itemp[7],itemp[8]);
        goto DONE ;
    }
    if( strcmp( &token[0][0],"tset") == 0 )
    {
        tset( op,echo,itemp[1],itemp[2],itemp[3],itemp[4],ftemp[5]*3.141592653589793/180.);
        goto DONE ;
    }
    if( strcmp( &token[0][0],"tmin") == 0 )
    {
        tmin( op,echo,itemp[1],itemp[2],itemp[3],itemp[4],itemp[5],potentials,nused);
        goto DONE ;
    }
    if( strcmp( &token[0][0],"tmap") == 0 )
    {
        tmap( op,echo,potentials,nused,
              itemp[1],itemp[2],itemp[3],itemp[4],
              itemp[5],itemp[6],itemp[7],itemp[8],itemp[9],itemp[10]
            );
        goto DONE ;
    }
    if( strcmp( &token[0][0] ,"rdebye") == 0)
    {
        debye_rdebye(itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp(&token[0][0],"normal") == 0)
    {
        FDnormal( forces,nused,echo,op,ftemp[1] );
        goto DONE;
    }

    if( strcmp( &token[0][0], "mompar" )== 0)
    {
        mom_param( itemp[1],ftemp[2],ftemp[3] );
        goto DONE;
    }
    if( strcmp( &token[0][0], "momadd" )== 0)
    {
        mom_add( itemp[1],itemp[2] );
        goto DONE;
    }
    if( strcmp( &token[0][0], "mom" )== 0)
    {
        mom( op,ftemp[1],itemp[2] );
        goto DONE;
    }
    if( strcmp( &token[0][0], "monitor" )== 0)
    {
        monitor( potentials,forces,nused,op );
        goto DONE;
    }
    if( strcmp( &token[0][0], "nzinactive" ) == 0)
    {
        inactivate_non_zero( itemp[1],itemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "inactive" ) == 0 )
    {
        if( strcmp( &token[1][0] ,"tether") != 0 )
            inactivate( itemp[1],itemp[2]);
        else
            inactive_tether();
        goto DONE;
    }
    if( strcmp( &token[0][0], "active" ) == 0)
    {
        activate( itemp[1],itemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "init4d" ) == 0)
    {
        init_fourd(ftemp[1] );
        goto DONE;
    }
    if( strcmp( &token[0][0], "analyze" )== 0)
    {
        analyze( potentials,nused,itemp[1],itemp[2],op );
        goto DONE;
    }
    if( strcmp( &token[0][0], "dipole" ) == 0)
    {
        dipole( op, itemp[1],itemp[2]);
        goto DONE;
    }
    if( strcmp(&token[0][0] , "tailor" ) == 0)
    {
        if( strcmp(&token[1][0], "qab") == 0 )
        {
            tailor_qab( itemp[2], ftemp[3],ftemp[4],ftemp[5]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "include") == 0 )
        {
            tailor_include( itemp[2],itemp[3]);
            goto DONE;
        }
        if( strcmp(&token[1][0], "exclude") == 0 )
        {
            tailor_exclude( itemp[2],itemp[3]);
            goto DONE;
        }
        aaerror(" undefined tailor option "); aaerror(&token[1][0]);
        goto DONE;
    }

    if( strcmp( &token[0][0], "read" ) == 0 )
    {
        newfile = fopen( &token[1][0],"r");
        if( newfile == NULL )
        { aaerror(" cannot open file for read "); aaerror(&token[1][0]); }
        else
        { read_eval_do(newfile,op); fclose(newfile); }
        goto DONE;
    }
    if( strcmp( &token[0][0], "output" ) == 0 )
    {
        /* if a non-zero version then write it out */
        if( itemp[2] > 0)
        {
            sprintf( errmes,"%s.%d",&token[1][0],itemp[2]);
            newfile = fopen( errmes,"w");
        } else {
            newfile = fopen( &token[1][0],"w");
        }
        if( newfile == NULL )
        { aaerror(" cannot open file for write "); aaerror(&token[1][0]); }
        else
        { read_eval_do(ip,newfile); }
        goto DONE;
    }
    if( strcmp( &token[0][0], "dump" ) == 0 )
    {
        for( itoken=1; itoken<MAXTOKEN; itoken++)
        {if( token[itoken][0] == '\0') goto DONE;
            if( strcmp(&token[itoken][0],"atom") == 0) dump_atoms( op);
            if( strcmp(&token[itoken][0],"bond") == 0) dump_bonds( op);
            if( strcmp(&token[itoken][0],"central") == 0) dump_central( op);
            if( strcmp(&token[itoken][0],"dill") == 0) dump_dills( op);
            if( strcmp(&token[itoken][0],"step") == 0) dump_steps( op);
            if( strcmp(&token[itoken][0],"noel") == 0) dump_noels( op);
            if( strcmp(&token[itoken][0],"abc") == 0) dump_abcs( op);
            if( strcmp(&token[itoken][0],"av5") == 0) dump_av5s( op);
            if( strcmp(&token[itoken][0],"morse") == 0) dump_morse( op);
            if( strcmp(&token[itoken][0],"angle") == 0) dump_angles( op);
            if( strcmp(&token[itoken][0],"torsion") == 0) dump_torsions( op);
            if( strcmp(&token[itoken][0],"ttarget") == 0) dump_ttargets( op);
            if( strcmp(&token[itoken][0],"hybrid") == 0) dump_hybrids( op);
            if( strcmp(&token[itoken][0],"restrain") == 0) dump_restrains( op);
            if( strcmp(&token[itoken][0],"orbit") == 0) dump_orbit( op);
            if( strcmp(&token[itoken][0],"pdb") == 0) dump_pdb( op,100);
            if( strcmp(&token[itoken][0],"variable") == 0) dump_variable(op);
            if( strcmp(&token[itoken][0],"velocity") == 0) dump_velocity(op);
            if( strcmp(&token[itoken][0],"force") == 0) dump_force(op);
            if( strcmp(&token[itoken][0],"tether") == 0) dump_tethers(op);
            if( strcmp(&token[itoken][0],"tgroup") == 0) dump_tgroup(op);
            if( strcmp(&token[itoken][0],"gtgroup") == 0) dump_gtor(op);
            if( strcmp(&token[itoken][0],"table") == 0) dump_table(op);
            if( strcmp(&token[itoken][0],"tbond") == 0) dump_tbond(op);
            if( strcmp(&token[itoken][0],"swarm") == 0) dump_swarms(op);
            if( strcmp(&token[itoken][0],"ghost") == 0) dump_ghost(op);
        }
        goto DONE;
    }
    if( strcmp( &token[0][0], "use" ) == 0 )
    {
        for( itoken=1; itoken<MAXTOKEN; itoken++)
        {
            if( token[itoken][0] == '\0') goto DONE;
            if( strcmp(&token[itoken][0],"none") == 0) nused = 0;
            if( strcmp(&token[itoken][0],"nonbon") == 0)
            {forces[nused] = u_f_nonbon; potentials[nused++] = u_v_nonbon;}
            if( strcmp(&token[itoken][0],"debye") == 0)
            {forces[nused] = f_debye; potentials[nused++] = v_debye;}
            if( strcmp(&token[itoken][0],"shadow") == 0)
            {forces[nused] = f_shadow; potentials[nused++] = v_shadow;}
            if( strcmp(&token[itoken][0],"fourd") == 0)
            {forces[nused] = f_fourd; potentials[nused++] = v_fourd;}
            if( strcmp(&token[itoken][0],"react") == 0)
            {forces[nused] = f_react; potentials[nused++] = v_react;}
            if( strcmp(&token[itoken][0],"screen") == 0)
            {forces[nused] = f_screen; potentials[nused++] = v_screen;}
            if( strcmp(&token[itoken][0],"gauss") == 0)
            {forces[nused] = f_gauss; potentials[nused++] = v_gauss;}
            if( strcmp(&token[itoken][0],"periodic") == 0)
            {forces[nused] = f_periodic; potentials[nused++] = v_periodic;}
            if( strcmp(&token[itoken][0],"bond") == 0)
            {forces[nused] = f_bond; potentials[nused++] = v_bond;}
            if( strcmp(&token[itoken][0],"mmbond") == 0)
            {forces[nused] = f_mmbond; potentials[nused++] = v_mmbond;}
            if( strcmp(&token[itoken][0],"hobond") == 0)
            {forces[nused] = f_ho_bond; potentials[nused++] = v_ho_bond;}
            if( strcmp(&token[itoken][0],"tether") == 0)
            {forces[nused] = f_tether; potentials[nused++] = v_tether;}
            if( strcmp(&token[itoken][0],"hotether") == 0)
            {forces[nused] = f_ho_tether; potentials[nused++] = v_ho_tether;}
            if( strcmp(&token[itoken][0],"restrain") == 0)
            {forces[nused] = f_restrain; potentials[nused++] = v_restrain;}
            if( strcmp(&token[itoken][0],"morse") == 0)
            {forces[nused] = f_morse; potentials[nused++] = v_morse;}
            if( strcmp(&token[itoken][0],"angle") == 0)
            {forces[nused] = f_angle; potentials[nused++] = v_angle;}
            if( strcmp(&token[itoken][0],"abc") == 0)
            {forces[nused] = f_abc; potentials[nused++] = v_abc;}
            if( strcmp(&token[itoken][0],"av5") == 0)
            {forces[nused] = f_av5; potentials[nused++] = v_av5;}
            if( strcmp(&token[itoken][0],"hoav5") == 0)
            {forces[nused] = f_ho_av5; potentials[nused++] = v_ho_av5;}
            if( strcmp(&token[itoken][0],"hoangle") == 0)
            {forces[nused] = f_ho_angle; potentials[nused++] = v_ho_angle;}
            if( strcmp(&token[itoken][0],"mmangle") == 0)
            {forces[nused] = f_mmangle; potentials[nused++] = v_mmangle;}
            if( strcmp(&token[itoken][0],"cangle") == 0)
            {forces[nused] = f_c_angle; potentials[nused++] = v_c_angle;}
            if( strcmp(&token[itoken][0],"torsion") == 0)
            {forces[nused] = f_torsion; potentials[nused++] = v_torsion;}
            if( strcmp(&token[itoken][0],"ttarget") == 0)
            {forces[nused] = f_ttarget; potentials[nused++] = v_ttarget;}
            if( strcmp(&token[itoken][0],"hybrid") == 0)
            {forces[nused] = f_hybrid; potentials[nused++] = v_hybrid;}
            if( strcmp(&token[itoken][0],"hard") == 0)
            {forces[nused] = f_hard; potentials[nused++] = v_hard;}
            if( strcmp(&token[itoken][0],"honoel") == 0)
            {forces[nused] = f_ho_noel; potentials[nused++] = v_ho_noel;}
            if( strcmp(&token[itoken][0],"noe_lin") == 0)
            {forces[nused] = f_noe_lin; potentials[nused++] = v_noe_lin;}
            if( strcmp(&token[itoken][0],"noel") == 0)
            {forces[nused] = f_noel; potentials[nused++] = v_noel;}
            if( strcmp(&token[itoken][0],"trace") == 0)
            {forces[nused] = f_trace; potentials[nused++] = v_trace;}
            if( strcmp(&token[itoken][0],"tbond") == 0)
            {forces[nused] = f_tbond; potentials[nused++] = v_tbond;}
            if( strcmp(&token[itoken][0],"scf") == 0)
            {forces[nused] = f_scf; potentials[nused++] = v_scf;}
            if( strcmp(&token[itoken][0],"image") == 0)
            {forces[nused] = f_image; potentials[nused++] = v_image;}
            if( strcmp(&token[itoken][0],"box") == 0)
            {forces[nused] = f_box; potentials[nused++] = v_box;}
            if( strcmp(&token[itoken][0],"damp") == 0)
            {forces[nused] = f_damp; potentials[nused++] = v_damp;}
            if( strcmp(&token[itoken][0],"cube") == 0)
            {forces[nused] = f_cube; potentials[nused++] = v_cube;}
            if( strcmp(&token[itoken][0],"dill") == 0)
            {forces[nused] = f_dill; potentials[nused++] = v_dill;}
            if( strcmp(&token[itoken][0],"step") == 0)
            {forces[nused] = f_step; potentials[nused++] = v_step;}
            if( strcmp(&token[itoken][0],"swarm") == 0)
            {forces[nused] = f_swarm; potentials[nused++] = v_swarm;}
            if( strcmp(&token[itoken][0],"central") == 0)
            {forces[nused] = f_central; potentials[nused++] = v_central;}
        }
        goto DONE;
    }
    if( strcmp( &token[0][0], "close" ) == 0 )
    {
        if( op != stdout )
        {
            fclose(op);
            return -1;
        }   goto DONE;
    }
    if( strcmp( &token[0][0], "seti" ) == 0 )
    {
        if( token[1][0] == '\0')
        {aaerror("seti requires a variable name: seti <name> value");
            goto DONE;
        }
        set_i_variable( &token[1][0], itemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "setf" ) == 0 )
    {
        if( token[1][0] == '\0')
        {aaerror("setf requires a variable name: setf <name> value");
            goto DONE;
        }
        set_f_variable( &token[1][0], ftemp[2]);
        goto DONE;
    }
    if( math( token,ftemp,itemp,ip,op,echo ) > 0 ) goto DONE;
    if( strcmp( &token[0][0], "trunc" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        tncnewt( potentials,forces,nused,itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "v_maxwell") == 0)
    {
        v_maxwell( ftemp[1],ftemp[2],ftemp[3],ftemp[4]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "v_rescale") == 0)
    {
        v_rescale( ftemp[1]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "verlet") == 0)
    {
        verlet( forces,nused, itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "pac") == 0)
    {
        pac( forces,nused, itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "tpac") == 0)
    {
        tpac( forces,nused, itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "ppac") == 0)
    {
        ppac( forces,nused, itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "ptpac") == 0)
    {
        ptpac( forces,nused, itemp[1],ftemp[2],ftemp[3],ftemp[4]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "hpac") == 0)
    {
        hpac( forces,potentials,nused, itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "pacpac") == 0)
    {
        pacpac( forces,nused, itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0],"richard") == 0)
    {
        richard( op,echo, potentials,forces,nused, &token[1][0],itemp[2],
                 ftemp[3],ftemp[4],ftemp[5]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "doubletime") == 0)
    {
        doubletime( forces,nused, itemp[1],ftemp[2],ftemp[3],ftemp[4]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "trust" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        trust( potentials,forces,nused,itemp[1],ftemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "steep" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        steep( potentials,forces,nused,itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "bfgs" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        bfgs( potentials,forces,nused,itemp[1],ftemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "dgeom" ) == 0)
    {
        if( nused <= 0) goto DONE;
        dgeom(echo,op, potentials,nused,itemp[1],itemp[2],ftemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "noegen") == 0)
    {
        noel_generate(ftemp[1],ftemp[2]);
        goto DONE ;
    }

    if( strcmp( &token[0][0], "bell" ) == 0)
    {
        if( nused <= 0) goto DONE;
        bellman( potentials,nused,itemp[1],itemp[2],itemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0],"kohonen") == 0)
    {
        if( nused <= 0) goto DONE;
        kohonen(potentials,forces,nused,
                itemp[1],ftemp[2],itemp[3],ftemp[4],ftemp[5],ftemp[6]);
        goto DONE;
    }
    if( strcmp( &token[0][0],"kdock") == 0)
    {
        if( nused <= 0) goto DONE;
        kdock(potentials,forces,nused,
              itemp[1],ftemp[2],itemp[3],ftemp[4],ftemp[5],ftemp[6]);
        goto DONE;
    }

    if( strcmp( &token[0][0], "gsdg" ) == 0)
    {
        if( nused <= 0) goto DONE;
        gsdg( potentials,nused,itemp[1],itemp[2],itemp[3]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "abuild" ) == 0)
    {
        if( nused <= 0) goto DONE;
        a_build( nused,potentials,forces,itemp[1],itemp[2],itemp[3],op,echo);
        goto DONE;
    }
    if( strcmp( &token[0][0], "dscf") == 0)
    {
        direct_scf(op,itemp[2],ftemp[3],&token[1][0]);
        goto DONE;
    }
    if( strcmp( &token[0][0],"grasp") == 0)
    {
        if( nused <=0) goto DONE;
        grasp(op,echo, potentials,forces,nused,itemp[1],0,itemp[2],itemp[3],
              itemp[4],&token[5][0]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "gtor" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        gtor_gene(op,echo, potentials,nused,itemp[1],itemp[2],itemp[3],
                  itemp[4],ftemp[5]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "genetic" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        gene(op, potentials,forces,nused,itemp[1],itemp[2],ftemp[3],
             ftemp[4],itemp[5]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "cngdel" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        cngdel( potentials,forces,nused,itemp[1],itemp[2],ftemp[3],echo);
        goto DONE;
    }
    /*      polytope imin,imax,niter,vstart,vfinal  polytope of a range
    */
    if( strcmp( &token[0][0], "polytope" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        simplex(ftemp[5],itemp[3],ftemp[4], potentials,nused,itemp[1],itemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "rigid" ) == 0 )
    {
        if( nused <= 0) goto DONE;
        rigid(ftemp[5],itemp[3],ftemp[4], potentials,nused,itemp[1],itemp[2]);
        goto DONE;
    }
    if( strcmp( &token[0][0] ,"ghost") == 0)
    {
	ghost( itemp[1],itemp[2],itemp[3]);
        goto DONE;
	}
    if( strcmp( &token[0][0],"buster") == 0)
	{
	buster(potentials,nused, itemp[1],op);
        goto DONE;
	}
    if( strcmp( &token[0][0] ,"gdock") == 0)
    {
        if( nused <= 0) goto DONE;
        gdock( ftemp[1],itemp[2],itemp[3],ftemp[4],ftemp[5],
               potentials,nused,itemp[6],itemp[7]);
        goto DONE;
    }
    if( strcmp( &token[0][0], "fdfield") == 0)
    {
        fd_field( op,op, ftemp[1],ftemp[2],ftemp[3],ftemp[4],ftemp[5],ftemp[6]);
        goto DONE;
    }
#ifdef TIME 
    if( strcmp( &token[0][0],"time") == 0)
    {
        fprintf( op," %f CPU \n",((float)clock())/CLOCKS_PER_SEC);
        goto DONE;
    }
#endif
    if( strcmp( &token[0][0],"echo" ) == 0)
    {
        echo = 1;
        if( strcmp( &token[1][0],"off") == 0) echo = 0;
        goto DONE;
    }
    if( strcmp( &token[0][0],"exit") ==0 ) exit(0);
    /* looping stuff */
    if( strcmp( &token[0][0],"loopi")  == 0)
    {
        if( token[1][0] == '\0')
        { aaerror(" must have a label to loop to "); goto DONE;}
        if( itemp[4] == 0) itemp[4] = 1;  /* must loop  */
        newfile = tmpfile();
        if( newfile == NULL )
        { aaerror(" cannot open temporary file in loopi"); goto DONE; }
        /* scan the input data until the label is found */
        loadloop( ip,newfile, &token[1][0]);
        /*  now do the loop */
        if( itemp[4] > 0)
        {
            for( itemp[0] = itemp[2];itemp[0]< itemp[3];itemp[0]+=itemp[4])
            {
                inloop = -1;
                if( tisvariable(&token[2][0]))
                    set_i_variable( &token[2][0], itemp[0]);
                rewind( newfile );
                read_eval_do(newfile,op);
            }
        } else{
            for( itemp[0] = itemp[2];itemp[0]< itemp[3];itemp[0]+=itemp[4])
            {
                inloop = -1;
                if( tisvariable(&token[2][0]))
                    set_i_variable( &token[2][0], itemp[0]);
                rewind( newfile );
                read_eval_do(newfile,op);
            }
        }
        inloop = 1;
        fclose( newfile);
        goto DONE;
    }
    if( strcmp( &token[0][0],"loopf")  == 0)
    {
        if( token[1][0] == '\0')
        { aaerror(" must have a label to loop to "); goto DONE;}
        if( ftemp[4] == 0.) ftemp[4] = 1.;  /* must loop  */
        newfile = tmpfile();
        if( newfile == NULL )
        { aaerror(" cannot open temporary file in loopi"); goto DONE; }
        /* scan the input data until the label is found */
        loadloop( ip,newfile, &token[1][0]);
        /*  now do the loop */
        if( ftemp[4] > 0.)
        {
            for( ftemp[0] = ftemp[2];ftemp[0]< ftemp[3];ftemp[0]+=ftemp[4])
            {
                inloop = -1;
                if( tisvariable(&token[2][0]))
                    set_f_variable( &token[2][0], ftemp[0]);
                rewind( newfile );
                read_eval_do(newfile,op);
            }
        } else  {
            for( ftemp[0] = ftemp[2];ftemp[0]> ftemp[3];ftemp[0]+=ftemp[4])
            {
                inloop = -1;
                if( tisvariable(&token[2][0]))
                    set_f_variable( &token[2][0], ftemp[0]);
                rewind( newfile );
                read_eval_do(newfile,op);
            }
        }
        inloop = 1;
        goto DONE;
    }
    /* check if its a label  and return */
    /* inloop returns -1 if in a loop which causes read_eval_do() to
    *  return and activates the loop routine */
    for( itemp[0]=0; itemp[0] < TOKENLENGTH; itemp[0]++)
    {
        if( token[0][itemp[0]] == '\0' || token[0][itemp[0]] == ' ')
        {
            if( itemp[0] == 0) break;
            if( token[0][itemp[0]-1] == ':') return inloop;
        }
    }
    /*  default unrecognized token */
    sprintf(&errmes[0]," unrecognized token >%s<",&token[0][0]);
    aaerror( errmes );
DONE: ;
    return 1;
}
/* aaerror is a general error call function */
aaerror( line )
char *line;
{
    fprintf(stderr ,"%s \n",line);
    return 0;
}
/* function tisvariable( char *p )
*
* returns 1 if the character string contains anything other than
* <+-><0-9>.<0-9><e><+-><0-9> 
* works on a tolowered string !!
*/
int tisvariable( p )
char *p;
{
    if( (*p != '+')&&(*p != '-')&& !(isdigit( (int) *p)) &&(*p != '.') )
        return 1;
    /* now for the rest we check until either '\0' or not a digit */
    p++;
    while( (*p != '\0') && (isdigit( (int) *p) ) ) p++;
    if( *p == '\0') return 0;
    if( (*p != '.') && (*p != 'e') ) return 1;
    p++;
    if( !(isdigit( (int) *p)) ){
        if( *p == '\0' ) return 0;
        if( (*p != '.') && (*p != 'e') ) return 1;
        p++;
    }
    if( *p == '\0') return 0;
    if( (*p != '+')&&(*p != '-')&& !(isdigit( (int) *p)) &&(*p != '.') )
        return 1;
    p++;
    if( *p == '\0') return 0;
    while( (*p != '\0') && ((isdigit( (int) *p))||(*p=='.')) ) p++;
    if( *p == '\0') return 0;
    return 1;
}
/* function tisint( char *p )
*
* check that a string is <+-><0-9>
* return 1 if true
* return 0 if not
*/
int tisint( p)
char *p ;
{
    char *pp;
    pp = p;
    while( *pp != '\0')
    { if( *pp == '.') return 0; pp++;}
    if( (*p != '+')&&(*p != '-')&& !(isdigit( (int) *p)) ) return 0;
    p++;
    while (*p != '\0')
    {
        if( !(isdigit( (int) *p )) ) return 0;
        p++;
    }
    return 1;
}

/* routine loadloop( FILE *ip, FILE *tp, char *label)
*
* read lines from ip and write to tp
* when the line begins with label  stop (after writing it )
*/
int loadloop( ip,tp,label)
FILE *ip,*tp;
char *label;
{
    char line[256], *fgets() ;
    char *sp,*wp;

    /*	printf( " the target label >%s<\n" , label);
    */
    while( fgets(  line,256,ip) != NULL )
    {
        fputs( line,tp );
        fputs("\n",tp);
        sp = line;
        while( *sp == ' ' && *sp != '\0') sp++;
        if( *sp != '\0' )
        {
            wp = sp;
            while(*wp != ';' && *wp != ' ' && *wp != '\0')
            { if( isupper(*wp)){*wp = (char)tolower((int)*wp);}
                wp++;}
            if( *wp == ' ' ) *wp = '\0';
            if( *wp == ';' ) *wp = '\0';
            if( strcmp(sp,label) == 0 ) return;
        }
    }
    aaerror(" must have a label for looping ");
    sprintf(line," where is >%s< label ?\n",label);
    aaerror( line );
    return ;
}
