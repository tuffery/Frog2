/* math.c
*  perform basic mathematical operations on data
* return values to memory locations
*  uses floating point contagion
*  defined ops
*   add  a b  a <= a+b
*   sub  a b  a <= a - b
*   mul  a b  a <= a*b
*   div  a b  a <= a/b
*   fix  a    a <= (int) a   (drops the fraction for atomdata)
*   sqrt a    a <= sqrt a  (becomes a float !!!)
*   nop  a    just echo the value
*   mov  a b  a <= b
*   randf a   put a random variable 0-<1. in a
*   max  a b   a <= max( a b )
*   min  a b   b <= min( a b )
*   serial a res atom   move the serial number of the atom
*                       into a  (res, atom are like 100 ca)
*   index a i    move the serial number of the i'th atom into a
*
*   linmin  search the atom.<dx,dy,dz> direction for a minimum
*   je a b label:  jump to label if equal
*   jne a b label:  jump to label if  not equal
*   jl a b label:  jump to label if  a < b
*   jg a b label:  jump to label if  a > b
*   jes a string label: jump to label: if a->label == string  
*   jnes a string label: jump to label: if a->label != string  
*        label:  is a label within the current script file
*                which may be within the current loop
*                the j<elg> commands will rewind the input file
*                and search for label: 
*
*
*   data types
*   imeadiate    e.g. a number like 3.14159
*   variable     e.g. a name  like pi
*   atomdata    serial.<x,y,z,fx,fy,fz,dx,dy,dz,vx,vy,vz,q,a,b,m,chi,jaa.serial>
*		serial may be a variable
*		valid format atoms with non-extant serial numbers will be
*		ignored and the operation will silently not happen
*		this allows sums over discontinous atom ranges
*
*
*  routines defined in this module
*
*	math()  does the work
* 	getatomdata() returns a float * or null for an atomdata
*	validatom()  returns nonzero when atomdata format is valid 
*   
*/
/* header from variable.c
*
* variable storage and retreival routines for AMMP
*
* using scope rules for structuring
*
*
*  variables are stored in linked list form
*
*   get_f_variable( char *name, float *fvalue )
*   get_i_variable( char *name,  int *ivalue )
*	returns variable who matches name (all lower case )
*   set_f_variable( char *name, float fvalue )
*   set_i_variable( char *name, int ivalue )
*	sets variable who matches name (all lower case )
*   match_variable( char *name ) returns pointer to name if there NULL if not
*   dump_variable(FILE *output  )
*	 dumps variables to  file 
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
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h"
#define NUM_SIG 6
#define NUM_TOT 7
#define INTEGER 0
#define FLOAT 1
typedef struct{
    int type;
    char name[NUM_TOT];
    union{  float f; int i;} value;
    void *next;
}  VARIABLE;
/* next define is from ammp.c */
#define TOKENLENGTH 80 

#define ATOMDATA 0
#define IMEADIATE 1
#define VARDATA 2

/* the potentials are defined in eval.c */
extern int (*potentials[10])(),(*forces[10])(),nused;

int math( tokens,fvalue,ivalue,ip,op,echo)
#ifdef ESV
#define tokens (*tokens)
#endif
char tokens[][TOKENLENGTH];
float fvalue[];
int ivalue[];
FILE *ip,*op;
int echo ;
{
    int adata,bdata;
    int atype, btype;
    float *foutpointer,*fp,fa,fb,randf();
    int *ioutpointer,ia,ib,*intp;
    VARIABLE *vp,*vos,*match_variable();
    float *getatomdata(); /* return the address of the atomdata if valid */
    int *getintatomdata(); /* return the address of the atomdata if valid */
    int intatomdata;
    int validatom(); /* return 0 iff not a valid atomdata format */
    int set_f_variable(), set_i_variable();
    int a_number(),tisint();
    ATOM *ap,*a_next(),*a_m_serial();
    int i,j;
    float a_l2_f(),a_max_f(),a_max_d();
    float linmin(); /* defined in optimist.c */
    /* set some default variables */
    set_i_variable( "numatm",a_number());
    set_f_variable( "l2f",a_l2_f());
    set_f_variable( "lmaxf",a_max_f());
    /*
    *   set_f_variable( char *name, float fvalue )
    *   set_i_variable( char *name, int ivalue )
    */
    /* type the variables */
    intatomdata = (1==0);
NEW_AVAR:
    if( (vp = match_variable(&tokens[1][0])) != NULL)
    {
        adata = VARDATA;
        atype = vp->type;
        vos = vp;
        if( atype == FLOAT){ fa = vp->value.f;}
        foutpointer  = &vp->value.f;
        if( atype == INTEGER){ia = vp->value.i; }
        ioutpointer = &vp->value.i;
        goto AFOUND;
    }
    if( (foutpointer = getatomdata(&tokens[1][0])) != NULL)
    {
        adata = ATOMDATA;
        atype = FLOAT;
        fa = *foutpointer;
        goto AFOUND;
    }
    if( (ioutpointer = getintatomdata(&tokens[1][0])) != NULL)
    {
        adata = ATOMDATA;
        atype = INTEGER;
        ia = *ioutpointer;
        intatomdata = (1==1);
        goto AFOUND;
    }
    if( validatom(&tokens[1][0]) != 0 ) return( 1);
    foutpointer = NULL; /* shouldn't change it but this is safe */
    adata = IMEADIATE;
    fa = fvalue[1];
    ia = ivalue[1];
    atype = FLOAT;
    if( tisint( &tokens[1][0]) == 1) atype = INTEGER;
    /*
    	if( tisint( &tokens[1][0]) == 1)printf(" INTEGER "); 
    	printf(">%s<\n",&tokens[1][0]);
    */
AFOUND:
    if( (vp = match_variable(&tokens[2][0])) != NULL)
    {
        /*	bdata = VARDATA;
        */
        btype = vp->type;
        if( btype == FLOAT) fb = vp->value.f;
        if( btype == INTEGER) ib = vp->value.i;
        goto BFOUND;
    }
    if( (fp = getatomdata(&tokens[2][0])) != NULL)
    {
        /*	bdata = ATOMDATA;
        */
        btype = FLOAT;
        fb = *fp;
        goto BFOUND;
    }
    if( (intp = getintatomdata(&tokens[2][0])) != NULL)
    {
        /*	bdata = ATOMDATA;
        */
        btype = INTEGER;
        ib = *intp;
        goto BFOUND;
    }
    if( validatom(&tokens[2][0]) != 0 ) return( 1);
    /*	bdata = IMEADIATE;
    */
    fb = fvalue[2];
    ib = ivalue[2];
    btype = FLOAT;
    if( tisint( &tokens[2][0]) == 1) btype = INTEGER;
    /*
    	if( tisint( &tokens[2][0]) == 1)printf(" INTEGER "); 
    	printf(">%s<\n",&tokens[2][0]);
    */
BFOUND:
    /*  make the converted data */
    if( atype == FLOAT) ia = (int)fa;
    if( atype == INTEGER) fa = (float)ia;
    if( btype == FLOAT) ib = (int)fb;
    if( btype == INTEGER) fb = (float)ib;
    /* is it a defined command ? */
    if( strcmp( &tokens[0][0], "add" ) == 0 )
    {
        if( atype == FLOAT) fa = fa + fb;
        if( atype == INTEGER) ia = ia + ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "sub" ) == 0 )
    {
        if( atype == FLOAT) fa = fa - fb;
        if( atype == INTEGER) ia = ia - ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "mul" ) == 0 )
    {
        if( atype == FLOAT) fa = fa * fb;
        if( atype == INTEGER) ia = ia * ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "div" ) == 0 )
    {
        if( atype == FLOAT) fa = fa / fb;
        if( atype == INTEGER) ia = ia / ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "fix" ) == 0 )
    {
        if( adata == ATOMDATA ) {
            ia = (int)fa;
            *foutpointer = (float)ia;
            if( echo  ) fprintf(op,"%d \n",ia);
            return(1);
        }
    if( atype == FLOAT) {atype = INTEGER; ia = (int)fa;}
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "sqrt" ) == 0 )
    {
        if( adata == ATOMDATA ) {
            if( fa > 0.)
                *foutpointer = sqrt(fa);
            else
                *foutpointer = -sqrt(-fa);
            if( echo  ) fprintf(op,"%f \n",*foutpointer);
            return(1);
        }
        atype = FLOAT;
        if( fa > 0) fa = sqrt(fa);
        else fa = -sqrt(-fa);
        goto GOOD_OP;
    }
    if( strcmp(&tokens[0][0],"linmin") == 0)
    {
        fa = 0.;
        fa = linmin( potentials,nused, sqrt(a_max_d()) );
        if( echo) fprintf(op,"%f step to minimum\n",fa);
        a_inc_d(fa);
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "nop" ) == 0 )
    {
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "mov" ) == 0 )
    {
        ia = ib;
        fa = fb;
        atype = btype;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "max" ) == 0 )
    {
        if( atype == FLOAT && fa < fb) fa = fb;
        if( atype == INTEGER && ia < ib) ia = ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "min" ) == 0 )
    {
        if( atype == FLOAT && fa > fb) fa = fb;
        if( atype == INTEGER && ia > ib) ia = ib;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "randf" ) == 0 )
    {
        atype = FLOAT;
        fa = randf( -1 );
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "serial" ) == 0 )
    {
        ia = 100*ib-1;
        j = ia + 100;
        i = -1;
        while( (ap = a_next(i)) != NULL )
        {
            i = 1;
            if( ap->serial > ia && ap->serial < j)
                if( math_match_atom(&tokens[3][0],ap) != 0 ){
                    atype = INTEGER;
                    ia = ap->serial;
                    goto GOOD_OP;
                }
        }
        ia = -1; /* never a serial number */
        atype = INTEGER;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "index" ) == 0 )
    {
        ap = a_next(-1);
        for( i=0; i< ib; i++)
            ap = a_next(i);
        ia = ap->serial;
        atype = INTEGER;
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "jes" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        ap = a_m_serial( ia ); if( ap == NULL) return 1;
        if( strcmp( &ap->name[0], &tokens[2][0]) == 0)
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "jnes" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        ap = a_m_serial( ia );
        if( ap == NULL) {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
            goto GOOD_OP;
        }
        if( strcmp( &ap->name[0], &tokens[2][0]) != 0)
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "jne" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        if( (atype == INTEGER && ia != ib ) ||
                (atype == FLOAT && fa != fb ) )
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "je" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        if( (atype == INTEGER && ia == ib ) ||
                (atype == FLOAT && fa == fb ) )
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "jg" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        if( (atype == INTEGER && ia > ib ) ||
                (atype == FLOAT && fa > fb ) )
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    if( strcmp( &tokens[0][0], "jl" ) == 0 )
    {
        if( tokens[3][0] == '\0')
        { aaerror("label: required for a jump \n"); goto GOOD_OP;}
        if( (atype == INTEGER && ia < ib ) ||
                (atype == FLOAT && fa < fb ) )
        {
            rewind(ip);
            math_findlabel( ip,&tokens[3][0]);
        }
        goto GOOD_OP;
    }
    return(-1);
    /* if found we jump here */
GOOD_OP:

    if( tisvariable(&tokens[1][0] ) && tokens[1][0] != '\0' && adata == IMEADIATE )
    {set_i_variable(&tokens[1][0],0); adata = VARDATA;
        vos = match_variable(&tokens[1][0]);
        foutpointer = &vos->value.f;
        ioutpointer = &vos->value.i;}

    if( adata != IMEADIATE )
    {
        if( adata == VARDATA) vos->type = atype;
        if( adata == ATOMDATA ){
            if( atype == INTEGER && !intatomdata){atype = FLOAT;}
        }
        if( atype == FLOAT && foutpointer != NULL) *foutpointer = fa;
        if( atype == INTEGER && ioutpointer != NULL) *ioutpointer = ia;
    }
    if( echo && atype == INTEGER ) fprintf(op,"%d \n",ia);
    if( echo && atype == FLOAT ) fprintf(op,"%f \n",fa);
    return( 1);
#ifdef ESV
#undef tokens
#endif
}
/* int validatom()
*  given a string return 0 if it is not of the form
*  stuff.<x y z fx fy fz dx dy dz vx vy vz q a b m chi jaa na serl>
*  these being the valid atomdata parameters 
*  returns 1 to 20 for x y z fx fy fz dx dy dz vx vy vz q a b m chi jaa na ser
*  so we don't have to lex again
*/
int validatom( who )
char *who;
{
    char *cp,*pp,*cp1,*cp2,*cp3;
    int i;
    cp = who;
    i = 0;
    while ( *cp != '\0')
    {
        if( *cp == '.') {i++;  pp = cp;}
        cp++;
    }
    if( i != 1) return (0 ); /* there can only be one '.' */
    cp = pp ; cp++;
    /* now we check the tailing characters */
    cp1 = cp; cp1++;
    cp2 = cp1; cp2++;
    cp3 = cp2; cp3++;
    *cp3 == '\0'; /* truncate the string here */
    if( *cp1 == '\0')
    {
        if( *cp == 'x') return (1);
        if( *cp == 'y') return (2);
        if( *cp == 'z') return (3);
        if( *cp == 'q') return (13);
        if( *cp == 'a') return (14);
        if( *cp == 'b') return (15);
        if( *cp == 'm') return (16);
        return (0);
    }
    if( *cp2 == '\0')
    {
        if( *cp == 'f') {
            if( *cp1 == 'x') return (4);
            if( *cp1 == 'y') return (5);
            if( *cp1 == 'z') return (6);
        }
        if( *cp == 'd' ){
            if( *cp1 == 'x') return (7);
            if( *cp1 == 'y') return (8);
            if( *cp1 == 'z') return (9);
        }
        if( *cp == 'v' ){
            if( *cp1 == 'x') return (10);
            if( *cp1 == 'y') return (11);
            if( *cp1 == 'z') return (12);
        }
        if( *cp == 'n' ) {
            if( *cp1 =='a') return 19;
        }
    }
    if( *cp3 == '\0')
    {
        if( *cp == 'c' && *cp1 == 'h' && *cp2 == 'i') return 17;
        if( *cp == 'j' && *cp1 == 'a' && *cp2 == 'a') return 18;
        if( *cp == 's' && *cp1 == 'e' && *cp2 == 'r') return 20;
    }
    return (0);
}
/* float *getatomdata() returns a float * or null for an atomdata
*/
float *getatomdata( who )
char *who ;
{
    int validatom();
    int i,j;
    char aser[TOKENLENGTH],*cp;
    ATOM *ap,*a_m_serial();
    VARIABLE *vp,*match_variable();
    static float fx,fy,fz;
    i = validatom( who );
    if( i == 0 ) return ( NULL ); /* if not the right format it aint one */

    cp = who; j = 0;
    while( *cp != '.')
    { aser[j++] = *cp; cp++; }
    aser[j] = '\0';
    if( (vp = match_variable(aser)) == NULL)
    { j = atoi(aser); } else {
        if( vp ->type == INTEGER) j = vp->value.i;
        if( vp ->type == FLOAT) j = (int)vp->value.f;
    }
    ap = a_m_serial(j);
    if( ap == NULL) return( NULL );
    if( i == 1 ) return ( &ap->x );
    if( i == 2 ) return ( &ap->y );
    if( i == 3 ) return ( &ap->z );
if( i == 4 ) { fx = ap->fx ; return( &fx);}
    if( i == 5 ) { fy = ap->fy ; return(&fy);}
    if( i == 6 ) { fz = ap->fz ; return(&fz);}
    if( i == 7 ) return ( &ap->dx );
    if( i == 8 ) return ( &ap->dy );
    if( i == 9 ) return ( &ap->dz );
    if( i == 10 ) return ( &ap->vx );
    if( i == 11 ) return ( &ap->vy );
    if( i == 12 ) return ( &ap->vz );
    if( i == 13 ) return ( &ap->q );
    if( i == 14 ) return ( &ap->a );
    if( i == 15 ) return ( &ap->b );
    if( i == 16 ) return ( &ap->mass );
    if( i == 17 ) return ( &ap->chi );
    if( i == 18 ) return ( &ap->jaa );
    if( i == 19 ) return ( &ap->na );
    if( i == 20 ) return ( NULL ); /* 20 is an int */
    return (NULL );
}
int *getintatomdata( who )
char *who ;
{
    int validatom();
    int i,j;
    char aser[TOKENLENGTH],*cp;
    ATOM *ap,*a_m_serial();
    VARIABLE *vp,*match_variable();
    static float fx,fy,fz;
    i = validatom( who );
    if( i == 0 ) return ( NULL ); /* if not the right format it aint one */

    cp = who; j = 0;
    while( *cp != '.')
    { aser[j++] = *cp; cp++; }
    aser[j] = '\0';
    if( (vp = match_variable(aser)) == NULL)
    { j = atoi(aser); } else {
        if( vp ->type == INTEGER) j = vp->value.i;
        if( vp ->type == FLOAT) j = (int)vp->value.f;
    }
    ap = a_m_serial(j);
    if( ap == NULL) return( NULL );
    if( i == 20 ) return ( &ap->serial ); /* 20 is an int */
    return (NULL );
}
/*	math_match_atom(char* atomname, ATOM * ap)
*
*  find if the atomname is part of the atom .name field
*  as in atomname = "ca" and ap->name = arg.ca
*
*/
int math_match_atom( who,ap)
char *who;
ATOM *ap;
{
    char *cp;
    cp = & ap->name[0];
    while( *cp != '.' )
    {
        if( *cp == '\0') return 0;
        cp++;
    }
    cp++;
    /*	printf(" >%s< >%s<\n",cp,who);
    */
    if( strcmp( who,cp ) == 0 )
    { return 1; }
    return 0;
}
/*
*	    math_findlabel(FILE * ip, char *label);
* search a file for label:
*  if label doesn't end in : add it
*
*/
math_findlabel( fp,label)
FILE *fp;
char *label;
{
    char *cp,*lp;
    char llabel[TOKENLENGTH];
    char myline[TOKENLENGTH]; /* since label is no longer than TOKENLENGTH we can skip such lines */
    int  inmyline;
    int i;
    char ac;

    cp = label;
    lp  = &llabel[0];
    while( *cp != '\0')
    {
        *lp = *cp ; lp++; cp ++;
    }
    cp = lp; cp--;
    if( *cp != ':'){ *lp = ':'; lp++;}
    *lp = '\0';
    /* now one char at a time scan ip for a blank */
    inmyline = 0;
    lp  = &llabel[0];
    while( (i= fgetc( fp )) != EOF )
    {
        ac = (char)i;
        if( !isspace((int) ac) )
        {
            if( ac == ';')
            {
                myline[inmyline] = '\0';
                /*	printf(">%s<\n",&myline[0]);
                */
                if( strcmp( lp,&myline[0]) == 0 ) return ;
                inmyline = 0;
            }else{
                if( inmyline > TOKENLENGTH) inmyline = 0;
                myline[inmyline++] = ac;
            }
        }
    }

}
