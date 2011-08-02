/* variable.c
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
#include <stdio.h>
#include <ctype.h>
#include <string.h>
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
VARIABLE *variableFIRST = NULL;
VARIABLE *variableLAST = NULL;
#define variableLONG sizeof(VARIABLE)

/****************************************************************/

/* function match_variable( char *name)
* tries to find the variable in the list returns NULL if not there *
*/
VARIABLE *match_variable( name )
char *name;
{
    VARIABLE *vp;
    vp = variableFIRST;
    while(1 )
    {
        if( vp == NULL ) return vp;
        if( strncmp( &(vp->name[0]), name,NUM_SIG) == 0) return vp;
        if( vp == vp->next) return NULL;
        vp = vp->next;
    }
}
/* function set_f_variable( char *name, float f)
*
*   allocates storage if needed
*   atempts to match name and update, otherwise adds to the list
*/
int set_f_variable( name,f )
char *name;
float f;
{
    VARIABLE *new,*match_variable();
    int i;

    new = match_variable( name);
    if( new == NULL)
    {
        if( ( new = malloc( variableLONG ) ) == NULL)
        {
            return 0;
        }
        new->next = NULL;
        /* initialize the pointers */
        if( variableFIRST == NULL) variableFIRST = new;
        if( variableLAST == NULL) variableLAST = new;
    }
    new->value.f = f;
    new->type = FLOAT;
    for( i=0; i< NUM_TOT; i++)
    {
        new->name[i] = *name;
        if( *name == '\0') break;
        name++;
    }
    if( new->next == NULL)
    {
        new -> next = new;
        variableLAST -> next = new;
        variableLAST = new;
    }
    return 1;
}
/* function set_i_variable( char *name, int f)
*
*   allocates storage if needed
*   atempts to match name and update, otherwise adds to the list
*/
int set_i_variable( name,f )
char *name;
int f;
{
    VARIABLE *new,*match_variable();
    int i;

    new = match_variable( name);
    if( new == NULL)
    {
        if( ( new = malloc( variableLONG ) ) == NULL)
        {
            return 0;
        }
        new->next = NULL;
        /* initialize the pointers */
        if( variableFIRST == NULL) variableFIRST = new;
        if( variableLAST == NULL) variableLAST = new;
    }
    new->value.i = f;
    new->type = INTEGER ;
    for( i=0; i< NUM_TOT; i++)
    {
        new->name[i] = *name;
        if( *name == '\0') break;
        name++;
    }
    if( new->next == NULL)
    {
        new -> next = new;
        variableLAST -> next = new;
        variableLAST = new;
    }
    return 1;
}
/* function get_f_variable(char * name )
*
* returns 0. or the float referenced or converted by name
*/
float get_f_variable( name )
char *name;
{
    VARIABLE *vp;
    VARIABLE *match_variable();

    vp = match_variable(name);
    if( vp == NULL ) return 0.;
    if( vp->type == FLOAT) return vp->value.f ;
    return (float) vp->value.i;
}
/* function get_i_variable(char * name )
*
* returns 0 or the integer referenced or converted by name
*/
int get_i_variable( name )
char *name;
{
    VARIABLE *vp;
    VARIABLE *match_variable();

    vp = match_variable(name);
    /* earlier errorneous version had return 0. */
    if( vp == NULL ) return 0;
    if( vp->type == INTEGER) return vp->value.i ;
    return (int) vp->value.f;
}
/*  dump_variable( FILE *where )
* 
* writes all of the variables out to the file in AMMP syntax 
*/
void dump_variable( where )
FILE *where;
{
    VARIABLE *vp;
    vp = variableFIRST;
    while(1)
    {
        if( vp->next == NULL) return;
        if( vp->type == INTEGER)
        {
            fprintf(where,"seti %s %d \;\n",&vp->name[0],vp->value.i);
        } else
        {
            fprintf(where,"setf %s %f \;\n",&vp->name[0],vp->value.f);
        }
        if( vp->next == vp) return;
        vp = vp->next;
    }
}
