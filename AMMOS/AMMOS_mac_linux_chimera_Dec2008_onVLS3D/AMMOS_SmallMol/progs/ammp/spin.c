/* spin.c
*
* routines to use floating gaussian spinals for QM calculations
*
*/
/*
*  copyright 1993,1994,1995 Robert W. Harrison
*  
*  This notice may not be removed
*  This program may be copied for scientific use
*  It may not be sold for profit without explicit
*  permission of the author(s) who retain any
*  commercial rights including the right to modify 
*  this notice
*/

#include <assert.h>

#define ANSI 1
/* misc includes - ANSI and some are just to be safe */
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#ifdef ANSI
#include <stdlib.h>
#endif
#include "ammp.h" /* this includes numeric.h */

#include "orbit.h" 


ORBIT *firstORBIT = NULL , *lastORBIT = NULL;
static int spinNUMBER = 0;
static int spinUPDATE = 0;

/* dummy function call  */
#define so_first() firstORBIT;


int spin( )
{
    static int highest = -1, lowest = -1;
    ORBIT *new,*so_m_serial();
    ATOM *ap,*a_m_serial();

    if( pair == 0 ) pair = 2;
    ap = a_m_serial(i1);
    if( ap == NULL ) {
        aaerror("No atom for an spinal\n"); return 0;
    }
    new = NULL;
    if( highest >= osn && lowest <= osn) new = o_m_serial(osn);
    if( highest < osn ) highest = osn;
    if( lowest > osn ) lowest = osn;
    if( new == NULL )
    {
        spinUPDATE = 1;
        new = malloc( sizeof(ORBIT));
        if( new == NULL )
        {
            aaerror("cannot allocate memory in spinal\n");
            return 0;
        }
        new->next = NULL;
    }
    if( firstORBIT == NULL ) {
        firstORBIT = new;
        highest = osn; lowest = osn; }
    if( lastORBIT == NULL ) lastORBIT = new;
    if( new->next == NULL)
    {
        new->gang = NULL;
        new->next = new;
        new->n = 0;
        new->active = (1==1);
        lastORBIT->next = new;
        lastORBIT = new;
    }

    new->type = type;
    new->myatom = ap;
    new->ncenter = 1;
    new->a1 = NULL;
    new->a2 = NULL;
    new->a3 = NULL;
    new->a4 = NULL;
    new->a5 = NULL;
    if( i2 >= 0 ) new->a1 = a_m_serial(i2);
    if( i3 >= 0 ) new->a2 = a_m_serial(i3);
    if( i4 >= 0 ) new->a3 = a_m_serial(i4);
    if( i5 >= 0 ) new->a4 = a_m_serial(i5);
    if( i6 >= 0 ) new->a5 = a_m_serial(i6);
    if( type == Orm ){ /*  multi atom center */
        if( i2 >= 0 ) new->ncenter += 1;
        if( i3 >= 0 ) new->ncenter += 1;
        if( i4 >= 0 ) new->ncenter += 1;
        if( i5 >= 0 ) new->ncenter += 1;
        if( i6 >= 0 ) new->ncenter += 1;
    }
    new->osn = osn;
    new->along = along;
    new->along2 = x;
    new->x = x;
    new->y = y;
    new->z = z;
    new->spin = spin;
    new->ipair = pair;
    return 1;
}

int dump_spin(outp)
FILE *outp;
{

}






/* useful functions from atom.c */
/* function o_number()
* returns number of spins defined
*  this is just spinNUMBER if spinUPDATE == 0
*  other wise just figure it out
*/
int so_number()
{
    ORBIT *op;
    if( spinUPDATE )
    {
        spinUPDATE = 0;
        spinNUMBER = 0;
        if( firstORBIT == NULL ) return 0 ;
        op = firstORBIT;
        while(1)
        {
            if( op->next == NULL) break;
            spinNUMBER++;
            if( op->next == op ) break;
            op = op->next;
        }
    }
    return spinNUMBER;
}


/* function o_m_serial( serial )
* returns NULL on error or returns the address of the ORBIT
* which matches serial
* cute?
*/
ORBIT *so_m_serial( serial )
int serial;
{
    static ORBIT *op = NULL;
    static ORBIT *lastmatched = NULL;
    int i , n, o_number();

    if( spinUPDATE) n= o_number();
    else n = spinNUMBER;

    op = firstORBIT;
    /* static pointer is hook for more efficient search */
    if( op == NULL) return NULL;
    if( lastmatched == NULL ) lastmatched = firstORBIT;

    if( serial == lastmatched->osn) return lastmatched;
    if( serial > lastmatched->osn) op = lastmatched;
    for( i=0; i< n; i++ )
    {
        if( op-> osn == serial) {lastmatched = op;return op;}
        if( op == op->next)op = firstORBIT ;
        else op = op->next;
    }
    return NULL;
}


/* function o_next( flag )
* returns NULL on error or last spinal
* then steps to the next
* cute?
* flag <= 0 starts it off
*/
ORBIT *so_next( flag )
int flag;
{
    static ORBIT *op = NULL;
    if( op == NULL) op = firstORBIT ;
    if( op == NULL) return NULL;
if( flag <= 0){ op = firstORBIT; return op;}
    if( op == op->next) return NULL;
    op = op->next;
    return op;
}
