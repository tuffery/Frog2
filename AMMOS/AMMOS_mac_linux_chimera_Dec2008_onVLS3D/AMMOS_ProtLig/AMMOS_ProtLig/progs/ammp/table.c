/* table.c
* support table driven potentials
*
*  ammp commands 
*	table id n ;
*   tableent id i r v;
*       dump table ;
*	use tbond;
*       tbond i j id scale;
*       dump tbond;
*
*  table ... is implemented as
*   create_table(), and add_pair_to_table;
*   tables are overwritable
*
*   alse have tbond(),v_tbond(), and f_tbond();
*
*
*
*/
/*
*  copyright 1993,1994,1995,1996 Robert W. Harrison
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

#define MAXENT 201 
typedef struct{
    int id, inme;
    float r[MAXENT],v[MAXENT];
    void *next;
} DATATABLE;
typedef struct{
    ATOM *a1,*a2;
    DATATABLE *which;
    float scale;
    void *next;
} TBOND;
DATATABLE *first_datatable = NULL,*last_datatable;
TBOND *first_tbond = NULL,*last_tbond;

void aaerror( const char *);
ATOM *a_m_serial( int );
int mom_solve( float *, float *, int, int);


#ifdef NEWCALL
int create_table( int id, int inme )
#else
int create_table( id,inme)
int id,inme;
#endif
{
    DATATABLE *new;
    int i;

    if( first_datatable == NULL )
    {
        new = malloc( sizeof( DATATABLE));
        if( new == NULL )
        { aaerror("cannot allocate memory in create_table");
            exit(0);}
        first_datatable = new;
        last_datatable = new;
        last_datatable->next = new;
        new->next = new;

    }else{
        new = first_datatable;
        while( new->id != id)
        {
            if( new->next == new) break;
            new = new->next;
        }
        if( new->next == new && new->id != id){
            new = malloc( sizeof( DATATABLE));
            if( new == NULL )
            { aaerror("cannot allocate memory in create_table");
                exit(0);}
            last_datatable->next = new;
            last_datatable = new;
            new->next = new;
        }
    }
    new->id = id;
    new->inme = inme;
    for( i=0; i< inme; i++)
    { new->r[i] = 0.; new->v[i] = 0.;}

}

#ifdef NEWCALL
int add_pair_to_table( int id,int count, float r, float v)
#else
int add_pair_to_table(  id, count,  r,  v)
int id,count;
float r,v;
#endif
{
    static DATATABLE *who = NULL;

    if( first_datatable == NULL ) return;
    /* this awkward structure avoids NULL->id random referencing */
    if( who == NULL   )
    {
        who = first_datatable;
    }
    if( who->id != id)
    {
        who = first_datatable;
        while( 1==1) {
            if( who->id == id) break;
            if( who->next == who ) return;
            who = who->next;}
    }
    who->r[count] = r;
    who->v[count] = v;

}
#ifdef NEWCALL
int dump_table( FILE *where)
#else
int dump_table( where)
FILE *where;
#endif
{
    DATATABLE *who;
    int i;

    who = first_datatable;
    if( who == NULL) return ;
    while( 1==1){
        fprintf( where,"table %d %d;\n",who->id,who->inme);

        for( i=0; i< who->inme; i++)
            fprintf(where,"tableent %d %d %f %f;\n",who->id,i,who->r[i],who->v[i]);
        fprintf(where,";\n");
        if( who->next == who ) return;
        who = who->next;
    }

}

#ifdef NEWCALL
int tbond( int s1, int s2, int id, float scale)
#else
int tbond(  s1, s2,  id,  scale)
int s1,s2,id;
float scale;
#endif
{
    TBOND *new;
    DATATABLE *table;
    ATOM *a1,*a2;

    a1 = a_m_serial(s1);
    a2 = a_m_serial(s2);
    if( a1 == NULL || a2 == NULL){
        aaerror("undefined atom in tbond"); return 1==0;}
    table = first_datatable;
    if( table == NULL ) return 1 == 0;
    while( 1==1 )
    {
        if( table->id == id ) break;
        if( table->next == table){
            aaerror("undefined bond table in tbond"); return 1==0;}
        table = table->next;
    }

    new = malloc( sizeof( TBOND));
    if( new == NULL )
    { aaerror("cannot allocate memory in tbond");
        exit(0);}
    if( first_tbond == NULL ) first_tbond = new;
    if( last_tbond == NULL ) last_tbond = new;
    last_tbond->next = new;
    last_tbond = new;
    new->next = new;
    new->a1 =a1;
    new->a2 = a2;
    new->scale = scale;
    new->which = table;
    return 1==1;
}


#ifdef NEWCALL
int dump_tbond( FILE *where)
#else
int dump_tbond(where)
FILE *where;
#endif
{
    TBOND *who;
    DATATABLE *which;
    ATOM *a1,*a2;
    who = first_tbond;
    if( who == NULL ) return;
    while( 1==1 )
    {
        a1 = who->a1;
        a2 = who->a2;
        which = who->which;
        fprintf( where,"tbond %d %d %d %f ;\n",
                 a1->serial,a2->serial, which->id, who->scale);
        if( who == who->next) return;
        who = who->next;
    }
}
#ifdef NEWCALL
int v_tbond( float *V, float lambda)
#else
int v_tbond( V,lambda)
float *V,lambda;
#endif
{
    TBOND *bp;
    DATATABLE *tp;
    float r,xt,yt,zt;
    float vl;
    float rred,rx;
    float rclose[3],vclose[3];  /* the three closest values */
    float matrix[3][3],vector[3];
    ATOM *a1,*a2;
    int i,it;

    bp = first_tbond;
    while( 1==1)
    {
        if( bp == NULL ) return 0;

        a1 = bp->a1;
        a2 = bp->a2;
        xt = a1->x - a2->x + lambda*(a1->dx-a2->dx);
        yt = a1->y - a2->y + lambda*(a1->dy-a2->dy);
        zt = a1->z - a2->z + lambda*(a1->dz-a2->dz);
        r = sqrt( xt*xt + yt*yt + zt*zt);
        tp = bp->which;
        if( tp->inme < 3) return tp->v[0]; /* catch really stupid errors */
        if( r < tp->r[0] || r > tp->r[ tp->inme-1])
        {
            if( r < tp->r[0] )      vl = tp->v[0];
            if( r > tp->r[ tp->inme-1] ) vl = tp->v[ tp->inme-1];
        } else {
            /*
            matrix[0][0] = 0.;
            matrix[1][0] = 0.;
            matrix[2][0] = 0.;
            matrix[0][1] = 0.;
            matrix[1][1] = 0.;
            matrix[2][1] = 0.;
            matrix[0][2] = 0.;
            matrix[1][2] = 0.;
            matrix[2][2] = 0.;
            vector[0] = 0.;
            vector[1] = 0.;
            vector[2] = 0.;
            */
            rclose[0] = 10.e10;
            rclose[1] = 10.e10;
            rclose[2] = 10.e10;
            for( i=0; i< tp->inme; i++)
            {rx = (tp->r[i]-r); rred = fabs(rx);
                if( rred <= fabs(rclose[0])){
                    rclose[2] = rclose[1];
                    rclose[1] = rclose[0];
                    rclose[0] = rx;
                    vclose[2] = vclose[1];
                    vclose[1] = vclose[0];
                    vclose[0] = tp->v[i];}
                else if( rred <= fabs(rclose[1])){
                    rclose[2] = rclose[1];
                    rclose[1] = rx;
                    vclose[2] = vclose[1];
                    vclose[1] = tp->v[i];}
                else if( rred <= fabs(rclose[2])){
                    rclose[2] = rx;
                    vclose[2] = tp->v[i];}
            }
            /*
            	for( i=0; i<3; i++)
            	{
            	rred =  rclose[i];
            	{
            		matrix[0][0] += 1.;
            		vector[0] += vclose[i];
            		matrix[0][1] += rred;
            		vector[1] += vclose[i]*rred;
            		rx = rred*rred;
            		matrix[0][2] += rx;
            		matrix[1][1] += rx;
            		vector[2] += rx*vclose[i];
            		rx *= rred;
            		matrix[1][2] += rx;
            		matrix[2][2] += rx*rred;
            	}
            	}

            	matrix[1][0] = matrix[0][1];
            	matrix[2][0] = matrix[0][2];
            	matrix[2][1] = matrix[1][2];

            	mom_solve( &matrix[0][0],&vector[0],3,3);
            	vl = vector[0];
            	*/
            /* lagrange polynomial differences */
            vl = vclose[0]*rclose[1]*rclose[2]/(rclose[0]-rclose[2])/(rclose[0]-rclose[1]);
            vl += vclose[1]*rclose[0]*rclose[2]/(rclose[1]-rclose[2])/(rclose[1]-rclose[0]);
            vl += vclose[2]*rclose[1]*rclose[0]/(rclose[2]-rclose[0])/(rclose[2]-rclose[1]);

            vector[1] = 0.;
            /*	printf(" in v_table %f %f %f %f\n", vl,vclose[0],vclose[1],vclose[2]);
            	printf(" %f %f %f\n",rclose[0],rclose[1],rclose[2]);
            	*/
        } /* end if */
        /*	if( vl < 1.e-5) vl = 1.e-5;*/
        *V -= bp->scale* (vl);
        if( bp == bp->next) return 1==1;
        bp = bp->next;
    }
}
#ifdef NEWCALL
int f_tbond(  float lambda)
#else
int f_tbond(lambda)
float lambda;
#endif
{
    TBOND *bp;
    DATATABLE *tp;
    float r,xt,yt,zt;
    float vl,dvl;
    ATOM *a1,*a2;
    int i,it;
    float rred,rx;
    float rclose[3],vclose[3];
    float matrix[3][3],vector[3];

    bp = first_tbond;
    while( 1==1)
    {
        if( bp == NULL ) return 0;

        a1 = bp->a1;
        a2 = bp->a2;
        xt = a1->x - a2->x + lambda*(a1->dx-a2->dx);
        yt = a1->y - a2->y + lambda*(a1->dy-a2->dy);
        zt = a1->z - a2->z + lambda*(a1->dz-a2->dz);
        r = sqrt( xt*xt + yt*yt + zt*zt);
        if( r < 2.5){
            dvl = 1.;
            if( r < 1.e-4) r = 1.e-4;
            goto SKIP;
        }
        tp = bp->which;
        if( tp->inme < 3) return 0;
        if( r < tp->r[0] || r > tp->r[ tp->inme-1])
        {
            if( r < tp->r[0] )     dvl = 0;
            if( r > tp->r[ tp->inme-1] ) dvl = 0.;
            vl = 1.;
        } else {
            /*
            matrix[0][0] = 0.;
            matrix[1][0] = 0.;
            matrix[2][0] = 0.;
            matrix[0][1] = 0.;
            matrix[1][1] = 0.;
            matrix[2][1] = 0.;
            matrix[0][2] = 0.;
            matrix[1][2] = 0.;
            matrix[2][2] = 0.;
            vector[0] = 0.;
            vector[1] = 0.;
            vector[2] = 0.;
            */
            rclose[0] = 10.e10;
            rclose[1] = 10.e10;
            rclose[2] = 10.e10;
            for( i=0; i< tp->inme; i++)
            {rx = (tp->r[i]-r); rred = fabs(rx);
                if( rred <= fabs(rclose[0])){
                    rclose[2] = rclose[1];
                    rclose[1] = rclose[0];
                    rclose[0] = rx;
                    vclose[2] = vclose[1];
                    vclose[1] = vclose[0];
                    vclose[0] = tp->v[i];}
                else if( rred <= fabs(rclose[1])){
                    rclose[2] = rclose[1];
                    rclose[1] = rx;
                    vclose[2] = vclose[1];
                    vclose[1] = tp->v[i];}
                else if( rred <= fabs(rclose[2])){
                    rclose[2] = rx;
                    vclose[2] = tp->v[i];}
            }

            /*
            	for( i=0; i< 3; i++)
            	{
            	rred = rclose[i];
            	{
            		matrix[0][0] += 1.;
            		vector[0] += vclose[i];
            		matrix[0][1] += rred;
            		vector[1] += vclose[i]*rred;
            		rx = rred*rred;
            		matrix[0][2] += rx;
            		matrix[1][1] += rx;
            		vector[2] += rx*vclose[i];
            		rx *= rred;
            		matrix[1][2] += rx;
            		matrix[2][2] += rx*rred;
            	}
            	}

            	matrix[1][0] = matrix[0][1];
            	matrix[2][0] = matrix[0][2];
            	matrix[2][1] = matrix[1][2];

            	mom_solve( &matrix[0][0],&vector[0],3,3);
            	vl = vector[0];
            	dvl = vector[1];
            	vl = vclose[0]*rclose[1]*rclose[2]/(rclose[0]-rclose[2])/(rclose[0]-rclose[1]);
            	vl += vclose[1]*rclose[0]*rclose[2]/(rclose[1]-rclose[2])/(rclose[1]-rclose[0]);
            	vl += vclose[2]*rclose[1]*rclose[0]/(rclose[2]-rclose[0])/(rclose[2]-rclose[1]);

            	*/
            dvl = vclose[0]*(rclose[1]+rclose[2])/(rclose[0]-rclose[2])/(rclose[0]-rclose[1]);
            dvl += vclose[1]*(rclose[0]+rclose[2])/(rclose[1]-rclose[2])/(rclose[1]-rclose[0]);
            dvl += vclose[2]*(rclose[1]+rclose[0])/(rclose[2]-rclose[0])/(rclose[2]-rclose[1]);

            vector[1] = 0.;

        } /* end if */

        /*	if( vl < 1.e-5) vl = 1.e-5; */
        /*	if( vector[2] > 0 ) dvl = -dvl; */
        dvl = dvl*bp->scale;
SKIP: ;
        xt = -dvl*xt/r;
        yt = -dvl*yt/r;
        zt = -dvl*zt/r;

        if( a1->active ){
            a1->fx += xt;
            a1->fy += yt;
            a1->fz += zt;
        }
        if( a2->active ){
            a2->fx -= xt;
            a2->fy -= yt;
            a2->fz -= zt;
        }


        if( bp == bp->next) return 1==1;
        bp = bp->next;
    }
}
