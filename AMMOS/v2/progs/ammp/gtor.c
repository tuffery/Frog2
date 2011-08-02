/* gtor.c
* 
*  genetic algorithm /torsion search support
*
*  works sort of like tgroup and gene/gdock
*  define a bunch of groups 
*  they do not need a serial number
*  and can be re-written
*
*  gtor will search all groups where 
*   1) atom2 or atom 3 is in range 
*   2) the moving atoms are active 
*
*  it will use the tset routines to do the building
*  it will use a superstud correction 
*/

/* copyright 2000 robert w harrison
*  GPL license granted 
*/

#include  <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "ammp.h"


typedef struct {
    ATOM *a1,*a2,*a3,*a4;
    int nstep;
    float base;
    void *next;
}  GTGROUP;

GTGROUP *first_GTGROUP = NULL;
GTGROUP *last_GTGROUP = NULL;
#define GLONG sizeof(GTGROUP);
#define MAX_PROGENY 4

/*
int gtor_tset(  ATOM *, ATOM *, ATOM*, ATOM *, float );
int gtor_define( int ,int, int, int,  float,int);
int dump_gtor( FILE *);
int ator_gene( void *, int,int (*())[],int , int, int, int, int, float);
float randf();
float get_torsion_value( ATOM *, ATOM *, ATOM *, ATOM*);
int tset_bond_build( ATOM *, ATOM *, ATOM *, ATOM *);
int set_torsion(ATOM *, ATOM *,ATOM *,ATOM *,float);
*/

#ifdef NEWCALL
int gtor_gene( FILE *op, int echo, int (*pots[])(), int npot,
               int nstep, int ngene, int low, int high, float toler)
#else
int gtor_gene( op,echo,pots,npot,nstep,ngene,low,high,toler)
FILE *op;
int echo,npot,nstep,ngene,low,high;
float toler;
int (*pots[])();
#endif
{
    int (*studbook)[];
    int (*genome)[];
    int (*g2)[];
    float (*fvals)[];
    float (*fitness)[];
    GTGROUP *(*groups)[],*gp;
    int (*bestgene)[];
    ATOM *a1,*a2,*a3,*a4;
    int ingroups,iter,i,j,best,igene;
    float Vbest;
    float vmin,vmax,vmean;
    float randf();
#ifdef GTK_VERSION
    void AMMP_whiz_bang();
#endif
    /* first we have to find howmany gtgroups will be used  (well
    don't we? */
    Vbest = 10.e20;
    if( ngene < 1) ngene = 1;
    ingroups = 0;
if( low > high) { i = low; low = high; high = i;}
    gp = first_GTGROUP;
    if( gp == NULL ) return 1==1; /* nothing defined for us */
    while(1==1)
    {
        a1 = gp->a1;
        a2 = gp->a2;
        a3 = gp->a3;
        a4 = gp->a4;
        if( (a2->serial >low && a2->serial < high) ||
                (a3->serial >low && a3->serial < high) ) ingroups += 1;
        if( gp == gp->next) break;
        gp = gp->next;
        if( gp == NULL) break;
    }
    if( ingroups == 0 ) return 1==1;
    /* so now malloc up the memory for the lists */
    groups = malloc( ingroups * sizeof( GTGROUP *) );
    i = 0;
    gp = first_GTGROUP;
    while( 1==1)
    {
        a1 = gp->a1;
        a2 = gp->a2;
        a3 = gp->a3;
        a4 = gp->a4;
        if( (a2->serial >low && a2->serial < high) ||
                (a3->serial >low && a3->serial < high) )
        {
            (*groups)[i++] = gp;
        }
        if( gp == gp->next) break;
        gp = gp->next;
        if( gp == NULL) break;
    }
    /* so now we have the active groups defined
    *  so let's allocate the other memory
    */
    genome = malloc( (ngene)*ingroups*sizeof(int));
    g2 = malloc( (ngene)*ingroups*sizeof(int));
    bestgene = malloc(ingroups*sizeof(int));
    studbook = malloc( ngene*sizeof(int));
    fvals = malloc( ngene*sizeof(float));
    fitness = malloc( ngene*sizeof(float));
    /* initiallize the genome */
    for( i=0; i< ingroups; i++)
    {
        gp = (*groups)[i];
        for(j=0; j< ngene; j++)
        {
            (*genome)[j*ingroups+i]=(int)(gp->nstep*randf());;
            if( (*genome)[j*ingroups+i] < 0) (*genome)[j*ingroups+i] = 0;
            if( (*genome)[j*ingroups+i] >= nstep)
                (*genome)[j*ingroups+i] = nstep-1;
        }
    }/*i */
    /* now do the work
    *   calculate the scores
    * select the best 
    * recombine the genome
    */
    for( iter=0; iter< nstep; iter++)
    {
        for( j=0;j< ngene; j++)
        {
            gtor_build( groups, genome,j*ingroups , ingroups);
            (*fvals)[j] = 0.;
            for( i=0; i< npot; i++)
                (*pots[i])( &(*fvals)[j], 0.);
        }
        for( j=0;j< ngene; j++)
        {
            if( Vbest > (*fvals)[j])
            {
                Vbest = (*fvals)[j];
                for(i=0; i< ingroups; i++)
                    (*bestgene)[i] = (*genome)[j*ingroups+i];
            }
            (*studbook)[j] = 0; /* why not ?*/
        }
#ifdef GTK_VERSION

        gtor_build( groups, bestgene,0, ingroups);
        AMMP_whiz_bang();
#endif
        for( j=0; j< ngene*ingroups; j++)
            (*g2)[j] = (*genome)[j];
        /* now assign a fitness
        *  given by (normalize)( 1.- (v-vmin)/(vmax-vmin));
        */
STUD_AGAIN: ;
        vmin = 10.e10; vmax = -vmin;
        for( j=0; j< ngene; j++)
        {
            if((*studbook)[j] < MAX_PROGENY) break;
        }
        if( j == ngene)
        {
            for( j=0; j< ngene; j++)
                (*studbook)[j] -= 1;
            goto STUD_AGAIN;
        }

        for( j=0;j< ngene; j++)
        {
            if( (*studbook)[j] < MAX_PROGENY)
            {
                if( vmax < (*fvals)[j]) vmax = (*fvals)[j];
                if( vmin > (*fvals)[j]) vmin = (*fvals)[j];
            }
        }
        if( fabs(vmax-vmin) < 1.e-2 ) break; /* this will skip to the return */
        if( echo) fprintf(op," %f %f %f\n", vmax, vmin, Vbest);
        for( j=0; j< ngene; j++)
        {
            if( (*studbook)[j] < MAX_PROGENY)
                (*fitness)[j] = (1.- ((*fvals)[j]-vmin)/(vmax-vmin));
        }
        vmin = 0.;
        for( j=0; j< ngene; j++)
            if( (*studbook)[j] < MAX_PROGENY) vmin += (*fitness)[j];
        vmin = 1./vmin;
        for( j=0; j< ngene; j++)
            if( (*studbook)[j] < MAX_PROGENY) (*fitness)[j] *= vmin;

        /* this is the selection phase */

        for(best = 0; best< ngene; best++)
        {
RE_SELECT: ;
            for( igene=-4; igene < 0; igene++)
            {
                for( j=0; j< ngene; j++)
                    if( (*studbook)[j] < MAX_PROGENY && randf() < (*fitness)[j])
                    {igene = j; goto SELECTED; (*studbook)[j]+= 1;}
            }
            if( j == ngene) goto DO_IT_ANYWAY;
            if( igene == 0) goto DO_IT_ANYWAY;
SELECTED: ;
            if( igene > ngene) goto RE_SELECT ;

            for( i=0; i< ingroups; i++)
            {
                (*g2)[best*ingroups + i] = (*genome)[igene*ingroups+i];
                if( randf() < toler) (*g2)[best*ingroups+i] += 1;
                else if( randf() < toler) (*g2)[best*ingroups+i] -= 1;
                if( (*g2)[best*ingroups+i] > (*groups)[i]->nstep)
                    (*g2)[best*ingroups+i] -= (*groups)[i]->nstep;
                if( (*g2)[best*ingroups+i] < 0)
                    (*g2)[best*ingroups+i] += (*groups)[i]->nstep;
                if( randf() < toler*toler) goto RE_SELECT;
            }

        } /* best */
        /* this is after the selection phase */
DO_IT_ANYWAY: ; /* g2 is initialized to genome */
        for( j=0; j< ngene*ingroups; j++)
            (*genome)[j] = (*g2)[j];

    }/* iter */

    gtor_build( groups, bestgene,0, ingroups);
#ifdef GTK_VERSION
    AMMP_whiz_bang();
#endif

    free( fitness); /* last malloc first free */
    free( fvals);
    free( studbook);
    free(bestgene);
    free( g2);
    free( genome);
    free( groups); /* first malloc last free */

    return 1==1;
}

/* build the molecule from the given genome */
#ifdef NEWCALL
/*int gtor_build( GTGROUP *(*groups)[], int gene[], int ingene )
*/
#else
int gtor_build( groups,  gene, start, ingene )
GTGROUP *(*groups)[];
int (*gene)[],ingene,start;
#endif
{/* ingene is what is called ingroups in gtor */
    int i;
    GTGROUP *gp;
    float alpha;
    int gtor_tset();

    for( i=0; i< ingene; i++)
    {
        /*
        int gtor_tset(ATOM *ap1,ATOM *ap2,ATOM *ap3,ATOM *ap4,float  alpha )
        */
        gp = (*groups)[i];
        alpha = gp->base + (*gene)[i+start]*TWOPI/(gp->nstep);
        gtor_tset( gp->a1,gp->a2,gp->a3,gp->a4,alpha);

    }

    return 1==1;
}


#ifdef NEWCALL
int gtor_define( int i1,int i2,int i3,int i4,  float base, int nstep)
#else
int gtor_define( i1,i2,i3,i4,base,nstep)
int i1,i2,i3,i4,nstep;
float base;
#endif
{
    GTGROUP *new;
    ATOM *a1,*a2,*a3,*a4,*a_m_serial();
    int a_number();

    if( a_number() < 4) return 1==0;
    a1 = a_m_serial(i1);
    a2 = a_m_serial(i2);
    a3 = a_m_serial(i3);
    a4 = a_m_serial(i4);
    if( a1 == NULL || a2 == NULL || a3 == NULL || a4 == NULL) return 1==0;

    new = first_GTGROUP;
    if( new != NULL)
    {
        while( 1== 1)
        {
            if( new->a1 == a1 && new->a2 == a2 && new->a3 == a3 && new->a4 ==  a4)
                goto USE_ME;
            if( new == new->next) break;
            new = new->next;
            if( new == NULL) break;
        }

    }
    /* if here we need a new element */
    new = malloc( sizeof(GTGROUP));
    if( new == NULL ) return 1==0;
    if( first_GTGROUP == NULL)
    {
        first_GTGROUP = new;
        last_GTGROUP = new;
    }
    last_GTGROUP->next = new;
    last_GTGROUP = new;
    new->next = new;
USE_ME: ; /* if we jump here, we're just rewriting the element */

    new->a1 = a1;
    new->a2 = a2;
    new->a3 = a3;
    new->a4 = a4;
    new->base = base*PI/180.;
    new->nstep = nstep;

    return 1==1;
}

int dump_gtor( FILE *op)
{
    GTGROUP *gp;
    ATOM *a1,*a2,*a3,*a4;
    gp = first_GTGROUP;
    if( gp == NULL ) return 1==1;
    while( 1==1)
    {
        a1 = gp->a1;
        a2 = gp->a2;
        a3 = gp->a3;
        a4 = gp->a4;
        fprintf(op,"gtgroup %d %d %d %d %f %d;\n",
                a1->serial,a2->serial,a3->serial,a4->serial,
                gp ->base*180./PI, gp->nstep);

        if( gp == gp->next) break;
        gp = gp->next;
        if( gp == NULL) break;
    }

    return 1==1;
}



#ifdef NEWCALL
int gtor_tset(ATOM *ap1,ATOM *ap2,ATOM *ap3,ATOM *ap4,float  alpha )
#else
int gtor_tset(ap1,ap2,ap3,ap4,alpha )
ATOM *ap1,*ap2,*ap3,*ap4;
float alpha;
#endif
{
    ATOM  *a_m_serial(),*a_next();
    float get_torsion_value();
    float original,delta;
    int numatom,a_number();
    int tset_bond_build();

    /* there must be atoms and they must be specified */
    numatom = a_number();
    if( numatom <= 0) return ;

    /* get the torsion value and check if we have to do the work */

    original = get_torsion_value( ap1,ap2,ap3,ap4);
    if( fabs(original  -alpha) < 1.e-3 ) return;

    /* set up the bond structure flags in ap->gx */
    tset_bond_build( ap1,ap2,ap3,ap4);

    delta = alpha -original;
    set_torsion( ap1,ap2,ap3,ap4,delta );

}/* end of gtor_tset */
