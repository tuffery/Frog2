/*
*  copyright 1994 Robert W. Harrison
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

/* routine list
*
*  xtalcell (a,b,c,alpha,beta,gamma) - set the cell 
*  dump_xtalcell( FILE *fp )  output the cell
*  scatter( char *atom, int id , float a1,b1,...a4,b4,c) read the scatter 
*  dump_scatter( FILE *fp) output the scatter factor
*  symmetry()
*  symmetry_parse()
*  dump_symmetry()
*  reflection()
*  dump_reflection() 
*
*  float getsthol( h,k,l) (returns sthol**2 )
*  int tofrac( x,y,z, &xf,&yf,&zf )  returns 0 if no unit cell
*  float scatter_factor( id, sthol^2) return scattering factor
*  
*/


typedef struct { float a,b,c,alpha,beta,gamma;
    float rtens[3][3], ftens[3][3];} UNIT_CELL;
typedef struct { int h,k,l; float fo,sigma,fca,fcb,sthol;
    float s[10];
    void *next;} REFLECTION;
typedef struct { char atom; int id; float a1,b1,a2,b2,a3,b3,a4,b4,c;
    void *next;} SCATTER;
#define MAXSYMOPER 50
#define MAXSCATT 20 
typedef struct {int id; char operation[MAXSYMOPER]
    ;float matrix[3][3], vector[3]; void *next; } SYMMETRY;

UNIT_CELL  theunitcell = {-1.,-1,-1}; /* the one and only one !!! */

REFLECTION *firstREFL = NULL;
REFLECTION *lastREFL = NULL;
SCATTER *firstSCATTER = NULL;

SYMMETRY IDENTITY = { /*always define x,y,z */
    0,"X,Y,Z",1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,NULL};
SYMMETRY *firstSYMMETRY = &IDENTITY;
/*
IDENTITY.id = 0;
IDENTITY.operation = "X,Y,Z";
IDENTITY.matrix[0][0] = 1.;
IDENTITY.matrix[0][1] = 0.;
IDENTITY.matrix[0][2] = 0.;
IDENTITY.matrix[1][0] = 0.;
IDENTITY.matrix[1][1] = 1.;
IDENTITY.matrix[1][2] = 0.;
IDENTITY.matrix[2][0] = 0.;
IDENTITY.matrix[2][1] = 0.;
IDENTITY.matrix[2][2] = 1.;
IDENTITY.vector[0] = 0.;
IDENTITY.vector[1] = 0.;
IDENTITY.vector[2] = 0.;
IDENTITY.next = NULL;
*/

/* define the cell parameters */
int xtalcell( a,b,c,alpha,beta,gamma)
float a,b,c,alpha,beta,gamma;
{
    float s[3][3],d1,d2;

    theunitcell.a = a;
    theunitcell.b = b;
    theunitcell.c = c;
    if( alpha <= 0.) alpha = 90.;
    if( beta <= 0.) beta = 90.;
    if( gamma <= 0.) gamma = 90.;
    theunitcell.alpha = alpha*3.14159265/180.;
    theunitcell.beta = beta*3.14159265/180.;
    theunitcell.gamma = gamma*3.14159265/180.;
    /* the forward tensor */
    s[0][0] = a*a;
    s[0][1] = a*b*cos(theunitcell.gamma);
    s[0][2] = a*c*cos(theunitcell.beta);
    s[1][1] = b*b;
    s[1][2] = b*c*cos(theunitcell.alpha);
    s[2][2] = c*c;
    s[2][0] = s[0][2]; /* being pedantic aren't we? */
    s[2][1] = s[1][2];
    s[1][0] = s[0][1];

    d1  = s[0][0]*( s[1][1]*s[2][2]- s[1][2]*s[2][1]);
    d1 -= s[0][1]*( s[1][0]*s[2][2]- s[1][2]*s[2][0]);
    d1 += s[0][2]*( s[1][0]*s[2][1]- s[1][1]*s[2][0]);
    if( d1 <= 1.e-7)
    {
        aaerror(" this unit cell has zero volume and is illegal \n");
        return;
    }
    d2 = 1./d1;

    /* cofactor expansion of the inverse */
    theunitcell.rtens[0][0] = d2*( s[1][1]*s[1][2]-s[1][2]*s[1][2]);
    theunitcell.rtens[0][1] = -d2*( s[1][0]*s[2][2]-s[2][0]*s[1][2]);
    theunitcell.rtens[0][2] = d2*( s[1][0]*s[2][1]-s[2][0]*s[1][1]);
    theunitcell.rtens[1][0] = -d2*( s[1][0]*s[2][2]-s[2][0]*s[1][2]);
    theunitcell.rtens[2][0] = d2*( s[1][0]*s[2][1]-s[2][0]*s[1][1]);
    theunitcell.rtens[1][1] = d2*( s[0][0]*s[2][2] - s[2][0]*s[0][2]);
    theunitcell.rtens[1][2] = -d2*( s[0][0]*s[2][1] - s[2][0]*s[0][1]);
    theunitcell.rtens[2][1] = -d2*( s[0][0]*s[2][1] - s[2][0]*s[0][1]);
    theunitcell.rtens[2][2] = d2*( s[0][0]*s[1][1] - s[1][0]*s[0][1]);

    /* copied from VECTOR.F */
    theunitcell.ftens[0][0] = 1.;
    theunitcell.ftens[0][1] = 0.;
    theunitcell.ftens[0][2] = 0.;
    theunitcell.ftens[1][0] = cos(theunitcell.gamma);
    theunitcell.ftens[1][1] = sin(theunitcell.gamma);
    theunitcell.ftens[1][2] = 0.;
    theunitcell.ftens[2][0] = cos(theunitcell.beta);
    theunitcell.ftens[2][1] = cos(theunitcell.alpha)-theunitcell.ftens[2][0]*
                              theunitcell.ftens[1][0]/theunitcell.ftens[1][1];
    theunitcell.ftens[2][2] = sqrt(theunitcell.ftens[2][0]- theunitcell.ftens[2][1])
                              *(theunitcell.ftens[2][0]+theunitcell.ftens[2][1]);

    theunitcell.ftens[0][0] *= 1./a;
    theunitcell.ftens[1][0] *= 1./b;
    theunitcell.ftens[1][1] *= 1./b;
    theunitcell.ftens[2][0] *= 1./c;
    theunitcell.ftens[2][1] *= 1./c;
    theunitcell.ftens[2][2] *= 1./c;

}

int tofrac( x,y,z,xf,yf,zf)
float x,y,z;
float *xf,*yf,*zf;
{
    if( theunitcell.a <=0 )
    {
        aaerror("There is no unit cell defined\n");
        return 0;
    }
    *xf = theunitcell.ftens[0][0]*x + theunitcell.ftens[1][0]*y + theunitcell.ftens[2][0]*z;
    *yf = theunitcell.ftens[1][1]*y + theunitcell.ftens[1][2]*z;
    *zf = theunitcell.ftens[2][2]*z ;

    return 1;
}

float getsthol( h,k,l )
int h,k,l;
{

    float x;

    x = theunitcell.rtens[0][0]*h*h + theunitcell.rtens[0][1]*h*k  +theunitcell.rtens[0][2]*h*l ;
    x += theunitcell.rtens[1][0]*k*h + theunitcell.rtens[1][1]*k*k  +theunitcell.rtens[1][2]*k*l ;
    x += theunitcell.rtens[2][0]*l*h + theunitcell.rtens[2][1]*l*k  +theunitcell.rtens[2][2]*l*l ;

    return x;
}

int dump_xtalcell( fp )
FILE *fp;
{
    fprintf( fp,"UNITCELL %f %f %f %f %f %f ;\n",
             theunitcell.a,theunitcell.b,theunitcell.c,
             theunitcell.alpha*180./3.14159265, theunitcell.beta*180./3.14159265,
             theunitcell.gamma*180./3.14159265 );
}


int scatter( atom, number, a1,b1,a2,b2,a3,b3,a4,b4,c)
char *atom;
int number;
float a1,b1,a2,b2,a3,b3,a4,b4,c;
{
    SCATTER *sp,*spold,*spnew;

    sp = firstSCATTER;
    spold = NULL;
    spnew = NULL;
    while( sp != NULL)
    {
        if( sp ->id == number){ spnew = sp; break;}
        spold = sp;
        sp = sp->next;
    }
    if( spnew == NULL)
    {
        spnew = malloc( sizeof( SCATTER) );
        if( spnew == NULL){
            aaerror(" cannot allocate memory for SCATTER\n");
            exit(0);
        }
        if( spold == NULL) {firstSCATTER = spnew;}
        else{  spold->next = spnew;}
        spnew->next = NULL;
    }

    spnew->atom = *atom;
    spnew->id = number;
    spnew->a1 = a1;
    spnew->a2 = a2;
    spnew->a3 = a3;
    spnew->a4 = a4;
    spnew->b1 = b1;
    spnew->b2 = b2;
    spnew->b3 = b3;
    spnew->b4 = b4;
    spnew->c = c;
}

int dump_scatter( fp )
FILE *fp;
{
    SCATTER *sp;

    sp = firstSCATTER;
    while( sp != NULL)
    {
        fprintf(fp,"SCATTER %c %d %f %f %f %f %f %f %f %f %f ;\n",
                sp->atom, sp->id,
                sp->a1,sp->b1,sp->a2,sp->b2,sp->a3,sp->b3,sp->a4,sp->b4,sp->c);
    }
}

float scatter_factor( id, sthol)
float sthol; /* sthol ^2 !!! */
int id;
{
    SCATTER *sp;
    float f;

    sp = firstSCATTER;
    while( sp != NULL)
    {
        if( sp->id == id) {
            /* int tables old Vol IV section 2.2 equation 3
            * sthol is sin(theta)/lambda squared    */
            f = sp->c + sp->a1*exp( -sp->b1*sthol);
            f +=  sp->a2*exp( -sp->b2*sthol);
            f +=  sp->a3*exp( -sp->b3*sthol);
            f +=  sp->a4*exp( -sp->b4*sthol);
            return f;
        }
    }
    return 0.;
}

int symmetry(id,oper,fp,echo )
int id;
char *oper;
int echo;
FILE *fp;
{
    SYMMETRY *sp,*spold,*spnew;
    char *cp,*cp2;
    int i;

    sp = firstSYMMETRY;
    spold = NULL;
    spnew = NULL;
    while( sp != NULL)
    {
        if( sp ->id == id){ spnew = sp; break;}
        spold = sp;
        sp = sp->next;
    }
    if( spnew == NULL)
    {
        spnew = malloc( sizeof( SYMMETRY) );
        if( spnew == NULL){
            aaerror(" cannot allocate memory for SYMMETRY\n");
            exit(0);
        }
        if( spold == NULL) {firstSYMMETRY = spnew;}
        else{  spold->next = spnew;}
        spnew->next = NULL;
    }

    spnew->id = id;
    cp = &spnew->operation[0];
    cp2 = oper;
    i = 0;
    while(*cp2 != '\0' && i < MAXSYMOPER)
    {  *(cp++) = *(cp2++);  i++;}
    if( i == MAXSYMOPER ){
        aaerror(oper);
        aaerror(" is too long to be a valid symmetry operator\n");
        exit(0);
    }
    /* now need to initialize the matrix and vector */

    if(
        symmetry_parse( oper, &spnew->matrix[0][0],&spnew->vector[0],fp,echo)
    )
    { if(spold != NULL)spold->next = NULL; free(spnew); }

}


int symmetry_parse( oper,matrix,vector,fp,echo)
char *oper;
float matrix[3][3],vector[3];
int echo;
FILE *fp;
{
    char local[MAXSYMOPER+1],*cp;
    int i,j,incomma,k,l,m;
    int comma[5];
    float x;
    /* copy to local, tolower, and remove spaces, flag for 3 operators
    * and keep the ',' places */
    cp = oper;
    i = 0;
    local[0] = ','; /* this is always ok if start parse from 1 */
    j = 1;
    incomma = 1;
    comma[0] = 0;
    for( i=0 ; i< MAXSYMOPER; i++)
    {
        if( isupper(*cp) )
        {
            local[j++] = tolower(*cp);
        } else if( *cp == ',')
        { comma[incomma++] = j; local[j++] = *cp; }
        else if( *cp != ' ') local[j++] = *cp;
        local[j] = '\0';
        cp ++;
        if( incomma > 4){
            aaerror(oper);
            aaerror(" has too many commas\n");
            return 1;
        }
        if( *cp == '\0'){ comma[incomma++] = j; break;}
    }
    if( incomma < 4)
    { aaerror(oper); aaerror(" ill formed symmetry operator\n"); return 1; }

    /* local is now ,<>,<>,<><,>\0
    * leading + will only enter iff local[1] == (x,y,z) 
    *  comma points to the (max) 4 commas or 3 commas and \0
    */
    /* symmetry outer loop */
    for( i=0;i<3; i++)
    {
        matrix[i][0] = 0;
        matrix[i][1] = 0;
        matrix[i][2] = 0;
        vector[i] = 0;
        /* find x */
        for( j=comma[i]+1; j< comma[i+1]; j++)
        {
            if( local[j] == 'x')
            {
                if( local[j-1] == '-') { matrix[i][0] = -1;}
                else {matrix[i][0] = 1.;}
            }
        }
        /* find y */
        for( j=comma[i]+1; j< comma[i+1]; j++)
        {
            if( local[j] == 'y')
            {
                if( local[j-1] == '-') { matrix[i][1] = -1;}
                else {matrix[i][1] = 1.;}
            }
        }
        /* find z  */
        for( j=comma[i]+1; j< comma[i+1]; j++)
        {
            if( local[j] == 'z')
            {
                if( local[j-1] == '-') { matrix[i][2] = -1;}
                else {matrix[i][2] = 1.;}
            }
        }
        /* find translation */
        for( j=comma[i]+1; j< comma[i+1]; j++)
        {
            if( local[j] == '/' || local[j] == '\\')
            {
                for( k=j-1; k> comma[i]; k--)
                {
                    if( !isdigit(local[k]) && local[k] != '-' && local[k] != '+') break;
                }
                /*
                	if( k == comma[i])
                	{ aaerror(oper); aaerror(" ill formed symmetry operator\n"); return 1; }
                */

                if( local[j] == '/')
                    sscanf(&local[k+1],"%d/%d",&l,&m);
                else
                    sscanf(&local[k+1],"%d\\%d",&l,&m);

                vector[i] = (float)l/(float)m;
                break;
            } /* if '/' */
        } /* for j */
    }/* end of outer loop (i) */

    if( echo )
    {
        fprintf(fp," operator %s parsed to :\n",oper);
        fprintf( fp," %f %f %f %f\n",matrix[0][0],matrix[0][1],matrix[0][2],vector[0]);
        fprintf( fp," %f %f %f %f\n",matrix[1][0],matrix[1][1],matrix[1][2],vector[1]);
        fprintf( fp," %f %f %f %f\n",matrix[2][0],matrix[2][1],matrix[2][2],vector[2]);
    }

    x  = matrix[0][0]*(matrix[1][1]*matrix[2][2]- matrix[2][1]*matrix[1][2]);
    x -= matrix[0][1]*(matrix[1][0]*matrix[2][2]- matrix[2][0]*matrix[1][2]);
    x += matrix[0][2]*(matrix[1][0]*matrix[2][1]- matrix[2][0]*matrix[1][1]);
    if( x < 0 ) x = -x;
    if( x < .999 || x > 1.01)
    { aaerror(oper); aaerror(" ill formed symmetry operator abs(det) != 1\n"); return 1; }
    return 0;

}


int dump_symmetry( fp )
FILE *fp;
{
    SYMMETRY *sp;
    sp = firstSYMMETRY;
    while ( sp != NULL)
    {
        fprintf(fp,"SYMMETRY %d %c%s%c;\n", sp->id,'"',sp->operation,'"');
        sp = sp->next;
    }
}



int reflection( h,k,l,fo,sigma,fa,fb )
int h,k,l;
float fo,sigma,fa,fb;
{
    REFLECTION *sp,*spnew;

    sp = lastREFL;
    spnew = malloc( sizeof( REFLECTION) );
    if( spnew == NULL){
        aaerror(" cannot allocate memory for a REFLECTION\n");
        exit(0);
    }
    if( sp == NULL) {firstREFL = spnew;}
    else{  sp->next = spnew;}
    lastREFL = spnew;
    spnew->next = NULL;
    spnew->h = h;
    spnew->k = k;
    spnew->l = l;
    spnew->fo = fo;
    if( sigma <= 1.e-6) sigma = 1.;
    spnew->sigma = sigma;
    spnew->fca = fa;
    spnew->fcb = fb;
    spnew->sthol = -1; /* negative sthol cannot be used */
}

int dump_reflection( fp )
FILE *fp;
{
    REFLECTION *rp;
    rp = firstREFL;
    while (rp->next != NULL)
    {
        fprintf(fp,"REFLCT %d %d %d %f %f %f %f;\n",
                rp->h,rp->k,rp->l,rp->fo,rp->sigma,rp->fca,rp->fcb);
        rp = rp->next;
    }
}

/* if dot notation is used assign the scatterers for each
* atom in the range ilow to ihi
*/
int assign_scatter( ilow,ihi)
int ilow,ihi;
{
    ATOM *ap,*a_next();
    SCATTER *sp;
    int i,j,numatm,a_number();
    char *cp;

    numatm = a_number();
    if( numatm <= 0 ) return ;
    if( firstSCATTER == NULL ) return;
    if( ilow > ihi)
    { i = ilow, ilow = ihi; ihi = i;}

    for( i=0; i< numatm; i++)
    {
        ap = a_next(i);
        if( ap->serial >= ilow && ap->serial <= ihi)
        {
            cp = &ap->name[0];
            while(*cp != '\0')
            {
                if( *(cp) == '.') break; cp++;
            }
            if( *cp == '.')
            {
                cp++;
                ap->iscatter = -1; /* don't let it match anyone */
                sp = firstSCATTER;
                j = 0;
                while( sp != NULL)
                {
                    if( sp->atom == *cp )
                    {
                        ap->iscatter = sp->id;
                        break;
                    }
                    sp = sp ->next;
                } /* while sp */
            } /* if *cp == . */
        }/*atom bounds */
    } /* atom loop */
}/*end of routine */

int v_xray( V,lambda )
float *V, lambda;
{
    float x,y,z,w;
    float  scatter_factor();
    int tofrac();
    /*  int tofrac( x,y,z, &xf,&yf,&zf )  returns 0 if no unit cell
    *  float scatter_factor( id, sthol^2) return scattering factor
    */
    ATOM *ap,*a_next();
    REFLECTION *rp;
    SYMMETRY *sp;
    int numatom,a_number();
    int i, j,inftype;
    int ftype[MAXSCATT];
    float fvalue[MAXSCATT];
    float fca,fcb;

    float vt,v1,v2;


    numatom = a_number();
    if( numatom <= 0 ) return;

    /* skip identity */
    sp = firstSYMMETRY->next;

    /* before doing anything calculate the
    * current fractional coordinates 
    */
    for( i=0; i< numatom; i++)
    {
        ap = a_next(i);
        x = ap->x + lambda*ap->dx;
        y = ap->y + lambda*ap->dy;
        z = ap->z + lambda*ap->dz;
        tofrac( x,y,z , &ap->xf,&ap->yf,&ap->zf);
    }

    vt = 0.;
    v1 = 0.;
    v2 = 0.;

    /* first calculate the fc and do the data */
    rp = firstREFL;
    while( rp != NULL )
    {
        if( rp->sthol > 0. )
        {
            fca = 0.; fcb = 0.;
            inftype = 0;
            /* loop over the atoms */
            for( i=0; i< numatom; i++)
            {
                ap = a_next(i);
                /* loop over the scattering types which have been calculated */
                for( j=0; j< inftype; j++)
                {
                    if( ftype[j] == ap->iscatter) break;
                }
                if( j == inftype)
                {
                    inftype ++;
                    ftype[j] = ap->iscatter;
                    fvalue[j] = scatter_factor( ap->iscatter, rp->sthol);
                }
                x = ap->xf; y = ap->yf; z = ap->zf;

                fca +=
                    fcb +=

                    }
                } /* checking if for doing the reflection */
            } /* while loop over the reflections*/

            /* now do the r2factor  */
            rp = firstREFL;
    while( rp != NULL )
    {
        if( rp->sthol > 0. )
        {
            fca = 0.; fcb = 0.;


        } /* checking if for doing the reflection */
    } /* while loop over the reflections*/

}/* end of v_xray */
