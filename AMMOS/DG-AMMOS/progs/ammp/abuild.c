/* abuild.c
*
* analytic building routine
*
* use resultants to build atoms based on 3 or more distances
* 
* use distance/planarity to build based on 1 or 2 distances
*
*  requires at least one atom not 0 0 0 
* and active.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ammp.h"
#include "numeric.h"

#define MAX_LOCAL 800
int a_build(nused,pot,force, ncycles,low_atom,high_atom,op,echo)
int nused, (*pot[])(), ncycles,low_atom,high_atom,echo;
int (*force[])();
FILE *op;
{
    ATOM *ap,*bp,*bps[MAX_LOCAL],*a_next();
    int numatom,a_number();
    int i,j,k,icycle;
    int iatom, inbps;
    float x,y,z,r,rc,randf();
    float a,b,c,d;
    float honker[MAX_LOCAL][3];
    float m[82],v[MAX_LOCAL];
    float m2[82],v2[10],m3[10];
    float rvec[3],svec[3];
    int v_bond(),v_angle(),v_noel(),v_noe_lin();
    int v_mmbond(),v_mmangle(),v_c_angle();
    int v_hybrid();
    int v_step();
    int v_hard(),v_nonbon(),u_v_nonbon();
    int f_bond(),f_angle(),f_noel();
    int kohonen_minimizer();
    int (*lp[3])(),(*lf[3])(),ilocal;
    int use_noel, use_bond,use_angle, use_hybrid, use_nonbon;
    void rand3(),gsdg();
#ifdef WINDOZ
    void force_screen_update(); /* in animate.c */
#endif
#ifdef GTK_VERSION
    void AMMP_whiz_bang();
#endif
    float abuild_line_search();
    /* gsdg_xxx() put r**2 in ap->vx  gsdg_hybrid sets hybrids*/
    int gsdg_angle(),gsdg_bond(),gsdg_noel(), gsdg_hybrid();

    numatom = a_number();
    if( numatom < 2) return 0; /* nothing to build if 0 or 1 atom */

    if( low_atom <= 0) low_atom = 0;
    if( high_atom <=0 ) high_atom = 10000000;
    if( low_atom > high_atom)
    { i=low_atom; low_atom= high_atom; high_atom=i;}
    if( ncycles < 1) ncycles = 1;

    use_bond = (1==0);
    use_angle = use_bond;
    use_noel = use_bond;
    use_hybrid= use_bond;
    use_nonbon = use_bond;
    for( i=0; i< nused; i++)
    {
        if( pot[i] == v_bond) use_bond = (1==1);
        if( pot[i] == v_mmbond) use_bond = (1==1);
        if( pot[i] == v_angle) use_angle = (1==1);
        if( pot[i] == v_c_angle) use_angle = (1==1);
        if( pot[i] == v_mmangle) use_angle = (1==1);
        if( pot[i] == v_noel) use_noel = (1==1);
        if( pot[i] == v_noe_lin) use_noel = (1==1);
        if( pot[i] == v_hybrid) use_hybrid = (1==1);
        if( pot[i] == v_hard) use_nonbon = (1==1);
        if( pot[i] == v_nonbon) use_nonbon = (1==1);
        if( pot[i] == u_v_nonbon) use_nonbon = (1==1);
    }
    ilocal = 0;
if( use_bond) {lf[ilocal] = f_bond; lp[ilocal++] = v_bond;}
    if( use_angle) {lf[ilocal] = f_angle; lp[ilocal++] = v_angle;}
    if( use_noel) {lf[ilocal] = f_noel; lp[ilocal++] = v_noel;}
    if( ilocal == 0) return 0; /* no distances will be generated*/

    for( icycle=0; icycle< ncycles; icycle ++)
    {
        ap = a_next(-1);

        for( iatom =0; iatom < numatom; iatom ++)
        {
            if( !ap->active) goto SKIP;
            /* zero the gsdg return registers */
            for( j=0; j< numatom; j++)
            { bp = a_next(j);
                bp->vx = 0.;
                bp->vy = 0.;
            }

            if( use_bond) gsdg_bond(ap);
            if( use_angle) gsdg_angle(ap);
            if( use_noel) gsdg_noel(ap);
            /*	not-so abortive attempt to add contacts (oh well) */
            if( use_nonbon){
                for( j=0; j< numatom; j++)
                {
                    bp = a_next(j);
                    if( bp != ap){
                        if( bp->vx == 0.)
                        {
                            /*			if( bp->active || bp->x != 0. || bp->y != 0. || bp->z != 0.)
                            			{*/
                            x = ap->x -bp->x;
                            y = ap->y -bp->y;
                            z = ap->z -bp->z;
                            r = x*x + y*y + z*z;
                            if( r < 8.){ /* was 8 */
                                bp->vy = -10.;
                                bp->vx = 16.;
                            }
                            /*			} */
                        }
                    }
                }
            } /* if use_nonbon */
            inbps = 0;
            for( j=0; j< numatom ; j++)
            {
                bp = a_next(j);
                /*printf("%d %d  %f\n",ap->serial,bp->serial,bp->vx); */
                if( bp->vx > 0. && bp != ap
                        &&( bp->x != 0.|| bp->y != 0. || bp->z != 0. || !bp->active))
                { bps[inbps++] = bp; if( inbps == MAX_LOCAL) break;}
                else
                { /* it didn't make the cut so kill it */ bp->vx = -1.; bp->vy = -1.;}
            }
            if(echo) fprintf(op," %s %d inbps %d\n",ap->name,ap->serial,inbps);
            switch( inbps){
            case 0: break;/* nada no distances to use */
            case 1:
                rand3(&x,&y,&z);
                bp = bps[0]; r = sqrt(bp->vx);
                ap->x = bp->x + r*x;
                ap->y = bp->y + r*y;
                ap->z = bp->z + r*z;
                break;
            case 2: /* numerical solution is simplest because not unique */
            case 3:
                if( ap->x == 0. && ap->y == 0. && ap->z == 0.)
                {
                    ap->z = bps[0]->z;
                    ap->x = bps[1]->x;
                    ap->y = bps[1]->y - bps[0]->y;
                }
                /* cut from gsdg */
                for( j=0; j< 10; j++){
                    rvec[0] = 0;
                    rvec[1] = 0;
                    rvec[2] = 0;
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
#ifdef WINDOZ
                    BE_NICE();
#endif
                    x = abuild_line_search( rvec,&y,ap,bps,inbps);

                    ap->x += y*rvec[0];
                    ap->y += y*rvec[1];
                    ap->z += y*rvec[2];
                }
                break;
                /*
                		a = bps[0]->vx;
                		x = bps[0]->x;
                		y = bps[0]->y;
                		z = bps[0]->z;
                		a -= x*x + y*y + z*z;
                		b = bps[1]->vx;
                		x = bps[1]->x;
                		y = bps[1]->y;
                		z = bps[1]->z;
                		b -= x*x + y*y + z*z;
                		c = bps[2]->vx;
                		x = bps[2]->x;
                		y = bps[2]->y;
                		z = bps[2]->z;
                		c -= x*x + y*y + z*z;
                		for( j=0; j< 81; j++) m[j] = 0.;
                		for( j=3; j< 9; j++) v[j] = 0;
                		v[0] = a;
                		v[1] = b-a;
                		v[2] = c-a;
                		m[0] = 1.;
                		m[1] = 1.;
                		m[2] = 1.;
                		m[6] = -2.*bps[0]->x;
                		m[7] = -2.*bps[0]->y;
                		m[8] = -2.*bps[0]->z;
                		m[9+6] = -2.*(bps[1]->x - bps[0]->x);
                		m[9+7] = -2.*(bps[1]->y - bps[0]->y);
                		m[9+8] = -2.*(bps[1]->z - bps[0]->z);
                		m[18+6] = -2.*(bps[2]->x - bps[0]->x);
                		m[18+7] = -2.*(bps[2]->y - bps[0]->y);
                		m[18+8] = -2.*(bps[2]->z - bps[0]->z);
                		m[27] = m[9+6];
                		m[27+3] = m[9+7];
                		m[27+4] = m[9+8];
                		m[27+6] = -v[1];
                		m[36+1] = m[9+7];
                		m[36+3] = m[9+6];
                		m[36+5] = m[9+8];
                		m[36+7] = -v[1];
                		m[45+2] = m[9+8];
                		m[45+4] = m[9+6];
                		m[45+5] = m[9+7];
                		m[45+8] = -v[1];
                		m[54] = m[18+6];
                		m[54+3] = m[18+7];
                		m[54+4] = m[18+8];
                		m[54+6] = -v[2];
                		m[63+1] = m[18+7];
                		m[63+3] = m[18+6];
                		m[63+5] = m[18+8];
                		m[63+7] = -v[2];
                		m[72+2] = m[18+8];
                		m[72+4] = m[18+6];
                		m[72+5] = m[18+7];
                		m[72+8] = -v[2];

                		mom_solve(m ,v,9 ,9);
                	
                		for( j=0; j< 9; j++)
                			printf("%d %f\n",j,v[j]);
                		ap->x = v[6] ;
                		ap->y = v[7] ;
                		ap->z = v[8] ;
                //		gsdg(lp,ilocal,10,ap->serial,ap->serial);
                		kohonen_minimizer( lp,lf,ilocal,10);
                		break;
                */		
                /*		case 3:*/
            case 4:
            default:
                a = bps[0]->vx;
                x = bps[0]->x;
                y = bps[0]->y;
                z = bps[0]->z;
                a -= x*x + y*y + z*z;
                for( i=1; i< inbps; i++)
                {
                    d = bps[i]->vx;
                    x = bps[i]->x;
                    y = bps[i]->y;
                    z = bps[i]->z;
                    d -= x*x + y*y + z*z;
                    v[i-1] = d-a;
                    honker[i-1][0] = -2.*(bps[i]->x - bps[0]->x);
                    honker[i-1][1] = -2.*(bps[i]->y - bps[0]->y);
                    honker[i-1][2] = -2.*(bps[i]->z - bps[0]->z);
                }
                for( i=0; i< 9; i++)
                    m2[i] = 0.;
                for( i=0; i< 3; i++)
                    v2[i] = 0.;
                for( i=0; i<3; i++)
                    for(j=0; j< 3; j++)
                    {
                        for( k=0; k<inbps-1; k++)
                            m2[i*3 +j] += honker[k][i]*honker[k][j];
                    }
                for( i=0; i< 3; i++)
                    for( k=0; k< inbps-1; k++)
                        v2[i] += honker[k][i]*v[k];

                jacobi( m2,m,3,100,1.e-10);
                if( echo)fprintf(op," eig^2 %f %f %f\n", m2[0],m2[4],m2[8]);
                if( m2[0] > 0.1) m2[0] = 1./m2[0]; else m2[0] = 0.;
                if( m2[4] > 0.1) m2[4] = 1./m2[4]; else m2[4] = 0.;
                if( m2[8] > 0.1) m2[8] = 1./m2[8]; else m2[8] = 0.;
                a = m2[0]; b= m2[4]; c = m2[8];
                for( i=0; i< 9; i++)
                {m2[i] = 0.; m3[i] = 0.;}

                for(i=0; i<3; i++)
                    for( j=0; j<3; j++)
                        m2[i] += m[j*3+i]*v2[j];
                m2[0]*=a; m2[1] *=b; m2[2] *=c;
                for( i=0; i< 3; i++)
                    for(j=0; j<3; j++)
                        m3[i] += m[i*3+j]*m2[j];

                ap->x = m3[0];
                ap->y = m3[1];
                ap->z = m3[2];
                /* cut from gsdg */
                if( use_hybrid) gsdg_hybrid(ap);
                for( j=0; j< 10; j++){
                    rvec[0] = 0;
                    rvec[1] = 0;
                    rvec[2] = 0;
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
                    rand3( &svec[0],&svec[1],&svec[2]);

                    x = abuild_line_search( svec, &y,ap,bps,inbps);
                    rvec[0] += y*svec[0];
                    rvec[1] += y*svec[1];
                    rvec[2] += y*svec[2];
#ifdef WINDOZ
                    BE_NICE();
#endif
                    x = abuild_line_search( rvec,&y,ap,bps,inbps);

                    ap->x += y*rvec[0];
                    ap->y += y*rvec[1];
                    ap->z += y*rvec[2];
                }


                break;
            } /* end of switch*/


SKIP: ;
            ap = ap->next;
        }/* iatom */
#ifdef WINDOZ
        force_screen_update();
#endif
#ifdef GTK_VERSION
        AMMP_whiz_bang();
#endif

    }/* icycle */

}/* end of abuild */

float abuild_line_search( vect, step,who ,bps,inbps)
float vect[3],*step;
ATOM *who;
ATOM *bps[];
int inbps;
{
    float val;
    float vt,lam;
    int i,j;
    float dstep;

    float abuild_dgeom();

    val = abuild_dgeom(vect,0.,who,bps,inbps);
    lam = 0;
    *step = 0;
    dstep = -.5;
    for( i=0; i< 3; i++)
    {
        dstep *= -.5;
        for( j = 0; j< 200 ; j++)
        {
            lam += dstep;
            vt =  abuild_dgeom(vect,lam,who, bps,inbps);
            if( vt < val ){ *step = lam; val = vt;} else {break;}
        }
        if( j == 200) dstep *= -2;
    }
    return val;
}/*end of routine */

float abuild_dgeom( vect,lam,who, bps,inbps)
ATOM *who;
float vect[3],lam;
ATOM *bps[];
int inbps;
{
    int numatm,a_number();
    int i;
    float x,y,z;
    ATOM *ap,*a_next();
    float dt;
    float dsum;

    numatm = a_number();
    x = who->x + vect[0]*lam;
    y = who->y + vect[1]*lam;
    z = who->z + vect[2]*lam;

    dsum = 0.;

    for( i=0; i< inbps; i++)
    {
        ap = bps[i];
        if( ap != who )
        {
            if(ap->vx > 0.){
                if( ap->vy > 0 )
                {
                    dt = (x -ap->x)*(x-ap->x);
                    dt += (y -ap->y)*(y-ap->y);
                    dt += (z -ap->z)*(z-ap->z);

                    dsum += ap->vy*(ap->vx -dt)*(ap->vx -dt);
                } else {

                    dt = (x -ap->x)*(x-ap->x);
                    dt += (y -ap->y)*(y-ap->y);
                    dt += (z -ap->z)*(z-ap->z);

                    if( ap->vx > dt)
                        dsum -= ap->vy*(ap->vx -dt)*(ap->vx -dt);
                }
            }
        }
    }
    return dsum;
}/* end of routine */
