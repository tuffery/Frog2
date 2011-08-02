/* routine to superimpose two sets of coordinates
*  
* cadged from Larry Andrews' Newvector.f library
*
C  The algorithm is from Ferro and Hermans, Acta Cryst., A33,345-347
C  (1977).  Other algorithms give the same result, such as:
C
C  Kearsley, Acta Cryst., A45, 208-210 (1989)
C  Diamond,  Acta Cryst., A44, 211-216 (1988)
C  Kabsch,   Acta Cryst., A34, 827-828 (1978)
C  Kabsch,   Acta Cryst., A32, 922-923 (1976)
C
C  The algorithm of Ferro and Hermans has the advantage that the
C  rotational components are removed iteratively, so that residual
C  numerical errors are correctly removed.  This leaves no
C  possibility of poorly orthogonalized matrices.  The explicit
C  iteration is not a problem.  The other methods simply hide the
C  iteration in the eigenvector solver.  Ferro and Hermans
C  simply build what is a simple, special purpose version of the QR
C  method.
*/

#include <math.h>
/*void matmul(float[],float[],float[],int,int);
*/

/* return the RMS error, -1 on internal error */
float bstrot( x,y,z,xx,yy,zz, na, matrix,delta)
int na;
float x[],y[],z[];
float xx[],yy[],zz[];
float matrix[3][3], delta[3];
{

    float tensor[3][3],cx,cy,cz,cxx,cyy,czz;
    float tx,ty,tz,txx,tyy,tzz;
    float rms;
    float sx[3][3],sy[3][3],sz[3][3];
    float sq[3][3];
    void cpyvec();
    void matmul();
    int i,j,ipass;


    if( na < 1) return -1.;
    /* find the centers of mass */
    cx = 0.;
    cy = 0.;
    cz = 0.;
    cxx = 0.;
    cyy = 0.;
    czz = 0.;
    for( i=0; i< na; i++)
    {
        cx += x[i];
        cy += y[i];
        cz += z[i];
        cxx += xx[i];
        cyy += yy[i];
        czz += zz[i];
    }
    cx /= na;
    cy /= na;
    cz /= na;
    cxx /= na;
    cyy /= na;
    czz /= na;
    /* make the metric tensor */
    for( i=0; i< 3; i++)
        for( j=0; j< 3; j++)
        {
            tensor[i][j] = 0.;
            matrix[i][j] = 0.;
            sx[i][j] = 0.;
            sy[i][j] = 0.;
            sz[i][j] = 0.;
        }
    matrix[0][0] = 1.;
    matrix[1][1] = 1.;
    matrix[2][2] = 1.;
    sx[0][0] = 1.;
    sy[1][1] = 1.;
    sz[2][2] = 1.;
    for( i=0; i<na; i++)
    {
        tx = x[i] - cx;
        ty = y[i] - cy;
        tz = z[i] - cz;
        txx = xx[i] - cxx;
        tyy = yy[i] - cyy;
        tzz = zz[i] - czz;
        tensor[0][0] += tx*txx;
        tensor[0][1] += tx*tyy;
        tensor[0][2] += tx*tzz;
        tensor[1][0] += ty*txx;
        tensor[1][1] += ty*tyy;
        tensor[1][2] += ty*tzz;
        tensor[2][0] += tz*txx;
        tensor[2][1] += tz*tyy;
        tensor[2][2] += tz*tzz;
    }
    /* now find the linear orthogonal transformation which symetrizes
    * the metric tensor */

    for( ipass = 0; ipass < 20; ipass ++)
    {
        rms = 0.;
        /* x */
        tx = atan2( tensor[2][1]-tensor[1][2],
                    tensor[1][1]+tensor[2][2]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sx[1][1] = ty;
        sx[2][1] = -tz;
        sx[1][2] = tz;
        sx[2][2] = ty;
        matmul( sx,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sx,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /* y */
        tx = atan2( tensor[2][0]-tensor[0][2],
                    tensor[0][0]+tensor[2][2]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sy[0][0] = ty;
        sy[2][0] = -tz;
        sy[0][2] = tz;
        sy[2][2] = ty;
        matmul( sy,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sy,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /* z */
        tx = atan2( tensor[0][1]-tensor[1][0],
                    tensor[1][1]+tensor[0][0]);
        rms += fabs(tx);
        ty = cos(tx); tz = sin(tx);
        sz[1][1] = ty;
        sz[0][1] = -tz;
        sz[1][0] = tz;
        sz[0][0] = ty;
        matmul( sz,tensor,sq,3,3);
        cpyvec(sq,tensor,9);
        matmul( sz,matrix,sq,3,3);
        cpyvec( sq,matrix,9);
        /*  termination critereon here */
        if( rms < 1.e-7) break;
    }
    rms = 0.;
    for(i=0; i< na; i++)
    {
        txx = xx[i] - cxx;
        tyy = yy[i] - cyy;
        tzz = zz[i] - czz;
        tx = matrix[0][0]*txx + matrix[1][0]*tyy + matrix[2][0]*tzz;
        ty = matrix[0][1]*txx + matrix[1][1]*tyy + matrix[2][1]*tzz;
        tz = matrix[0][2]*txx + matrix[1][2]*tyy + matrix[2][2]*tzz;
        tx += cx - x[i];
        ty += cy - y[i];
        tz += cz - z[i];
        rms += tx*tx + ty*ty + tz*tz;
    }
    tx = matrix[0][0]*cxx + matrix[1][0]*cyy + matrix[2][0]*czz;
    ty = matrix[0][1]*cxx + matrix[1][1]*cyy + matrix[2][1]*czz;
    tz = matrix[0][2]*cxx + matrix[1][2]*cyy + matrix[2][2]*czz;
    delta[0] = cx - tx;
    delta[1] = cy - ty;
    delta[2] = cz - tz;
    return sqrt(rms/na);


}/* end of routine */
/* copy a vector into another */
void cpyvec(orig,copy,n)
float orig[],copy[];
int n;
{
    int i;
    for( i=0; i< n; i++)
        copy[i] = orig[i];
}
/*
multiply c[n][n] =  a[n][m] b[m][n];
*/
void matmul( a,b,c,n,m)
float a[],b[],c[];
int n,m;
{
    int i,j,k,ioff,koff;

    for( i=0; i< n*n; i++)
        c[i] = 0.;
    for( i=0; i< n; i++)
    {
        ioff = i*n;
        for( j=0; j< n; j++)
        {
            koff = 0.;
            for( k=0; k<m; k++)
            {
                c[ ioff +j] += a[ioff + k] *b[ j +koff];
                koff += m;
            }
        }
    }

}
