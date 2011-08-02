/* dscf_line is called by dscf_krymin()
*
* it searches for the closest local minimum
*
*/
float dscf_line( search,where,params,n, what )
int n;
float search[],params[],*where[];
float (*what)();
{
    float step,dstep,myval,bestval;
#define MAXBUF 100
    float steps[MAXBUF],vals[MAXBUF];
    int inbuf;
    int i,irun,itry;;

    inbuf = 1;
    steps[0] = 0.;
    vals[0] = (*what)();
    bestval = vals[0];

    step = 0.;
    dstep = 1.;
    for( itry = 0; itry < 8; itry ++)
    {
        /* use irun to search  in terms of dstep */
        for( irun = 0; irun < 100; irun ++)
        {
            /* first look in the table */

            for( i=0; i< inbuf; i++)
            {
                if( dstep == steps[i] ) {
                    myval = vals[i];
                    goto STEP_FOUND;
                }
            }
            /* if here then we have to calculate the function */
            for( i=0; i< n; i++)
            {
                *where[i] = params[i] + search[i]*dstep;
            }
            myval = (*what)();
            if( inbuf < MAXBUF)
            {vals[inbuf] = myval;steps[inbuf++] = dstep;}
STEP_FOUND:
            if( myval > bestval)
            {
                dstep = -dstep*.5;
                break;
            } else {
                bestval = myval;
                step = dstep;
            }

        }/* irun */
    }/* itry */

    /* update to the best found */
    for( i=0; i< n; i++)
    {
        *where[i] = params[i] + search[i]*step;
    }
    return step;
}

