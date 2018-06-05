/* $Id: spm_update_Upsilon.c 7172 2017-09-21 16:31:30Z john $ */
/* (c) Mikael Brudfors (2018) */

#include "mex.h"
#include <math.h>
#define MAXCLASSES 1024

static void update_Upsilon(mwSize dm[], float R[], float w[], int ix_k, int ix_l, double val[])
{
    mwSize i0, i1, i2, m, n;
    float a = 0;
    float *R0 = NULL, *R1 = NULL;
    int it = 0;        
    double val1 = 0;
    
    m = dm[0]*dm[1]*dm[2];        
    
    mwSize i2start = it%2;
    for(i2=0; i2<dm[2]; i2++) /* Inferior -> Superior */
    {
        mwSize i1start = (i2start == (i2%2));
        for(i1=0; i1<dm[1]; i1++) /* Posterior -> Anterior */
        {
            mwSize i0start = (i1start == (i1%2));                
            R1 = R + dm[0]*(i1+dm[1]*i2);

            for(i0=i0start; i0<dm[0]; i0+=2) /* Left -> Right */
            {                    
                float *RR = NULL;

                /* Pointers to current voxel in first volume */                    
                R0 = R1 + i0;

                /* Initialise neighbour counts to zero */
                a = 0;

                /* Count neighbours of each class */
                if(i2>0)       /* Inferior */
                {

                    RR = R0 - dm[0]*dm[1];
                    a  += RR[ix_l*m]*w[2];
                }

                if(i2<dm[2]-1) /* Superior */
                {
                    RR = R0 + dm[0]*dm[1];
                    a  += RR[ix_l*m]*w[2];
                }

                if(i1>0)       /* Posterior */
                {
                    RR = R0 - dm[0];
                    a  += RR[ix_l*m]*w[1];
                }

                if(i1<dm[1]-1) /* Anterior */
                {
                    RR = R0 + dm[0];
                    a  += RR[ix_l*m]*w[1];
                }

                if(i0>0)       /* Left */
                {
                    RR = R0 - 1;
                    a  += RR[ix_l*m]*w[0];
                }

                if(i0<dm[0]-1) /* Right */
                {
                    RR = R0 + 1;
                    a  += RR[ix_l*m]*w[0];
                }
                                  
                val1 += (double)R0[ix_k*m]*a;

            } /* < Loop over dm[0]*/
        } /* < Loop over dm[1]*/
    } /* < Loop over dm[2]*/
    
    val[0] = val1;
    
} /* update_Upsilon */
        
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Initialise variables */
    mwSize i;
    mwSize dm[4];
    float w[3];
    float *R = NULL;
    
    /* Check correct number of inputs and outputs */
    if (nrhs!=4 || nlhs!=1)
        mexErrMsgTxt("Incorrect usage");

    /* Check correct data-type of inputs */
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("First arg must be numeric, real, full and single.");

    /* Get dimensions */
    if (mxGetNumberOfDimensions(prhs[0])>4)
        mexErrMsgTxt("First arg has wrong number of dimensions.");

    for(i=0; i<mxGetNumberOfDimensions(prhs[0]); i++)
        dm[i] = mxGetDimensions(prhs[0])[i];

    for(i=mxGetNumberOfDimensions(prhs[0]); i<4; i++)
        dm[i] = 1;

    if (dm[3]>MAXCLASSES) mexErrMsgTxt("Too many classes.");
    
    /* Adjustment for anisotropic voxel sizes.  w should contain
       the square of each voxel size. */
    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsSingle(prhs[1]))
        mexErrMsgTxt("Third arg must be numeric, real, full and single.");

    if (mxGetNumberOfElements(prhs[1]) != 3)
        mexErrMsgTxt("Third arg must contain three elements.");

    /* Read input */
    R = (float *)mxGetData(prhs[0]);
        
    for(i=0; i<3; i++) 
        w[i] = ((float *)mxGetData(prhs[1]))[i];

    int k = mxGetScalar(prhs[2]);
    k     = k - 1;
    
    int l = mxGetScalar(prhs[3]);
    l     = l - 1;
    
    /* Allocate output */
    mwSize dm1[2];
    for(i=0; i<2; i++)
        dm1[i] = 1;
    
    double *val;
    plhs[0] = mxCreateNumericArray(2,dm1, mxDOUBLE_CLASS, mxREAL);
    val     = (double *)mxGetData(plhs[0]);

    /* mexPrintf("Enter::update_Upsilon\n"); */
    
    update_Upsilon(dm,R,w,k,l,val);
}