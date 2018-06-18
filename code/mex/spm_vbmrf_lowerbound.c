/* $Id: spm_vbmrf_lowerbound.c 7172 2017-09-21 16:31:30Z john $ */
/* (c) Mikael Brudfors (2018) */

#include "mex.h"
#include <math.h>
#define MAXCLASSES 1024

static void compute_lowebound(mwSize dm[], float R[], float ln_upsilon[], float w[], int code, double val[])
{
    mwSize i0, i1, i2, k, m, n;
    float a[MAXCLASSES], e[MAXCLASSES];
    float *R0 = NULL, *R1 = NULL;
    int it = 0;
    double val1 = 0;
    int cnt_nb = 0;
    
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
                for(k=0; k<dm[3]; k++) a[k] = 0.0;

                /* Count neighbours of each class */
                cnt_nb = 0;
                if(i2>0)       /* Inferior */
                {

                    RR = R0 - dm[0]*dm[1];
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[2];
                            cnt_nb++;
                        }
                    }
                }

                if(i2<dm[2]-1) /* Superior */
                {
                    RR = R0 + dm[0]*dm[1];
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[2];
                            cnt_nb++;
                        }
                    }
                }

                if(i1>0)       /* Posterior */
                {
                    RR = R0 - dm[0];
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[1];
                            cnt_nb++;
                        }
                    }
                }

                if(i1<dm[1]-1) /* Anterior */
                {
                    RR = R0 + dm[0];
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[1];
                            cnt_nb++;
                        }
                    }
                }

                if(i0>0)       /* Left */
                {
                    RR = R0 - 1;
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[0];
                            cnt_nb++;
                        }
                    }
                }

                if(i0<dm[0]-1) /* Right */
                {
                    RR = R0 + 1;
                    for(k=0; k<dm[3]; k++) 
                    {
                        if ((RR[k*m]!=RR[k*m]) == false && RR[k*m]!=0) /* check if NaN */
                        {   
                            a[k] += RR[k*m]*w[0];
                            cnt_nb++;
                        }
                    }
                }

                /* Responsibility data is uint8, so correct scaling.
                   Note also that data is divided by 6 (the number
                   of neighbours examined). */
                for(k=0; k<dm[3]; k++)
                    a[k]/=(255.0*((float)cnt_nb));

                if (code == 1) 
                {
                    /* Weights are in the form of a matrix,
                       shared among all voxels. */
                    float *g;
                    for(k=0, g=ln_upsilon; k<dm[3]; k++)
                    {
                        e[k] = 0;
                        for(n=0; n<dm[3]; n++, g++)
                            e[k] += (*g)*a[n];
                    }
                }
                else if (code == 2)
                {
                    /* Weights are assumed to be a diagonal matrix,
                       so only the diagonal elements are passed. */
                    for(k=0; k<dm[3]; k++)
                        e[k] = ln_upsilon[k]*a[k];
                }
                else if (code == 3)
                {
                    /* Separate weights for each voxel, in the form of
                       the full matrix (loads of memory). */
                    float *g;
                    g = ln_upsilon + i0+dm[0]*(i1+dm[1]*i2);
                    for(k=0; k<dm[3]; k++)
                    {
                        e[k] = 0.0;
                        for(n=0; n<dm[3]; n++, g+=m)
                            e[k] += (*g)*a[n];
                    }
                }
                else if (code == 4)
                {
                    /* Separate weight matrices for each voxel,
                       where the matrices are assumed to be symmetric
                       with zeros on the diagonal. For a 4x4
                       matrix, the elements are ordered as
                       (2,1), (3,1), (4,1), (3,2), (4,2), (4,3).
                     */
                    float *g;
                    g = ln_upsilon + i0+dm[0]*(i1+dm[1]*i2);
                    for(k=0; k<dm[3]; k++) e[k] = 0.0;

                    for(k=0; k<dm[3]; k++)
                    {
                        for(n=k+1; n<dm[3]; n++, g+=m)
                        {
                            e[k] += (*g)*a[n];
                            e[n] += (*g)*a[k];
                        }
                    }
                }

                for(k=0; k<dm[3]; k++)
                {
                    if ((R0[k*m]!=R0[k*m]) == false && R0[k*m]!=0 && cnt_nb!=0) /* check if NaN */
                    {
                        val1 += (double)R0[k*m]*e[k];                        
                    }
                }

            } /* < Loop over dm[0]*/
        } /* < Loop over dm[1]*/
    } /* < Loop over dm[2]*/
    
    val[0] = val1;
    
} /* compute_lowebound */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Initialise variables */
    mwSize i;
    mwSize dm[4];
    float w[3];
    float *ln_upsilon = NULL;
    float *R = NULL;
    int code = 0;

    /* Check correct number of inputs and outputs */
    if (nrhs!=3 || nlhs!=1)
        mexErrMsgTxt("Incorrect usage");

    /* Check correct data-type of inputs */
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsSingle(prhs[0]))
        mexErrMsgTxt("First arg must be numeric, real, full and single.");

    if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]))
        mexErrMsgTxt("Second arg must be numeric, real and full.");

    /* Get dimensions */
    if (mxGetNumberOfDimensions(prhs[0])>4)
        mexErrMsgTxt("First arg has wrong number of dimensions.");

    for(i=0; i<mxGetNumberOfDimensions(prhs[0]); i++)
        dm[i] = mxGetDimensions(prhs[0])[i];

    for(i=mxGetNumberOfDimensions(prhs[0]); i<4; i++)
        dm[i] = 1;

    if (dm[3]>MAXCLASSES) mexErrMsgTxt("Too many classes.");

    /* Set code, which indicates type of weight array (ln_upsilon) */
    if (mxGetDimensions(prhs[1])[1] == 1)
    {
        code = 2;
        if (mxGetDimensions(prhs[1])[0] != dm[3])
            mexErrMsgTxt("Second arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[1])) mexErrMsgTxt("Second arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[1])==2)
    {
        code = 1;
        if (mxGetDimensions(prhs[1])[0] != dm[3] || mxGetDimensions(prhs[1])[1] != dm[3])
            mexErrMsgTxt("Second arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[1])) mexErrMsgTxt("Second arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[1])==5)
    {
        code = 3;
        for(i=0; i<4; i++)
            if (mxGetDimensions(prhs[1])[i] != dm[i])
                mexErrMsgTxt("Second arg has incompatible dimensions.");

        if (mxGetDimensions(prhs[1])[4] != dm[3])
            mexErrMsgTxt("Second arg has incompatible dimensions.");

        if (!mxIsSingle(prhs[1])) mexErrMsgTxt("Second arg must be single.");
    }
    else if (mxGetNumberOfDimensions(prhs[1])==4)
    {
        for(i=0; i<3; i++)
            if (mxGetDimensions(prhs[1])[i] != dm[i])
                mexErrMsgTxt("Second arg has incompatible dimensions.");

        if (mxGetDimensions(prhs[1])[3] != (dm[3]*(dm[3]-1))/2)
            mexErrMsgTxt("Second arg has incompatible dimensions.");

        code = 4;
    }
    else
        mexErrMsgTxt("Second arg has incompatible dimensions.");

    /* Get data */
    ln_upsilon = (float *)mxGetData(prhs[1]);
    
    /* Adjustment for anisotropic voxel sizes.  w should contain
       the square of each voxel size. */
    if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsSingle(prhs[2]))
        mexErrMsgTxt("Third arg must be numeric, real, full and single.");

    if (mxGetNumberOfElements(prhs[2]) != 3)
        mexErrMsgTxt("Third arg must contain three elements.");

    for(i=0; i<3; i++) w[i] = ((float *)mxGetData(prhs[2]))[i];

    /* Copy input to output */
    R = (float *)mxGetData(prhs[0]);    
    
    mwSize dm1[2];
    for(i=0; i<2; i++)
        dm1[i] = 1;
        
    double *val;
    plhs[0] = mxCreateNumericArray(2,dm1, mxDOUBLE_CLASS, mxREAL);
    val     = (double *)mxGetData(plhs[0]);
    
    compute_lowebound(dm,R,ln_upsilon,w,code,val);
}