/* 
[Qt,c] = mexTriDTriD(sDim,nDim,randSeed); 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string> 
#include <iostream>
#include <vector>
#include <cmath>
#include "mex.h"

using namespace std;

/* gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{    
    /* [A,c,Ks] = mexAggregateSDPcone(A0,c0,K0s,Ks) */
    /* ************************************************** */
    
        /* Check for proper number of input and output arguments */    
    if (nrhs != 3) {
        mexErrMsgTxt("3 input arguments required.");
    } 
//    if(nlhs != 2){
//        mexErrMsgTxt("2 output arguments required");
//    } 

    
    /* ---> */ 
    /* reading sDim */
    const mxArray* sDim_ptr = prhs[0];
    unsigned int sDim = 0;  
    sDim = (unsigned int) *mxGetPr(sDim_ptr);    
    std::cout << "sDim = " << sDim << std::endl;
            
    const mxArray* nDim_ptr = prhs[1];
    unsigned int nDim = 0;  
    nDim = (unsigned int) *mxGetPr(nDim_ptr);
    std::cout << "nDim = " << nDim << std::endl;

    const mxArray* randSeed_ptr = prhs[2];
    unsigned int randSeed = 0;  
    randSeed = (unsigned int) *mxGetPr(randSeed_ptr);
    std::cout << "randSeed = " << randSeed << std::endl;
    
    srand(randSeed);
    
    mwSize Qtrowsize = (1+sDim)*(1+sDim); 
    mwSize QtcolSize = nDim*nDim;
    mwIndex Qtnnz = (3*nDim-2)*(3*sDim+1); 
    
    plhs[0] = mxCreateSparse(Qtrowsize,QtcolSize,Qtnnz,mxREAL);
    mwIndex* Qtjc = mxGetJc(plhs[0]);
    mwIndex* Qtir = mxGetIr(plhs[0]);
    double*  Qtpr = mxGetPr(plhs[0]);

    vector<double> aVect(1+sDim); 
    vector<double> dVect(1+sDim); 
    
    mwIndex nnzPt = 0; 
    mwIndex cloPt = 0; 
    Qtjc[0] = 0; 
    
    for (mwIndex p=0; p < nDim; p++) {
        mwIndex qStart = p;
        if (p > 0) qStart = p-1;
        mwIndex qEnd = p+1;
        if (p == nDim-1) qEnd = nDim-1;
        for (mwIndex q=qStart; q <=qEnd; q++) {
            if (q >= p) {
                for (mwIndex i = 0; i <= sDim; i++) {
                    aVect[i] = (double) rand()/(1.0+RAND_MAX);
                    // std::cout << "aVect[" << i << "] = " << aVect[i] << std::endl;
                    dVect[i] = (double) rand()/(1.0+RAND_MAX);
                }
            } 
            // j = 0; 
            for (mwIndex i=0; i <= sDim; i++) {
                Qtir[nnzPt] = i;
                Qtpr[nnzPt] = aVect[i];
                nnzPt++;
            }
            for (mwIndex j=1; j <= sDim; j++) {
                Qtir[nnzPt] = j*(sDim+1);
                Qtpr[nnzPt] = aVect[j];
                nnzPt++;
                Qtir[nnzPt] = j*(sDim+1)+j;
                Qtpr[nnzPt] = dVect[j];
                nnzPt++;
            }
        }
        Qtjc[colPt+1] = Qtjc[colPtPt] + 3*sDim+1; 
        colPtPt++; 
        if (p < nDim -1) {
            for (mwIndex k = 0; k < nDim -2; k++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        }   
    }
    
    return;
}
