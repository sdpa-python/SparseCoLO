/* 
[Qt,c, sDim, Js] = mexDiagTriDQOP(sDim,nDim,randSeed); 
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
    // std::cout << "sDim = " << sDim << std::endl;
            
    const mxArray* nDim_ptr = prhs[1];
    unsigned int nDim = 0;  
    nDim = (unsigned int) *mxGetPr(nDim_ptr);
    // std::cout << "nDim = " << nDim << std::endl;

    const mxArray* randSeed_ptr = prhs[2];
    unsigned int randSeed = 0;  
    randSeed = (unsigned int) *mxGetPr(randSeed_ptr);
    // std::cout << "randSeed = " << randSeed << std::endl;
    
    // Constructing output Qt
    
    mwSize Qtrowsize = (1+sDim)*(1+sDim); 
    // std::cout << "Qtrowsize = " << Qtrowsize << std::endl;
    
    mwSize QtcolSize = nDim*nDim;
    // std::cout << "QtcolSize = " << QtcolSize << std::endl;

    mwIndex Qtnnz = (5*nDim-4)*sDim+nDim; 
    // std::cout << "Qtnnz     = " << Qtnnz << std::endl;

    plhs[0] = mxCreateSparse(Qtrowsize,QtcolSize,Qtnnz,mxREAL);
    mwIndex* Qtjc = mxGetJc(plhs[0]);
    mwIndex* Qtir = mxGetIr(plhs[0]);
    double*  Qtpr = mxGetPr(plhs[0]);

    vector<double> Lvect(1+sDim); 
    vector<double> offDVect(1+sDim); 
    vector<double> dVect(1+sDim); 
    
    mwIndex nnzPt = 0; 
    mwIndex colPt = 0; 
    Qtjc[0] = 0; 
    
    srand(randSeed);

    for (mwIndex p=0; p < nDim; p++) {
        mwIndex qStart = p;
        if (p > 0) qStart = p-1;
        mwIndex qEnd = p+1;
        if (p == nDim-1) qEnd = nDim-1;
        for (mwIndex q=qStart; q <=qEnd; q++) {
            if (q >= p) {
                for (mwIndex i = 0; i <= sDim; i++) {
                    Lvect[i] = (double) rand()/(1.0+RAND_MAX)-0.5;
                    dVect[i] = (double) (-rand()/(1.0+RAND_MAX));
                    dVect[0] = 1.0; 
                    // std::cout << "   dVect[" << i << "] = " << dVect[i] << std::endl;
                }
            } 
            if (q == p) {
                for (mwIndex j=0; j <= sDim; j++) {
                    Qtir[nnzPt] = (sDim+2)*j;
                    Qtpr[nnzPt] = dVect[j];
                    nnzPt++;                             
                }
                Qtjc[colPt+1] = Qtjc[colPt] + (sDim+1); 
                colPt++;                             
            }  
            else {
                for (mwIndex i=1; i<=sDim; i++) {
                    Qtir[nnzPt] = i;
                    Qtpr[nnzPt] = Lvect[i];
                    nnzPt++;               
                }
                for (mwIndex j=1; j <= sDim; j++) {
                    Qtir[nnzPt] = (sDim+1)*j;
                    Qtpr[nnzPt] = Lvect[j];
                    nnzPt++;                             
                }
                Qtjc[colPt+1] = Qtjc[colPt] + 2*sDim; 
                colPt++;             
            }
        }
        if (p < nDim-1) {
            for (mwIndex k = 0; k < nDim - 2; k++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        }   
    }
    
    for (mwIndex colPt = 0; colPt <= QtcolSize; colPt++)
        // std::cout << "Qtjc[" << colPt << "] = " << Qtjc[colPt] << std::endl;
    for (mwIndex nnzPt = 0; nnzPt < Qtnnz; nnzPt++) {
        // std::cout << "Qtir[" << nnzPt << "] = " << Qtir[nnzPt] << std::endl;
        // std::cout << "Qtpr[" << nnzPt << "] = " << Qtpr[nnzPt] << std::endl;        
    }    
    
    // Constructing output ct
    
    mwSize crowsize = (sDim+1)*(sDim+1); 
    mwSize ccolSize = 1; 
    mwSize cnnz = 2*sDim; 
    
    plhs[1] = mxCreateSparse(crowsize,ccolSize,cnnz,mxREAL);
    mwIndex* cjc = mxGetJc(plhs[1]);
    mwIndex* cir = mxGetIr(plhs[1]);
    double*  cpr = mxGetPr(plhs[1]);

    for (mwIndex i = 0; i <= sDim; i++) {
        Lvect[i] = (double) rand()/(1.0+RAND_MAX)-0.5;
    }
    cjc[0] = 0;
    nnzPt = 0; 
    for (mwIndex i=1; i<=sDim; i++) {
        cir[nnzPt] = i;
        cpr[nnzPt] = Lvect[i];
        nnzPt++;               
    }
    for (mwIndex j=1; j <= sDim; j++) {
        cir[nnzPt] = (sDim+1)*j;
        cpr[nnzPt] = Lvect[j];
        nnzPt++;                             
    }
    cjc[1] = cnnz;    
    
    // Constructing output sDim
    plhs[2]= mxCreateDoubleScalar(sDim); 

    // Constructing output nDim
    plhs[3]= mxCreateDoubleScalar(nDim); 

    return;
}
