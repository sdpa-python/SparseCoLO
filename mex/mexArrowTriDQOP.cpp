/* 
[Qt,c, sDim, Js] = mexArrowTriDQOP(sDim,kBlock,randSeed); 
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
            
    const mxArray* kblock_ptr = prhs[1];
    unsigned int kblock = 0;  
    kblock = (unsigned int) *mxGetPr(kblock_ptr);
    // std::cout << "nDim = " << nDim << std::endl;

    const mxArray* randSeed_ptr = prhs[2];
    unsigned int randSeed = 0;  
    randSeed = (unsigned int) *mxGetPr(randSeed_ptr);
    // std::cout << "randSeed = " << randSeed << std::endl;
    
    // nDim
    unsigned int nDim = 2*kblock + 1;
    
    // Constructing output Qt
    
    mwSize Qtrowsize = (1+sDim)*(1+sDim); 
    // std::cout << "Qtrowsize = " << Qtrowsize << std::endl;
    
    mwSize QtcolSize = nDim*nDim;
//    std::cout << "QtcolSize = " << QtcolSize << std::endl;

    mwIndex Qtnnz = nDim*(sDim+1)+6*kblock*(5*sDim-2); 
//    std::cout << "Qtnnz     = " << Qtnnz << std::endl;

    plhs[0] = mxCreateSparse(Qtrowsize,QtcolSize,Qtnnz,mxREAL);
    mwIndex* Qtjc = mxGetJc(plhs[0]);
    mwIndex* Qtir = mxGetIr(plhs[0]);
    double*  Qtpr = mxGetPr(plhs[0]);

    vector<double> Lvect(1+sDim); 
    vector<double> offDVect(1+sDim); 
    vector<double> dVect(1+sDim); 
    vector<double> bLvect(1+sDim); 
    vector<double> boffDVect(1+sDim); 
    vector<double> bdVect(1+sDim); 

    mwIndex nnzPt = 0; 
    mwIndex colPt = 0; 
    Qtjc[0] = 0; 
    
    srand(randSeed);
    
    for (mwIndex k=0; k < kblock; k++) {
        if (k > 0) {
            for (mwIndex j = 0; j < 2*k; j++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        } 
        // the first (left upper) diagonal element of Qt
        Qtir[nnzPt] = 0;
        Qtpr[nnzPt] = 1.0;
        nnzPt++;                             
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+2)*j;
            Qtpr[nnzPt] = -rand()/(1.0+RAND_MAX);
            nnzPt++;                             
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(sDim+1); 
        colPt++;
        //
        for (mwIndex i = 0; i <= sDim; i++) {
            Lvect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            dVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            offDVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
                    // std::cout << "   dVect[" << i << "] = " << dVect[i] << std::endl;
        }
        Lvect[0] = 0.0;
        // the first (or left below) off-diagonal element
        for (mwIndex i=1; i <= sDim; i++){
            Qtir[nnzPt] = i;
            Qtpr[nnzPt] = Lvect[i];
            nnzPt++;
        }
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+1)*j; 
            Qtpr[nnzPt] = Lvect[j];
            nnzPt++;
            if (j > 1) {
                Qtir[nnzPt] = (sDim+2)*j-1; 
                Qtpr[nnzPt] = offDVect[j-1];
                nnzPt++;
            }
            Qtir[nnzPt] = (sDim+2)*j; 
            Qtpr[nnzPt] = dVect[j];
            nnzPt++;
            if (j < sDim) {
                Qtir[nnzPt] = (sDim+2)*j+1; 
                Qtpr[nnzPt] = offDVect[j];
                nnzPt++;
            }                                       
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(5*sDim-2); 
        colPt++;  
        //
        if (k < kblock-1) {
            for (mwIndex j = 0; j < 2*(kblock-k-1); j++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        }
        // the left borderd element 
        for (mwIndex i = 0; i <= sDim; i++) {
            bLvect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            bdVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            boffDVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            // std::cout << "   dVect[" << i << "] = " << dVect[i] << std::endl;
        }
        bLvect[0] = 0.0;
        for (mwIndex i=1; i <= sDim; i++){
            Qtir[nnzPt] = i;
            Qtpr[nnzPt] = bLvect[i];
            nnzPt++;
        }
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+1)*j; 
            Qtpr[nnzPt] = bLvect[j];
            nnzPt++;
            if (j > 1) {
                Qtir[nnzPt] = (sDim+2)*j-1; 
                Qtpr[nnzPt] = boffDVect[j-1];
                nnzPt++;
            }
            Qtir[nnzPt] = (sDim+2)*j; 
            Qtpr[nnzPt] = bdVect[j];
            nnzPt++;
            if (j < sDim) {
                Qtir[nnzPt] = (sDim+2)*j+1; 
                Qtpr[nnzPt] = boffDVect[j];
                nnzPt++;
            }                                       
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(5*sDim-2); 
        colPt++;  
        // 
        if (k > 0) {
            for (mwIndex j = 0; j < 2*k; j++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        } 
        // the second (or right upper) off-diagonal element
        for (mwIndex i=1; i <= sDim; i++){
            Qtir[nnzPt] = i;
            Qtpr[nnzPt] = Lvect[i];
            nnzPt++;
        }
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+1)*j; 
            Qtpr[nnzPt] = Lvect[j];
            nnzPt++;
            if (j > 1) {
                Qtir[nnzPt] = (sDim+2)*j-1; 
                Qtpr[nnzPt] = offDVect[j-1];
                nnzPt++;
            }
            Qtir[nnzPt] = (sDim+2)*j; 
            Qtpr[nnzPt] = dVect[j];
            nnzPt++;
            if (j < sDim) {
                Qtir[nnzPt] = (sDim+2)*j+1; 
                Qtpr[nnzPt] = offDVect[j];
                nnzPt++;
            }                                       
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(5*sDim-2); 
        colPt++;  
        // the second (left lower) diagonal element of Qt
        Qtir[nnzPt] = 0;
        Qtpr[nnzPt] = 1.0;
        nnzPt++;                             
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+2)*j;
            Qtpr[nnzPt] = -rand()/(1.0+RAND_MAX);
            nnzPt++;                             
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(sDim+1); 
        colPt++;
        //
        if (k < kblock-1) {
            for (mwIndex j = 0; j < 2*(kblock-k-1); j++) {
                Qtjc[colPt+1] = Qtjc[colPt]; 
                colPt++; 
            }
        }
        // the right borderd element 
        for (mwIndex i = 0; i <= sDim; i++) {
            bLvect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            bdVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            boffDVect[i] = (double) 2*(rand()/(1.0+RAND_MAX)-0.5);
            // std::cout << "   dVect[" << i << "] = " << dVect[i] << std::endl;
        }
        bLvect[0] = 0.0;
        for (mwIndex i=1; i <= sDim; i++){
            Qtir[nnzPt] = i;
            Qtpr[nnzPt] = bLvect[i];
            nnzPt++;
        }
        for (mwIndex j=1; j <= sDim; j++) {
            Qtir[nnzPt] = (sDim+1)*j; 
            Qtpr[nnzPt] = bLvect[j];
            nnzPt++;
            if (j > 1) {
                Qtir[nnzPt] = (sDim+2)*j-1; 
                Qtpr[nnzPt] = boffDVect[j-1];
                nnzPt++;
            }
            Qtir[nnzPt] = (sDim+2)*j; 
            Qtpr[nnzPt] = bdVect[j];
            nnzPt++;
            if (j < sDim) {
                Qtir[nnzPt] = (sDim+2)*j+1; 
                Qtpr[nnzPt] = boffDVect[j];
                nnzPt++;
            }                                       
        }
        Qtjc[colPt+1] = Qtjc[colPt]+(5*sDim-2); 
        colPt++;          
     }
    
//     std::cout << "colPt = " << colPt << std::endl;
//     std::cout << "nnzPt = " << nnzPt << std::endl;
     
//     for (mwIndex colPt = 0; colPt <= QtcolSize; colPt++)
//         std::cout << "Qtjc[" << colPt << "] = " << Qtjc[colPt] << std::endl;
//     for (mwIndex nnzPt = 0; nnzPt < Qtnnz; nnzPt++) {
//         std::cout << "Qtir[" << nnzPt << "] = " << Qtir[nnzPt] << std::endl;
//         std::cout << "Qtpr[" << nnzPt << "] = " << Qtpr[nnzPt] << std::endl;        
//     }    
    
    // the last column of Qt        
    for (mwIndex q=nDim-1; q < nDim*(nDim-1); q += nDim) { 
        for (mwIndex i=Qtjc[q];  i < Qtjc[q+1]; i++) {
            Qtir[nnzPt] = Qtir[i];
            Qtpr[nnzPt] = Qtpr[i]; 
            nnzPt++;
        } 
        Qtjc[colPt+1] = Qtjc[colPt]+Qtjc[q+1] - Qtjc[q]; 
        colPt++; 
    }
        
    // the last (nDim,nDim) element of Qt
    Qtir[nnzPt] = 0;
    Qtpr[nnzPt] = 1.0;
    nnzPt++;                             
    for (mwIndex j=1; j <= sDim; j++) {
        Qtir[nnzPt] = (sDim+2)*j;
        Qtpr[nnzPt] = (double) -rand()/(1.0+RAND_MAX);
        nnzPt++;                             
    }
    Qtjc[colPt+1] = Qtjc[colPt]+(sDim+1); 
    
    
    // for (mwIndex colPt = 0; colPt <= QtcolSize; colPt++)
        // std::cout << "Qtjc[" << colPt << "] = " << Qtjc[colPt] << std::endl;
    // for (mwIndex nnzPt = 0; nnzPt < Qtnnz; nnzPt++) {
        // std::cout << "Qtir[" << nnzPt << "] = " << Qtir[nnzPt] << std::endl;
        // std::cout << "Qtpr[" << nnzPt << "] = " << Qtpr[nnzPt] << std::endl;        
    // }    
    
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
