#include "mex.h"
#include <cstring>

unsigned int np, ns;

void numerico_BC(double * res, double * ya, double * yb, double * P_bombeio, double *P_sinal)
{
    unsigned int index_sinal_f           = 2*np;
    unsigned int index_sinal_b           = 2*np + ns;
    unsigned int index_ase_f             = 2*np + 2*ns;
    unsigned int index_ase_b             = 2*np + 3*ns;  
    
    double * Ppump = new double[np];
    double * Sinal = new double[ns];
    
    unsigned int i = 0;
    for(i; i < np; i++)
        Ppump[i] = yb[i] - P_bombeio[i];
    
    for(i = 0; i < ns; i++)
        Sinal[i] = ya[index_sinal_f + i] - P_sinal[i];
        
    unsigned int S_np = np*sizeof(double);
    unsigned int S_ns = ns*sizeof(double);
    
    memcpy(res, Ppump, S_np);
    memcpy(res+np, ya + np, S_np);
    memcpy(res+index_sinal_f, Sinal, S_ns);
    memcpy(res+index_sinal_b,yb+index_sinal_b, S_ns);
    memcpy(res+index_ase_f, ya+index_ase_f, S_ns);
    memcpy(res+index_ase_b, yb+index_ase_b, S_ns);
    
    return;
}

//res = numerico_BC(ya,yb,bombeio,sinal)
void mexFunction(unsigned int nlhs,mxArray *plhs[],unsigned int nrhs,const mxArray *prhs[])
{
    mxArray *B_N, *B_P; //mxArray para bombeio
    mxArray *S_N, *S_P; //mxArray para sinal
    double *P_bombeio, *P_sinal;
    double * ya, * yb; //Entradas
    double * res; //Saída
      
    if(nrhs != 4)
        mexErrMsgTxt("A função deve ter 11 entradas");
    if(nlhs != 1)
        mexErrMsgTxt("A função deve ter uma única saída (dydx)");
    
    //Verifica parâmetros passados como entradas
    if(!mxIsStruct(prhs[2]))
        mexErrMsgTxt("A variável bombeio deve ser do tipo struct");
    else if(!mxIsStruct(prhs[3]))
        mexErrMsgTxt("A variável sinal deve ser do tipo struct");
    
    //Bombeio
    B_N = mxGetField(prhs[2], 0, "N_Bombeios");
    if (B_N == NULL)
        mexErrMsgTxt("Campo N_Bombeios de bombeio não encontrado.");
    B_P = mxGetField(prhs[2],0,"P");
    if (B_P == NULL)
        mexErrMsgTxt("Campo P de bombeio não encontrado.");
        
    np = (unsigned int)mxGetScalar(B_N);
    P_bombeio = mxGetPr(B_P);    
        
    //Sinal
    S_N = mxGetField(prhs[3],0,"N_Sinais");
    if (S_N == NULL)
        mexErrMsgTxt("Campo N_Sinais de sinal não encontrado.");
    S_P = mxGetField(prhs[3],0,"P");
    if (S_P == NULL)
        mexErrMsgTxt("Campo P de sinal não encontrado.");
    
    ns = (unsigned int)mxGetScalar(S_N);
    P_sinal = mxGetPr(S_P);
      
    ya = mxGetPr(prhs[0]);
    yb = mxGetPr(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(2*np + 4*ns, 1, mxREAL);
    res = mxGetPr(plhs[0]);
    
    numerico_BC(res,ya,yb,P_bombeio,P_sinal);
}