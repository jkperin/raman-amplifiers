#include "mex.h"
#include <cstring>

typedef struct Fibra {
    double * alfas;
    double * alfap;
} Fibra;

typedef struct Parametros {
    double * Nes;
    double * Nep;
    double * Rss;
    double * Rpp;
    double * gsps;
    double * gpps;
    double * Taseps;
    double * gspp;
    double * gppp;
    double * Tasepp;
    double * gsss;
    double * gpss;
    double * Tasess;
} Parametros;

unsigned int np, ns;

inline int ConvIndexp(int i, int j)
{
    return i + j*np;
}

inline int ConvIndexs(int i, int j)
{
    return i + j*ns;
}

void numerico_ODE(double *dydx, double * y, const Fibra& fibra, const Parametros& Param)
{
    unsigned int index_sinal_f           = 2*np;
    unsigned int index_sinal_b           = 2*np + ns;
    unsigned int index_ase_f             = 2*np + 2*ns;
    unsigned int index_ase_b             = 2*np + 3*ns;
      
    //Aloca memória
    double * Psf = new double[ns];
    double * Psb = new double[ns];
    double * Ppf = new double[np];
    double * Ppb = new double[np];
    double * Pasef = new double[ns];
    double * Paseb = new double[ns];
    
    double Nf1 = 0.0, Nf2 = 0.0, Nf3 = 0.0, Nb1 = 0.0, Nb2 = 0.0, Nb3 = 0.0;
    double Nase1 = 0.0, Nase2 = 0.0, Naseb1 = 0.0, Naseb2 = 0.0;
    double Case1 = 0.0, Case2 = 0.0, Caseb1 = 0.0, Caseb2 = 0.0;
    double Np1f = 0.0, Np2f = 0.0, Np3f = 0.0;
    double Np1b = 0.0, Np2b = 0.0, Np3b = 0.0; 
    
    //Monta EDO
    
    for(int j = 0; j < ns; j++) // %j=1:ns
    {
        for(int k = 0; k < np; k++) //k=1:np
        {
            //Ganho do sinal devido aos pumps +z e -z
            Nf1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*y[index_sinal_f + j]; 
            Nb1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*y[index_sinal_b + j];
            //Ganho da ASE devido aos bombeios +z e -z   
            Nase1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*y[index_ase_f + j]; 
            Naseb1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*y[index_ase_b + j];
            //%Geração da ASE em torno dos bombeios +z e -z
            Case1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*2*Param.Nes[j]*Param.Taseps[ConvIndexp(k,j)]; 
            Caseb1 += Param.gsps[ConvIndexp(k,j)]*(y[k]+y[k+np])*2*Param.Nes[j]*Param.Taseps[ConvIndexp(k,j)];
        }

        for(int m = 0; m < ns; m++) //m=1:ns  
        {
            //Ganho do sinal devido aos outros sinais +z e -z
            Nf2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*y[index_sinal_f + j]; 
            Nb2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*y[index_sinal_b + j];
            //mexPrintf("%.20f;%.20f\n", Nf2,Nb2);
            // Deplecao do sinal devido aos outros sinais
            Nf3 += Param.gpss[ConvIndexs(j,m)]*((y[m+index_sinal_f]+y[m+index_sinal_b]+ 4*Param.Nes[m]*Param.Tasess[ConvIndexs(j,m)])*y[index_sinal_f + j]); 
            Nb3 += Param.gpss[ConvIndexs(j,m)]*((y[m+index_sinal_f]+y[m+index_sinal_b]+ 4*Param.Nes[m]*Param.Tasess[ConvIndexs(j,m)])*y[index_sinal_b + j]);
            // Ganho da ASE devido aos bombeios +z e -z   
            Nase2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*y[index_ase_f + j];    
            Naseb2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*y[index_ase_b + j];
            // Geração da ASE em torno dos bombeios +z e -z
            Case2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*2*Param.Nes[j]*Param.Tasess[ConvIndexs(j,m)];
            Caseb2 += Param.gsss[ConvIndexs(j,m)]*(y[m+index_sinal_f]+y[m+index_sinal_b])*2*Param.Nes[j]*Param.Tasess[ConvIndexs(j,m)];
        }

        // Eq. 2.31 (Tese da Shirley) - Para sinais
        // -Atenuação na fibra + Espalhamento duplo de Rayleigh + Ganho devido a 
        // bombeios + Ganho devido a sinais - Depleção devido aos sinais
        Psf[j] = -fibra.alfas[j]*y[index_sinal_f + j] + Param.Rss[j]*y[index_sinal_b + j] + Nf1 + Nf2 - Nf3; // +z
        // Atenuação na fibra - Espalhamento duplo de Rayleigh - Ganho devido a 
        // bombeios - Ganho devido a sinais + Depleção devido aos sinais 
        Psb[j] = fibra.alfas[j]*y[index_sinal_b + j]- Param.Rss[j]*y[index_sinal_f + j] - Nb1 - Nb2 + Nb3; // -z

        // Eq. 2.31 (Tese da Shirley) - Para ASE 
        // -Atenuaçao da fibra + Espalhamento duplo de Rayleigh + Ganho da ASE 
        // devido aos bombeios + Ganho da ASE devido aos sinais + Amplificação
        // da ASE em torno dos bombeios + Amplificação da ASE em torno do Sinais
        Pasef[j]=-fibra.alfas[j]*y[index_ase_f + j] + Param.Rss[j]*y[index_ase_b + j] + Nase1 + Nase2 + Case1 + Case2; //+z
        // +Atenuaçao da fibra - Espalhamento duplo de Rayleigh - Ganho da ASE 
        // devido aos bombeios - Ganho da ASE devido aos sinais - Amplificação
        // da ASE em torno dos bombeios - Amplificação da ASE em torno dos Sinais
        Paseb[j]= +fibra.alfas[j]*y[index_ase_b + j] - Param.Rss[j]*y[index_ase_f + j] - Naseb1 - Naseb2 - Caseb1 - Caseb2; //-z

        Nf1 = Nf2 = Nf3 = Nb1 = Nb2 = Nb3 = 0.0;
        Nase1 = Nase2 = Naseb1 = Naseb2 = 0.0;
        Case1 = Case2 = Caseb1 = Caseb2 = 0.0; 
    }
  
    for(int k = 0; k < np; k++)
    {
        for(int j = 0; j < ns; j++)
        {
            //Depleção dos bombeios devido aos sinais(+z e -z), ASE(+z e -z) e
            //Depleção devido a amplificação da ASE nas menores frequências(+z e -z)
            Np1f += Param.gpps[ConvIndexp(k,j)]*(y[index_sinal_f + j]+y[index_sinal_b + j]+y[index_ase_f + j]+y[index_ase_b + j]+4*Param.Nes[j]*Param.Taseps[ConvIndexp(k,j)])*y[k]; 
            Np1b += Param.gpps[ConvIndexp(k,j)]*(y[index_sinal_f + j]+y[index_sinal_b + j]+y[index_ase_f + j]+y[index_ase_b + j]+4*Param.Nes[j]*Param.Taseps[ConvIndexp(k,j)])*y[k+np];
        }

        for(int n = 0; n < np; n++)
        {
            //Deplecao dos bombeios devido aos outro bombeios(+z e -z) e 
            //Depleção devido a amplificação da ASE nas menores frequências(+z e -z)
            Np2f += Param.gppp[ConvIndexp(k,n)]*(y[n]+y[n+np]+ 4*Param.Nep[n]*Param.Tasepp[ConvIndexp(k,n)])*y[k]; 
            Np2b += Param.gppp[ConvIndexp(k,n)]*(y[n]+y[n+np]+ 4*Param.Nep[n]*Param.Tasepp[ConvIndexp(k,n)])*y[k+np];

            //Ganho dos bombeios devidos aos outros bombeios (+z  e -z)
            Np3f += Param.gspp[ConvIndexp(k,n)]*(y[n]+y[n+np])*y[k]; 
            Np3b += Param.gspp[ConvIndexp(k,n)]*(y[n]+y[n + np])*y[k+np]; 
        }
        //Eq. 2.31 (Tese da Shirley) - Para Bombeios 
        // Os bombeios são contra-propagantes
        // Atenuaçao da fibra - Espalhamento duplo de Rayleigh + Depleção dos
        // bombeios devidos aos outro sinais + Depleção devido aos outros
        // bombeios - Ganho devido aos outros bombeios
        Ppf[k] = fibra.alfap[k]*y[k] - Param.Rpp[k]*y[k+np] + Np1f + Np2f - Np3f; //-z
        // -Atenuaçao da fibra + Espalhamento duplo de Rayleigh - Depleção dos
        // bombeios devidos aos outro sinais - Depleção devido aos outros
        // bombeios + Ganho devido aos outros bombeios
        Ppb[k] = -fibra.alfap[k]*y[k+np] + Param.Rpp[k]*y[k] - Np1b - Np2b + Np3b; //+z
                
        Np1f = Np2f = Np3f = Np1b = Np2b = Np3b = 0.0;
    }
    
    unsigned int S_np = np*sizeof(double);
    unsigned int S_ns = ns*sizeof(double);
       
    memcpy(dydx,Ppf,S_np);
    memcpy(dydx+np,Ppb,S_np);
    memcpy(dydx+index_sinal_f,Psf,S_ns);
    memcpy(dydx+index_sinal_b,Psb,S_ns);
    memcpy(dydx+index_ase_f,Pasef,S_ns);
    memcpy(dydx+index_ase_b,Paseb,S_ns);
    
    //Limpa Memória
    delete [] Psf;
    delete [] Psb;
    delete [] Ppf;
    delete [] Ppb;
    delete [] Pasef;
    delete [] Paseb;
}

/*************************************************************************/
/************************* GATEWAY ROUTINE *******************************/
/*************************************************************************/

//function dydx = numerico_ODE(x,y,bombeio,sinal,fibra,Param)

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mxArray *B_N; //mxArray para bombeio
    mxArray *S_N; //mxArray para sinal
    mxArray *F_alfas, *F_alfap; //mxArray para fibra
    mxArray *P_Nes, *P_Nep, *P_Rss, *P_Rpp, *P_gsps, *P_gpps, *P_Taseps, *P_gspp, *P_gppp, *P_Tasepp, *P_gsss, *P_gpss, *P_Tasess; //mxArray para Param
    double * y; //Entradas
    double * dydx; //Saída
    
    Fibra fibra;
    Parametros Param;
   
    if(nrhs != 6)
        mexErrMsgTxt("A função deve ter 11 entradas");
    if(nlhs != 1)
        mexErrMsgTxt("A função deve ter uma única saída (dydx)");
    
    //Verifica parâmetros passados como entradas
    if(!mxIsStruct(prhs[2]))
        mexErrMsgTxt("A variável bombeio deve ser do tipo struct");
    else if(!mxIsStruct(prhs[3]))
        mexErrMsgTxt("A variável sinal deve ser do tipo struct");
    else if(!mxIsStruct(prhs[4]))
        mexErrMsgTxt("A variável fibra deve ser do tipo struct");
    else if(!mxIsStruct(prhs[5]))
        mexErrMsgTxt("A variável Param deve ser do tipo struct");
    
    //Bombeio
    B_N = mxGetField(prhs[2], 0, "N_Bombeios");
    if (B_N == NULL)
        mexErrMsgTxt("Campo N_Bombeios de bombeio não encontrado.");
    
    np = (unsigned int)mxGetScalar(B_N);
        
    //Sinal
    S_N = mxGetField(prhs[3],0,"N_Sinais");
    if (S_N == NULL)
        mexErrMsgTxt("Campo N_Sinais de sinal não encontrado.");

    ns = (unsigned int)mxGetScalar(S_N);

    //Fibra
    F_alfas = mxGetField(prhs[4], 0, "alfas");
    if (F_alfas == NULL)
        mexErrMsgTxt("Campo alfas de fibra não encontrado.");
    F_alfap = mxGetField(prhs[4], 0, "alfap");
    if (F_alfap == NULL)
        mexErrMsgTxt("Campo alfap de fibra não encontrado.");
    
    fibra.alfas = mxGetPr(F_alfas);
    fibra.alfap = mxGetPr(F_alfap);
    
    P_Nes = mxGetField(prhs[5], 0, "Nes");
    if (P_Nes == NULL)
         mexErrMsgTxt("Campo Nes de fibra não encontrado.");
    P_Nep = mxGetField(prhs[5], 0, "Nep");
    if (P_Nep == NULL)
         mexErrMsgTxt("Campo Nep de Param não encontrado.");
    P_Rss = mxGetField(prhs[5], 0, "Rss");
    if (P_Rss == NULL)
         mexErrMsgTxt("Campo Rss de Param não encontrado.");
    P_Rpp = mxGetField(prhs[5], 0, "Rpp");
    if (P_Rpp == NULL)
         mexErrMsgTxt("Campo Rpp de Param não encontrado.");
    P_gsps = mxGetField(prhs[5], 0, "gsps");
    if (P_gsps == NULL)
         mexErrMsgTxt("Campo gsps de Param não encontrado.");
    P_gpps = mxGetField(prhs[5], 0, "gpps");
    if (P_gpps == NULL)
         mexErrMsgTxt("Campo gpps de Param não encontrado.");
    P_Taseps = mxGetField(prhs[5], 0, "Taseps");
    if (P_Taseps == NULL)
         mexErrMsgTxt("Campo Taseps de Param não encontrado.");
    P_gspp = mxGetField(prhs[5], 0, "gspp");
    if (P_gspp == NULL)
         mexErrMsgTxt("Campo gspp de Param não encontrado.");
    P_gppp = mxGetField(prhs[5], 0, "gppp");
    if (P_gppp == NULL)
         mexErrMsgTxt("Campo gppp de Param não encontrado.");
    P_Tasepp = mxGetField(prhs[5], 0, "Tasepp");
    if (P_Tasepp == NULL)
         mexErrMsgTxt("Campo Tasepp de Param não encontrado.");
    P_gsss = mxGetField(prhs[5], 0, "gsss");
    if (P_gsss == NULL)
         mexErrMsgTxt("Campo gsss de Param não encontrado.");
    P_gpss = mxGetField(prhs[5], 0, "gpss");
    if (P_gpss == NULL)
         mexErrMsgTxt("Campo gpss de Param não encontrado.");
    P_Tasess = mxGetField(prhs[5], 0, "Tasess");
    if (P_Tasess == NULL)
         mexErrMsgTxt("Campo Tasess de Param não encontrado.");
    
    Param.Nes = mxGetPr(P_Nes);
    Param.Nep = mxGetPr(P_Nep);
    Param.Rss = mxGetPr(P_Rss);
    Param.Rpp = mxGetPr(P_Rpp);
    Param.gsps = mxGetPr(P_gsps);
    Param.gpps = mxGetPr(P_gpps);
    Param.Taseps = mxGetPr(P_Taseps);
    Param.gspp = mxGetPr(P_gspp);
    Param.gppp = mxGetPr(P_gppp);
    Param.Tasepp = mxGetPr(P_Tasepp);
    Param.gsss = mxGetPr(P_gsss);
    Param.gpss = mxGetPr(P_gpss);
    Param.Tasess = mxGetPr(P_Tasess);
        
    y = mxGetPr(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(2*np + 4*ns, 1, mxREAL);
    dydx = mxGetPr(plhs[0]);
    
    numerico_ODE(dydx,y,fibra,Param);
}