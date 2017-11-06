#include "mex.h"
#include <math.h>

#define Vluz 2.99792458e8
#define res_L 100 //número de pontos do vetor z
#define pi 3.141592653589793

typedef struct Bombeio{
    int N_Bombeios;
    double * lambda;
    double * P;
    double FPL;
} Bombeio;

typedef struct Sinal {
    double * lambda;
    int N_Sinais;
} Sinal;

typedef struct Fibra {
    double L;
    double * alfas;
    double * alfap;
    double * Crpicopump;
    double * Crnormal;
    double * sepfreq;
} Fibra;

/*************************** vet.h - inicio ******************************/
void mul_vet(double * result, double * vet, double num, int size)
{
	for(int i = 0;i < size;i++)
		*(result+i) = *(vet+i)*num;
    return;
}

void div_vet(double * result, double * vet, double num, int size)
{
	for(int i=0; i < size;i++)
		*(result+i) = num/(*(vet+i));
    return;
}

void zera(double * vet, int size)
{
	for(int i=0;i<size;i++)
		*(vet+i) = 0.0;
    return;
}

/**************************** vet.h - fim ********************************/


/****************************  matlab.h - inicio *************************/

int length(const double * vet)
{
	int i = 0;
	while(vet[i] >= 0.0) i++;
	return ++i;
}

double trapz(const double * eml_x, const double * eml_y, int size)
{
	double * eml_b_x = new double [size];
	double result, eml_ylast;
	int eml_i, eml_iy, eml_k;

	for(eml_i = 0; eml_i < size; eml_i++) 
		eml_b_x[eml_i] = eml_x[eml_i];

	for(eml_i = 0; eml_i < size -1; eml_i++)
		eml_b_x[eml_i] = eml_b_x[eml_i + 1] - eml_b_x[eml_i];

	result = 0.0;
	eml_i = 0;
	eml_iy = 1;
	eml_ylast = eml_y[0];
	for(eml_k = 2; eml_k < 101; eml_k++) 
	{
		eml_iy++;
		eml_i++;
		result += eml_b_x[eml_i - 1] * ((eml_ylast + eml_y[eml_iy - 1]) / 2.0);
		eml_ylast = eml_y[eml_iy - 1];
	}
    delete [] eml_b_x;
	return result;
}

double min(const double * vet, int size)
{
	double tmp = vet[0];
	for(int i = 1; i < size;i++)
	{
		if(tmp > vet[i])
			tmp = vet[i];
		else
			continue;
	}
	return tmp;
}

double max(const double * vet, int size)
{
	double tmp = vet[0];
	for(int i = 1;i < size;i++)
	{
		if(tmp < vet[i])
			tmp = vet[i];
		else
			continue;
	}
	return tmp;
}

double mean(const double * vet, int size)
{
	double soma = 0;
	for(int i = 0; i < size;i++)
		soma += *(vet+i);
	return  soma/size;
}

void linspace(double z[], double ini, double fim)
{
	double deltaz = (fim - ini)/(res_L-1);
	z[0] = 0.0;
	for(int i = 1;i < res_L;i++)
		z[i] = deltaz + z[i-1];
    return;
} 

/****************************  matlab.h - fim ****************************/


/*************************************************************************/
/*************************** DRA_ANALITICO *******************************/
/*************************************************************************/

void DRA_Analitico(const Bombeio& bombeio, const Sinal& sinal, const Fibra& fibra, double& ripple, double& Ganho_Medio, double& Ganho_on_off_medio, double * Ganho_on_off, double * Ppump0)
{
    double alfap = mean(fibra.alfap,bombeio.N_Bombeios);  //(neper/m)
        
	int N_Cr = length(fibra.Crnormal);
	double * fp, *fs, *wp;
	fp = new double [bombeio.N_Bombeios];
	wp = new double [bombeio.N_Bombeios];
	fs = new double [sinal.N_Sinais];
	
	div_vet(fp,bombeio.lambda,Vluz*1e9,bombeio.N_Bombeios);
	div_vet(fs,sinal.lambda,Vluz*1e9,sinal.N_Sinais);
	mul_vet(wp,fp,2*pi,bombeio.N_Bombeios);

	int i,j,k,h;

/*************************************************************************/
/****************** CRNOVO ENTRE PUMPS E SINAIS **************************/
/*************************************************************************/

	double ** deltafps = new double* [bombeio.N_Bombeios];
	double ** Crnovops = new double* [bombeio.N_Bombeios];
	double ** Cr_upps = new double* [bombeio.N_Bombeios];
	for(i = 0;i < bombeio.N_Bombeios; i++)
	{
		deltafps[i] = new double[sinal.N_Sinais];
		Crnovops[i] = new double[sinal.N_Sinais];
		Cr_upps[i] = new double[sinal.N_Sinais];
	}

	for(i=0;i<sinal.N_Sinais;i++)
		for(j=0;j<bombeio.N_Bombeios;j++)
		{
			deltafps[j][i] = abs(fp[j]-fs[i]);
			for(k=0;k<N_Cr-1;k++)
				if((fibra.sepfreq[k] <= abs(deltafps[j][i]))&&(fibra.sepfreq[k+1] >= abs(deltafps[j][i])))
				{
					Crnovops[j][i] = fibra.Crnormal[k] + (abs(deltafps[j][i]) - fibra.sepfreq[k])*(fibra.Crnormal[k+1] - fibra.Crnormal[k])/(fibra.sepfreq[k+1] - fibra.sepfreq[k]);
					break;
				}
			Cr_upps[j][i] = fibra.Crpicopump[j]*Crnovops[j][i];
		}


/*************************************************************************/
/********************** CRNOVO ENTRE PUMPS *******************************/
/*************************************************************************/

	double ** deltafpp = new double* [bombeio.N_Bombeios];
	double ** Crnovopp = new double* [bombeio.N_Bombeios];
	double ** Cr_uppp = new double* [bombeio.N_Bombeios];
	for(i = 0;i < bombeio.N_Bombeios; i++)
	{
		deltafpp[i] = new double[bombeio.N_Bombeios];
		Crnovopp[i] = new double[bombeio.N_Bombeios];
		Cr_uppp[i] = new double[bombeio.N_Bombeios];
	}

	for(i=0;i<bombeio.N_Bombeios;i++)
		for(j=0;j<bombeio.N_Bombeios;j++)
		{
			deltafpp[j][i] = abs(fp[j]-fp[i]);
			for(k=0;k<N_Cr-1;k++)
				if((fibra.sepfreq[k] <= abs(deltafpp[j][i]))&&(fibra.sepfreq[k+1] >= abs(deltafpp[j][i])))
				{
					Crnovopp[j][i] = fibra.Crnormal[k] + (abs(deltafpp[j][i]) - fibra.sepfreq[k])*(fibra.Crnormal[k+1] - fibra.Crnormal[k])/(fibra.sepfreq[k+1] - fibra.sepfreq[k]);
					break;
				}
			Cr_uppp[j][i] = fibra.Crpicopump[j]*Crnovopp[j][i];
		}

/*************************************************************************/
/************** EQUACAO DOS PUMPS COM ITERACAO ENTRE ELES ****************/
/*************************************************************************/
       
	double z[res_L];
    linspace(z, 0.0, fibra.L);

/**************** PROPAGAÇAO DOS BOMBEIOS ****************************/

    double Pp_g = 0.0, Pp_l = 0.0;
	double Np_g = 0.0, Np_l = 0.0;
	double Pp_p[res_L];
	double Pg_p[res_L];
	double ** Pp = new double* [bombeio.N_Bombeios];
	for(j =0;j<bombeio.N_Bombeios;j++)
		Pp[j] = new double[res_L];

	for(j=0;j<bombeio.N_Bombeios;j++)
	{
		zera(Pp_p,res_L);
		zera(Pg_p,res_L);
        for(i=0;i<bombeio.N_Bombeios;i++)
		{
            if(i>j)
			{
				for(k=0;k<bombeio.N_Bombeios;k++)
				{
                    if(k>i)
                        Pp_l += (wp[i]/wp[k])*Cr_uppp[i][k]*bombeio.P[k];
					else if(k<i)
                        Pp_g += Cr_uppp[i][k]*bombeio.P[k];
				}

				for(h=0;h<res_L;h++)
					Pp_p[h] -= ((wp[j]/wp[i])*Cr_uppp[j][i]*bombeio.P[i]*(1-exp(-1/(alfap*bombeio.FPL)*(1-exp(-alfap*(fibra.L-z[h])))*(-Pp_g+Pp_l))))/(-Pp_g+Pp_l);
				
				Pp_g = 0.0, Pp_l = 0.0;
			}
            else if (i<j)
			{
				for(k=0;k<bombeio.N_Bombeios;k++)
				{
					if (k>i)
                        Np_l += (wp[i]/wp[k])*Cr_uppp[i][k]*bombeio.P[k];
                    else if (k<i)
                        Np_g += Cr_uppp[i][k]*bombeio.P[k];
				}

				for(h=0;h<res_L;h++)
					Pg_p[h] += (Cr_uppp[j][i]*bombeio.P[i]*(1-exp(-1/(alfap*bombeio.FPL)*(1-exp(-alfap*(fibra.L-z[h])))*(-Np_g+Np_l))))/(-Np_g+Np_l);
                
				Np_g=0.0, Np_l=0.0;
			}
		}
		for(h = 0;h<res_L;h++)
			Pp[j][h] = bombeio.P[j]*exp(-alfap*(fibra.L-z[h]))*exp(Pp_p[h] + Pg_p[h]);
	}
    
    //Ppump0 = transp(Pp(:,1));
    for(i = 0; i < bombeio.N_Bombeios; i++)
        Ppump0[i] = Pp[i][0];

/*************************************************************************/
/****** GANHO ANALITICO OBTIDO COM A SEGUNDA ITERAÇAO DOS BOMBEIOS *******/
/*************************************************************************/
	
	double * Ganho_sem_RamandB =  new double [sinal.N_Sinais];
    double * GA_sinaldB = new double [sinal.N_Sinais];
	double ** K = new double* [bombeio.N_Bombeios];
	for(h=0;h<bombeio.N_Bombeios;h++)
		K[h] = new double [sinal.N_Sinais];
	double ganho_pump;

	for(i=0;i<sinal.N_Sinais;i++)
	{
		ganho_pump = 0.0;
		for(j = 0;j<bombeio.N_Bombeios;j++)
		{
			K[j][i] = Cr_upps[j][i]/bombeio.FPL;
			ganho_pump += K[j][i]*trapz(z,Pp[j],res_L);
		}
		Ganho_sem_RamandB[i] = exp(-fibra.alfas[i]*fibra.L); //neper
		Ganho_sem_RamandB[i] = 10*log10(Ganho_sem_RamandB[i]);

         GA_sinaldB[i]= exp(-fibra.alfas[i]*fibra.L+ganho_pump); //neper
         GA_sinaldB[i] = 10*log10(GA_sinaldB[i]);

         Ganho_on_off[i] = GA_sinaldB[i]-Ganho_sem_RamandB[i];
	}

/*************************************************************************/
/***********************         SAÍDA           *************************/
/*************************************************************************/

	ripple                  = max(GA_sinaldB,sinal.N_Sinais)-min(GA_sinaldB,sinal.N_Sinais);
    Ganho_Medio             = mean(GA_sinaldB,sinal.N_Sinais);
    Ganho_on_off_medio      = mean(Ganho_on_off,sinal.N_Sinais);


/************************ LIMPA MEMÓRIA **********************************/

	delete [] Ganho_sem_RamandB;
    delete [] GA_sinaldB;
	delete [] fp;
	delete [] fs;
	delete [] wp;

	for(i = 0;i < bombeio.N_Bombeios; i++)
	{
		delete [] deltafps[i];
		delete [] Crnovops[i];
		delete [] Cr_upps[i];
		delete [] deltafpp[i];
		delete [] Crnovopp[i];
		delete [] Cr_uppp[i];
		delete [] K[i];
        delete [] Pp[i];
	}
    
    delete [] deltafps;
    delete [] Crnovops;
    delete [] Cr_upps;
    delete [] deltafpp;
    delete [] Crnovopp;
    delete [] Cr_uppp;
    delete [] K;
    delete [] Pp;
    
    return;
}

/*************************************************************************/
/*********************** GATEWAY ROUTINE *********************************/
/*************************************************************************/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mxArray *B_N, *B_lambda, *B_P, *B_FPL; //mxArray para bombeio
    mxArray *S_N, *S_lambda; //mxArray para sinal
    mxArray *F_alfas, *F_alfap, *F_Crp, *F_Crnorm, *F_sf, *F_L; //mxArray para fibra
    double ripple, Ganho_Medio, Ganho_on_off_medio, * Ganho_on_off, * Ppump0; //Saída
    
    Bombeio bombeio;
    Sinal sinal;
    Fibra fibra;
   
    if(nrhs != 3)
        mexErrMsgTxt("A função deve ter três entradas (bombeio,sinal,fibra)");
    if(nlhs != 5)
        mexErrMsgTxt("A função deve ter quatro saídas (ripple,Ganho_Medio,Ganho_on_off_medio, GA_sinaldB, Ppump0)");

    if((!mxIsStruct(prhs[0])) || (!mxIsStruct(prhs[1])) || (!mxIsStruct(prhs[2])))
        mexErrMsgTxt("As entradas devem ser do tipo struct");
        
    
    //Bombeio
    //Não utilizados: bombeio.Bwp
    B_N = mxGetField(prhs[0], 0, "N_Bombeios");
    if (B_N == NULL)
        mexErrMsgTxt("Campo N_Bombeios de bombeio não encontrado.");
    B_lambda = mxGetField(prhs[0],0,"lambda");
    if (B_lambda == NULL)
        mexErrMsgTxt("Campo lambda de bombeio não encontrado.");
    B_P = mxGetField(prhs[0],0,"P");
    if (B_P == NULL)
        mexErrMsgTxt("Campo P de bombeio não encontrado.");
    B_FPL = mxGetField(prhs[0],0,"FPL");
    if (B_FPL == NULL)
        mexErrMsgTxt("Campo FPL de bombeio não encontrado.");
       
    bombeio.N_Bombeios = (int)mxGetScalar(B_N);
    bombeio.lambda = mxGetPr(B_lambda);
    bombeio.P = mxGetPr(B_P);
    bombeio.FPL = mxGetScalar(B_FPL);
    
    //Sinal
    //Não utilizados: sinal.Bws
    S_N = mxGetField(prhs[1],0,"N_Sinais");
    if (S_N == NULL)
        mexErrMsgTxt("Campo N_Sinais de sinal não encontrado.");
    S_lambda = mxGetField(prhs[1],0,"lambda");
    if (S_lambda == NULL)
        mexErrMsgTxt("Campo lambda de sinal não encontrado.");
    
    sinal.N_Sinais = (int)mxGetScalar(S_N);
    sinal.lambda = mxGetPr(S_lambda);
    
	//Fibra
    //Não utilizados: fibra.lamb, fibra.alfasdBkm, fibra.alfapdBkm, fibra.Crpicosinal, fibra.Aeff
    F_L = mxGetField(prhs[2], 0, "L");
    if (F_L == NULL)
        mexErrMsgTxt("Campo L de fibra não encontrado.");
    F_alfas = mxGetField(prhs[2], 0, "alfas");
    if (F_alfas == NULL)
        mexErrMsgTxt("Campo alfas de fibra não encontrado.");
    F_alfap = mxGetField(prhs[2], 0, "alfap");
    if (F_alfap == NULL)
        mexErrMsgTxt("Campo alfap de fibra não encontrado.");
    F_Crp = mxGetField(prhs[2], 0, "Crpicopump");
    if (F_Crp == NULL)
        mexErrMsgTxt("Campo Crpicopump de fibra não encontrado.");
    F_Crnorm = mxGetField(prhs[2], 0, "Crnormal");
    if (F_Crnorm == NULL)
        mexErrMsgTxt("Campo Crnormal de fibra não encontrado.");
    F_sf = mxGetField(prhs[2], 0, "sepfreq");
    if (F_sf == NULL)
        mexErrMsgTxt("Campo sepfreq de fibra não encontrado.");

    fibra.L = mxGetScalar(F_L);
    fibra.alfas = mxGetPr(F_alfas);
    fibra.alfap = mxGetPr(F_alfap);
    fibra.Crpicopump = mxGetPr(F_Crp);
    fibra.Crnormal = mxGetPr(F_Crnorm);
    fibra.sepfreq = mxGetPr(F_sf);
    
    plhs[3] = mxCreateDoubleMatrix(1, sinal.N_Sinais, mxREAL);
    
    plhs[4] = mxCreateDoubleMatrix(1, bombeio.N_Bombeios, mxREAL);
    
    Ganho_on_off = mxGetPr(plhs[3]);
    Ppump0 = mxGetPr(plhs[4]);
    
    DRA_Analitico(bombeio,sinal,fibra,ripple,Ganho_Medio,Ganho_on_off_medio,Ganho_on_off,Ppump0);

    plhs[0] = mxCreateDoubleScalar(ripple);
    plhs[1] = mxCreateDoubleScalar(Ganho_Medio);
    plhs[2] = mxCreateDoubleScalar(Ganho_on_off_medio);
}
