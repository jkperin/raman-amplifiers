    /*B_lambda = mxGetField(prhs[2],0,"lambda");
    if (B_lambda == NULL)
        mexErrMsgTxt("Campo lambda de bombeio não encontrado.");
    B_P = mxGetField(prhs[2],0,"P");
    if (B_P == NULL)
        mexErrMsgTxt("Campo P de bombeio não encontrado.");
    B_FPL = mxGetField(prhs[2],0,"FPL");
    if (B_FPL == NULL)
        mexErrMsgTxt("Campo FPL de bombeio não encontrado.");
    B_Bwp = mxGetField(prhs[2],0,"Bwp");
    if (B_FPL == NULL)
        mexErrMsgTxt("Campo Bwp de bombeio não encontrado.");*/
       
    bombeio.N_Bombeios = (int)mxGetScalar(B_N);
    bombeio.lambda = mxGetPr(B_lambda);
    bombeio.P = mxGetPr(B_P);
    bombeio.FPL = mxGetScalar(B_FPL);
    bombeio.Bwp = mxGetScalar(B_Bwp);
    
    
        S_lambda = mxGetField(prhs[3],0,"lambda");
    if (S_lambda == NULL)
        mexErrMsgTxt("Campo lambda de sinal não encontrado.");
    S_P = mxGetField(prhs[3],0,"P");
    if (S_P == NULL)
        mexErrMsgTxt("Campo P de sinal não encontrado.");
    S_Bws = mxGetField(prhs[3],0,"Bws");
    if (S_Bws == NULL)
        mexErrMsgTxt("Campo Bws de sinal não encontrado.");
        
    sinal.lambda = mxGetPr(S_lambda);
    sinal.P = mxGetPr(S_P);
    sinal.Bws = mxGetScalar(S_Bws);
    
        F_L = mxGetField(prhs[4], 0, "L");
    if (F_L == NULL)
        mexErrMsgTxt("Campo L de fibra não encontrado.");
    
    
    F_asdBkm = mxGetField(prhs[4], 0, "alfasdBkm");
    if (F_asdBkm == NULL)
        mexErrMsgTxt("Campo alfasdBkm de fibra não encontrado.");
    F_apdBkm = mxGetField(prhs[4], 0, "alfapdBkm");
    if (F_apdBkm == NULL)
        mexErrMsgTxt("Campo alfapdBkm de fibra não encontrado.");   
    F_Crs = mxGetField(prhs[4],0, "Crpicosinal");
    if (F_Crs == NULL)
        mexErrMsgTxt("Campo Crpicosinal de fibra não encontrado.");
    F_Crp = mxGetField(prhs[4], 0, "Crpicopump");
    if (F_Crp == NULL)
        mexErrMsgTxt("Campo Crpicopump de fibra não encontrado.");
    F_Crnorm = mxGetField(prhs[4], 0, "Crnormal");
    if (F_Crnorm == NULL)
        mexErrMsgTxt("Campo Crnormal de fibra não encontrado.");
    F_sf = mxGetField(prhs[4], 0, "sepfreq");
    if (F_sf == NULL)
        mexErrMsgTxt("Campo sepfreq de fibra não encontrado.");
    F_Aeffp = mxGetField(prhs[4], 0, "Aeffp");
    if (F_Aeffp == NULL)
        mexErrMsgTxt("Campo Aeffp de fibra não encontrado.");
    F_Aeffs = mxGetField(prhs[4], 0, "Aeffs");
    if (F_Aeffs == NULL)
        mexErrMsgTxt("Campo Aeffs de fibra não encontrado.");
    F_KR = mxGetField(prhs[4], 0, "KR");
    if (F_KR == NULL)
         mexErrMsgTxt("Campo KR de fibra não encontrado.");
    F_NA = mxGetField(prhs[4], 0, "NA");
    if (F_NA == NULL)
         mexErrMsgTxt("Campo NA de fibra não encontrado.");
    F_no = mxGetField(prhs[4], 0, "no");
    if (F_no == NULL)
         mexErrMsgTxt("Campo no de fibra não encontrado.");
    
    
        fibra.L = mxGetScalar(F_L);
    fibra.alfas = mxGetPr(F_alfas);
    fibra.alfap = mxGetPr(F_alfap);
    fibra.alfasdBkm = mxGetPr(F_asdBkm);
    fibra.alfapdBkm = mxGetPr(F_apdBkm);
    fibra.Crpicosinal = mxGetPr(F_Crs);
    fibra.Crpicopump = mxGetPr(F_Crp);
    fibra.Crnormal = mxGetPr(F_Crnorm);
    fibra.sepfreq = mxGetPr(F_sf);
    fibra.Aeffp = mxGetPr(F_Aeffp);
    fibra.Aeffs = mxGetPr(F_Aeffs);
    fibra.KR = mxGetScalar(F_KR);
    fibra.NA = mxGetScalar(F_NA);
    fibra.no = mxGetScalar(F_no);
    
    
    
    
    
    
    
    
    
    