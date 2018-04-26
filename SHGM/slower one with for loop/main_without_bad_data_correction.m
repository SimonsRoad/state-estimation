clear all
%LOADING THE DATA
%load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;
%load Dataset_B.mat;
load Measurement.mat

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=20; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
c1=2.7;

%initial start point
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c,1);
iter=0;

a=size(meas.BusV,2);%getting the total number of measurements
meas_new=meas;
%SOLTUION ALGORITHM STARTS HERE
for i=1:a
    meas.BusInjP=meas_new.BusInjP(:,i);
    meas.BusInjQ=meas_new.BusInjQ(:,i);
    meas.BusV=meas_new.BusV(:,i);
    meas.LineP_FT=meas_new.LineP_FT(:,i);
    meas.LineP_TF=meas_new.LineP_TF(:,i);
    meas.LineQ_FT=meas_new.LineQ_FT(:,i);
    meas.LineQ_TF=meas_new.LineQ_TF(:,i);
    meas.std_BusInjP=meas_new.std_BusInjP(:,i);
    meas.std_BusInjQ=meas_new.std_BusInjQ(:,i);
    meas.std_BusV=meas_new.std_BusV(:,i);
    meas.std_LineP_FT=meas_new.std_LineP_FT(:,i);
    meas.std_LineP_TF=meas_new.std_LineP_TF(:,i);
    meas.std_LineQ_FT=meas_new.std_LineQ_FT(:,i);
    meas.std_LineQ_TF=meas_new.std_LineQ_TF(:,i);
    %getting Jacobian sparsity pattern (done only once)
    [topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2017(topo, meas);

    %obtaining the admittance matrix in the compressed row form (done only once)
    Y_bus = f_Y_bus_compressed_v2017(topo);

    %constructing vector z, which is a vector of available measurements
    %NOTE: different types of measurements should be in the same order as in the handout
    %STUDENT CODE 1
    z=[meas.BusV(ind_meas.ind_meas_V); meas.BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.BusInjQ(ind_meas.ind_meas_Qinj); meas.LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.LineP_TF(ind_meas.ind_meas_Pji); meas.LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.LineQ_TF(ind_meas.ind_meas_Qji)];
    N_total=length(z);
    H = f_measJac_H_v2017( Ve, thetae, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);%for calculating the non-zero values in each row of H
    for i=1:N_total
        m=full((sum(H(i,:)~=0)));%non-zero values of row i of H
        b(i)=chi2inv(0.975,m)^(0.5);
    end

    %NOTE: matrices R must be constructed in sparse format!!!
    %STUDENT CODE 2
    R_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
    R_s=sparse(1:N_total,1:N_total,(((R_diag).^(2)).^(-0.5)),N_total,N_total);%R_s=R^(-1/2)
    W=sparse(1:N_total,1:N_total,1./((R_diag).^2),N_total,N_total);
    Wsqrt=sparse(1:N_total,1:N_total,1./(R_diag),N_total,N_total);

    [ Ve, thetae, eps_all, time, convergence ]=f_SE_NR_algorithm_v2017 ( Ve, thetae, topo, Y_bus, z, W, R_s,... 
        R_diag, b, c1, ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, N_total);
    state(:,i)=[thetae;Ve];
end