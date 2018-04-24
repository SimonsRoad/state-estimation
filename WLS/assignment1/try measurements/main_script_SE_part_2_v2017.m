clc;clear;
%LOADING DATASET
load Dataset_B.mat;


%HERE COME USER-DEFINED PARAMETERS OF GN algorithm (using default options is recommended)
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=2;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid
Errortype=1;%make sure go into the loop
y=1;

%% DATA PREPARATION
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

%constructing the matrix of weights W using standard deviations provided in structure 'meas'
%NOTE: matrix W must be constructed in sparse format!!!
%STUDENT CODE 2
W_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
    meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
    meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
    meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
W=sparse(1:N_total,1:N_total,1./((W_diag).^2),N_total,N_total);
Wsqrt=sparse(1:N_total,1:N_total,1./(W_diag),N_total,N_total);


%% INITIAL SE SOLUTION
fprintf('INITIAL STATE ESTIMATION: \n')
%choosing the initial point
%STUDENT CODE 3
V=ones(topo.nBus,1);
theta=zeros(topo.nBus,1);
%calling the GN algorithm and recording initial solution of SE
[ V_SE0, theta_SE0, ~, ~, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
V_SE1=V_SE0;
theta_SE1=theta_SE0;
%loop for bad data detection
while Errortype==1
%% BAD DATA DETECTION


fprintf('SEARCHING FOR BAD DATA: \n')
%getting most recent H and h
[ H ] = f_measJac_H_v2017( V_SE1, theta_SE1, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
[ h ] = f_measFunc_h_v2017( V_SE1, theta_SE1, Y_bus, topo, ind_meas, N_meas);

%computing the required matrices
G = H'*W*H;
K = H*(G\(H'*W));
S = (eye(N_total) - K);
%computing vector of residuals
r=z-h;

% identifying critical measurements
i_critical = abs(diag(S)) < 1e-12;

%getting the diagonal elements of matrices S and W
S_diag = diag(S);
W_diag = diag(W);

% creating vector for normalized residuals
r_n = zeros(N_total,1);

% computing normalized residuals for redundant measurements
r_n(~i_critical) = abs(r(~i_critical)) ./ sqrt(S_diag(~i_critical)./W_diag(~i_critical));

%getting the value and index of largest normalized residual
[r_max,r_ind] = max(r_n);
fprintf('Measurement with index %d has the largest residual, the value of which is %.2f.\n',r_ind,r_max)

%checking the value of the largest normalized residual
if r_max > 4 %value 4 is taken from experience
    %residual too high - error in redundant measurement
    Errortype = 1;
    fprintf('There is a correctable error in the data set.\n\n\n')
else
    if meas.warning==1 %if there is a warning from SCADA system
        %critical bad measurement
        Errortype = 2;
        fprintf('A critical measurement is incorrect.\n\n\n')
    else %if there is no warning from SCADA system
        %good measurement
        Errortype = 0;
        fprintf('The data set is good.\n\n\n')
    end
end

%plotting residuals (and matrix S if need be)
%surf(S);
figure(y);
bar(1:N_total,r_n);
xlabel('Measurements');
ylabel('Residuals');


%% BAD DATA CORRECTION
fprintf('CORRECTING BAD DATA: \n')
%if there is bad data in the redundant measurement, it is possible to delete it.
if (Errortype ~= 1)
    if (Errortype == 2)
        fprintf('The error is in a critical measurement -> impossible to correct.\n\n')
    else
        fprintf('There is no bad data in the measurement set.\n\n')
    end
    return;
end

%deleting the faulty measurement from the dataset
if Errortype == 1
    if (r_ind<=N_meas.N_meas_V)
        meas.BusV(ind_meas.ind_meas_V(r_ind))=NaN;
    elseif (r_ind<=N_meas.N_meas_V+N_meas.N_meas_Pinj)
        meas.BusInjP(ind_meas.ind_meas_Pinj(r_ind-N_meas.N_meas_V))=NaN;
    elseif (r_ind<=N_meas.N_meas_V+N_meas.N_meas_Pinj+N_meas.N_meas_Qinj)
        meas.BusInjP(ind_meas.ind_meas_Qinj(r_ind-N_meas.N_meas_Pinj-N_meas.N_meas_V))=NaN;
    elseif (r_ind<=N_meas.N_meas_V+N_meas.N_meas_Pinj+N_meas.N_meas_Qinj+N_meas.N_meas_Pij)
        meas.LineP_FT(ind_meas.ind_meas_Pij(r_ind-N_meas.N_meas_Pinj-N_meas.N_meas_Qinj-N_meas.N_meas_V))=NaN;
    elseif (r_ind<=N_meas.N_meas_V+N_meas.N_meas_Pinj+N_meas.N_meas_Qinj+N_meas.N_meas_Pij+N_meas.N_meas_Pji)
        meas.LineP_TF(ind_meas.ind_meas_Pji(r_ind-N_meas.N_meas_Pinj-N_meas.N_meas_Qinj-N_meas.N_meas_V-N_meas.N_meas_Pij))=NaN;
    elseif (r_ind<=N_meas.N_meas_V+N_meas.N_meas_Pinj+N_meas.N_meas_Qinj+N_meas.N_meas_Pij+N_meas.N_meas_Pji+N_meas.N_meas_Qij)
        meas.LineQ_FT(ind_meas.ind_meas_Qij(r_ind-N_meas.N_meas_Pinj-N_meas.N_meas_Qinj-N_meas.N_meas_V-N_meas.N_meas_Pij-N_meas.N_meas_Pji))=NaN;
    else
        meas.LineQ_TF(ind_meas.ind_meas_Qji(r_ind-N_meas.N_meas_Pinj-N_meas.N_meas_Qinj-N_meas.N_meas_V-N_meas.N_meas_Pij-N_meas.N_meas_Pji-N_meas.N_meas_Qij))=NaN;
    end
    
        %recomputing Jacobian sparsity pattern
    [topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2017(topo, meas);

    %reconstructing vector z
    z=[meas.BusV(ind_meas.ind_meas_V); meas.BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.BusInjQ(ind_meas.ind_meas_Qinj); meas.LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.LineP_TF(ind_meas.ind_meas_Pji); meas.LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.LineQ_TF(ind_meas.ind_meas_Qji)];
    N_total=length(z);

    %reconstructing the matrix of weights W
    W_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
        meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
        meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
        meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
    W=sparse(1:N_total,1:N_total,1./((W_diag).^2),N_total,N_total);
    Wsqrt=sparse(1:N_total,1:N_total,1./W_diag,N_total,N_total);

    %selecting the initial point
    V=ones(topo.nBus,1);
    theta=zeros(topo.nBus,1);

    %resolving SE and recording the new solution
    [ V_SE1, theta_SE1, eps, ~, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
        ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
    if (convergence==1)
        fprintf('Bad data is corrected.\n\n')
    else
        fprintf('No convergence after removal of bad data.\n\n')
    end
end


%plotting V and theta before and after bad data correction in one plot
figure(y+1);
subplot(2,1,1)  
plot(1:topo.nBus, V_SE0, 1:topo.nBus, V_SE1);
   legend('With error', 'Without error');
title('Voltage Magnitudes in p.u.')
subplot(2,1,2)   
plot(1:topo.nBus, theta_SE0*180/pi, 1:topo.nBus, theta_SE1*180/pi);
title('Voltage Angles in degrees')
y=y+2;
end


