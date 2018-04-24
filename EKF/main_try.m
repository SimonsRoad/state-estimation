clc;clear;
%LOADING THE DATA
load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;
%load Dataset_B.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=4;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid
alpha=0.775;%smoothing parameter
beta=0.1;%smoothing parameter



%initial start point
Pe=0.1*eye(topo.nBus*2-1,topo.nBus*2-1);
Q=10^(-4)*eye(topo.nBus*2-1,topo.nBus*2-1);
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c,1);
Vf=ones(c,1);
thetaf=zeros(c-1,1);
a=zeros(c*2-1,1);
b=zeros(c*2-1,1);
i=0;
Errortype=1;%make sure go into the loop


%SOLTUION ALGORITHM STARTS HERE
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

%constructing the matrix of measurements covariance matrix R and state vector covariance matrix Pe using standard deviations provided in structure 'meas'
%NOTE: matrices R must be constructed in sparse format!!!
%STUDENT CODE 2
R_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
    meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
    meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
    meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
R=sparse(1:N_total,1:N_total,((R_diag).^2),N_total,N_total);
W=sparse(1:N_total,1:N_total,1./((R_diag).^2),N_total,N_total);

%build a loop to reduce the gap between the initial estimation and real
%measremnts(due to the unprecise inital points)
while (i<30)
%TASK 1: forcasted state vectors and measurements
%Vf 14*1, thetaf 13*1
[ Vf, thetaf, a, b, zf, Tf, H] = forcast (a, b, alpha, beta, Vf, thetaf, Ve, thetae, Pe, Q,...
    topo, Y_bus, ind_meas, N_meas, H_decoupled, H_sparse );
xf = [thetaf;Vf];

%TASK 2: state correction
[ Ve, thetae, xe] = f_SE_NR_algorithm_v2017 (topo, z, zf, H, R, xf);
i=i+1;
end



%TASK 3: anomaly detection
%STUDENT CODE 6
while Errortype==1
    fprintf('SEARCHING FOR BAD DATA: \n')
    %getting most recent H and h
    ze=f_measFunc_h_v2017( Ve, thetae, Y_bus, topo, ind_meas, N_meas);%estimated measurements
    [ H ] = f_measJac_H_v2017( Ve, thetae, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);

    %innovation analysis
    N=R+Tf;%innovation error covariance matrix
    omega_f=sqrt(abs(diag(N)));
    v = z-zf;%innovation vector
    v_N=abs(v)./omega_f;


    %residual analysis
    % identifying critical measurements    
    r=z-ze;%residual
    %computing the required matrices
    G = H'*W*H;
    K = H*(G\(H'*W));
    S = (eye(N_total) - K);
    S_diag = diag(S);
    W_diag = diag(W);
    % creating vector for normalized residuals
    r_N = zeros(N_total,1);
    i_critical = abs(diag(S)) < 1e-12;
    r_N(~i_critical) = abs(r(~i_critical)) ./ sqrt(S_diag(~i_critical)./W_diag(~i_critical));
    %getting the value and index of largest normalized residual
    [r_max,r_ind] = max(r_N);
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
        R_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
                meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
                meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
                meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
        R=sparse(1:N_total,1:N_total,((R_diag).^2),N_total,N_total);
        W=sparse(1:N_total,1:N_total,1./((R_diag).^2),N_total,N_total);

        %TASK 2: state correction
        [ Vf, thetaf, a, b, zf, Tf, H] = forcast (a, b, alpha, beta, Vf, thetaf, Ve, thetae, Pe, Q,...
        topo, Y_bus, ind_meas, N_meas, H_decoupled, H_sparse );
        v=z-zf;
        Pe = inv(H.'/R*H);%state vector covariance matrix at time k+1 
        K = Pe*H.'/R;
        xf = [thetaf;Vf];
        xe = xf+K*v;%xe is 27*1 here
        thetae=[0;xe(1:topo.nBus-1,1)];%add reference theta inside
        Ve=xe(topo.nBus:topo.nBus*2-1,1);
    end
    
end

