clc;clear;
%LOADING DATASET
load Dataset_B.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm (using default options is recommended)
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %use full Jacobian matrix
H_sparse=1; %use matrix in a sparse format
linsolver=2;  %use Cholesky


%% DATA PREPARATION 
%getting Jacobian sparsity pattern (done only once)

%deleting the NaN rows in meas
%V
t = any(isnan(meas.BusV),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.BusV(v,:) = [];
meas.std_BusV(v,:) = [];
%Pinj
t = any(isnan(meas.BusInjP),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.BusInjP(v,:) = [];
meas.std_BusInjP(v,:) = [];
%Qinj
t = any(isnan(meas.BusInjQ),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.BusInjQ(v,:) = [];
meas.std_BusInjQ(v,:) = [];
%LinePFT
t = any(isnan(meas.LineP_FT),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.LineP_FT(v,:) = [];
meas.std_LineP_FT(v,:) = [];
%LinePTF
t = any(isnan(meas.LineP_TF),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.LineP_TF(v,:) = [];
meas.std_LineP_TF(v,:) = [];
%LineQFT
t = any(isnan(meas.LineQ_FT),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.LineQ_FT(v,:) = [];
meas.std_LineQ_FT(v,:) = [];
%LineQTF
t = any(isnan(meas.LineQ_TF),2);
[i,~] = ind2sub(size(t),find(t));
v = unique(i);
v = v(:);
meas.LineQ_TF(v,:) = [];
meas.std_LineQ_TF(v,:) = [];

x=1;
y=0;
while x==1

[topo, ind_meas, N_meas]=f_meas_indices_and_H_sparsity_pattern_v2017(topo, meas);

%obtaining the admittance matrix in the compressed row form (done only once)
Y_bus = f_Y_bus_compressed_v2017(topo);



%constructing vector z, which is a vector of available measurements
%NOTE: different types of measurements should be in the same order as in the handout
%STUDENT CODE 1
Com=[meas.BusV,meas.std_BusV;
     meas.BusInjP,meas.std_BusInjP;
     meas.BusInjQ,meas.std_BusInjQ;
     meas.LineP_FT,meas.std_LineP_FT;
     meas.LineP_TF,meas.std_LineP_TF;
     meas.LineQ_FT,meas.std_LineQ_FT;
     meas.LineQ_TF,meas.std_LineQ_TF;
    ];
Com = Com(all(~isnan(Com),2),:); % delete nan rows
z = Com(:,1);


%constructing the matrix of weights W and its square root Wsqrt using standard deviations provided in structure 'meas'
%NOTE: matrices W and Wsqrt must be constructed in sparse format!!!
%STUDENT CODE 2
W=diag(Com(:,2).^(-2));
W=sparse(W);
Wsqrt=diag(Com(:,2).^(-1));




%% INITIAL SE SOLUTION
fprintf('INITIAL STATE ESTIMATION: \n')
%choosing the initial point
%STUDENT CODE 3
a=size(topo.busNumbers,1);
V=ones(a,1);
theta=zeros(a,1);

%calling the GN algorithm and recording initial solution of SE (for
%subsequent plotting of voltage magnitudes and angles)
[ V_SE0, theta_SE0, ~, ~, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );



%% BAD DATA DETECTION
%STUDENT CODE 4
[ H ] = f_measJac_H_v2017( V_SE0, theta_SE0, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
[ h ] = f_measFunc_h_v2017( V_SE0, theta_SE0, Y_bus, topo, ind_meas, N_meas);
Htran=H.';
G=Htran*W*H;
K=(H/G)*Htran*W;
S=1.-K;
omega=S/W;
num_r=size(omega,1);
r=z-h;
omega_s=sqrt(abs(diag(omega)));
r_N=abs(r)./omega_s;
r_N(~isfinite(r_N))=0;
%%detect bad data 
i=(1:num_r)';
%scatter(i,r);
threshold=4;
x=0;
for i=1:num_r
    if r_N(i)>threshold
        x=1;
    else
        x=0;
    end
end
if x==1
    disp('have bad data')
else
    disp('no bad data')
end
%% BAD DATA CORRECTION
%STUDENT CODE 5
maxmea=0;
if x==1
   [val, idx] = max(r_N);
end

%deleting the related measurements in mea
a1=size(meas.BusV,1);
a2=size(meas.BusV,1)+size(meas.BusInjP,1);
a3=size(meas.BusV,1)+size(meas.BusInjP,1)+size(meas.BusInjQ,1);
a4=size(meas.BusV,1)+size(meas.BusInjP,1)+size(meas.BusInjQ,1)+size(meas.LineP_FT,1);
a5=size(meas.BusV,1)+size(meas.BusInjP,1)+size(meas.BusInjQ,1)+size(meas.LineP_FT,1)+size(meas.LineP_TF,1);
a6=size(meas.BusV,1)+size(meas.BusInjP,1)+size(meas.BusInjQ,1)+size(meas.LineP_FT,1)+size(meas.LineP_TF,1)+size(meas.LineQ_FT,1);
a7=size(meas.BusV,1)+size(meas.BusInjP,1)+size(meas.BusInjQ,1)+size(meas.LineP_FT,1)+size(meas.LineP_TF,1)+size(meas.LineQ_FT,1)+size(meas.LineQ_TF,1);

%BusV
if idx <= a1
    d=idx;
    meas.BusV(d,:) = [];
    meas.std_BusV(d,:) = [];
%BusInjP    
elseif idx>a1 && idx<=a2
    d=idx-a1;
    meas.BusInjP(d,:) = [];
    meas.std_BusInjP(d,:) = [];
%BusInjQ    
elseif idx>a2 && idx<=a3
    d=idx-a2;
    meas.BusInjQ(d,:) = [];
    meas.std_BusInjQ(d,:) = [];
%LineP_FT    
elseif idx>a3 && idx<=a4
    d=idx-a3;
    meas.LineP_FT(d,:) = [];
    meas.std_LineP_FT(d,:) = [];
%LineP_TF    
elseif idx>a4 && idx<=a5
    d=idx-a4;
    meas.LineP_TF(d,:) = [];
    meas.std_LineP_TF(d,:) = [];
elseif idx>a5 && idx<=a6
    d=idx-a5;
    meas.LineQ_FT(d,:) = [];
    meas.std_LineQ_FT(d,:) = [];
elseif idx>a6 && idx<=a7
    d=idx-a6;
    meas.LineQ_TF(d,:) = [];
    meas.std_LineQ_TF(d,:) = [];    
end
%iteration times
y=y+1;
end


