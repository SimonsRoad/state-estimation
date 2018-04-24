clear all;
%LOADING THE DATA
load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;
%load Dataset_B.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=20; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=1;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid
c1=2.7;

%initial start point
c=size(topo.busNumbers,1);
Ve=ones(c,1);
thetae=zeros(c,1);
iter=0;
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
b=9*ones(N_total,1);%should be changed to checked table

%NOTE: matrices R must be constructed in sparse format!!!
%STUDENT CODE 2
R_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
    meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
    meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
    meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
R=sparse(1:N_total,1:N_total,((R_diag).^2),N_total,N_total);
R_s=sparse(1:N_total,1:N_total,(((R_diag).^(2)).^(-0.5)),N_total,N_total);%R_s=R^(-1/2)
W=sparse(1:N_total,1:N_total,1./((R_diag).^2),N_total,N_total);
Wsqrt=sparse(1:N_total,1:N_total,1./(R_diag),N_total,N_total);

eps=eps_tol+1; %making sure we'll enter the loop
while iter<=Max_iter && eps>eps_tol
%TASK 1: getting H(x0):
 H = f_measJac_H_v2017( Ve, thetae, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);

%TASK 2: calculating PS:

L=R_s*H;%L=[l1,l2,...lm]T; example:l1=L(1,:)'
omega=full(L*L');

if N_total>9
    if rem(N_total,2)==1
        fm=N_total/(N_total-0.9);
    else
        fm=1;
    end
else 
    fm=1;%simplify fm, if N_total is smaller than 9, fm should be determined based on the table in the paper
end    

for k=1:N_total
    for i=1:N_total
        for j=1:N_total
            if(i==j)
                %do nothing
            else
                x(j)=abs(omega(i,k)+omega(j,k));
            end
        end
        x(x==0) = [];%deleting the 0 in x, not sure if it is correct to do so?? 
        y(i)=lomed(x);
    end
    sm(k)=1.1926*fm*lomed(y);
end

for i=1:N_total
    PS(i)=max(omega(i,:)./sm);
end

%calculate w(i)
for i=1:N_total
    if (b(i)/PS(i))^2<1
        if (b(i)/PS(i))^2<0.001
            w(i)=0.001;
        else
            w(i)=(b(i)/PS(i))^2;
        end
    else
        w(i)=1;
    end
end

h=f_measFunc_h_v2017( Ve, thetae, Y_bus, topo, ind_meas, N_meas);
r=z-h;

for i=1:N_total
    rS(i)=r(i)/(R_diag(i)*w(i));
    if abs(rS(i)<=c1)
        phi(i)=rS(i);
    else
        phi(i)=c1*sign(rS(i));
    end    
end    
q=phi./rS;
q(isnan(q))=1;%when rS(i)=0, q(rS(i))=rS(i)/rS(i)=1
Q=diag(q);
G=H'*W*Q*H;
rhs=H'*W*Q*(z-h);
Delta_X=G\rhs;
eps=max(abs(Delta_X));
iter=iter+1;
thetae(2:topo.nBus)=thetae(2:topo.nBus)+Delta_X(1:topo.nBus-1);
Ve=Ve+Delta_X(topo.nBus:end);
end


%TASK 2: state correction


