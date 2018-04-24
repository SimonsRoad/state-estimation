function [ Ve, thetae, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( Ve, thetae, topo, Y_bus, z, W, ...
    R_s, R_diag, b, c1, ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, N_total)
%**************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function solves the State Estimation problem cast as a WLS problem
%using the Newton-Raphson method. The possible solution options include
%choosing different algorithms for solving a linear system at each
%iteration as well as using either a full GN method or fast decoupled GN
%method. NOTE: For more information about this function and reasons behind
%choosing the algorithms that you see here, please refer to the TA manual.

%INPUTS
%V           - vector of voltage magnitudes of ALL buses in the system
%theta       - vector of voltage phase angles for ALL buses in the system
%topo        - structure containing the topology of the system and
%              parameters of system elements. Used for computing H and h
%Y_bus       - structure that contains real and imag parts of admittances
%              of structurally nonzero elements in the admittance matrix
%z           - vector containing all available measurements in the system
%W           - diagonal matrix containing weights assigned to measurements
%Wsqrt       - diagonal matrix, elements are square roots of elements of W   
%ind_meas    - structure that for each type of measurement contains info
%              on what are the indexes of available measurements in a list
%              of all possible measurements for this type
%N_meas      - structure that contains the number of available measurements
%              for each type of measurements 
%eps_tol     - stopping criterion. GN stops if max(|delta_x|)<eps_tol
%Max_iter    - stopping criterion, maximum number of GN iterations
%H_decoupled - control variable. If 0, full Jacobian is constructed and 
%              'full' GN algorithm is used. If 1, decoupled Jacobian is 
%              constructed and fast decoupled algorithm is used
%H_sparse    - control variable. If 1, matrix H will be in a sparse form. 
%              If 0, matrix H will be in a dense form
%linsolver   - control variable for choosing a linear system solver:
%              1 - matrix inverse of G
%              2 - Cholesky factorization of G (done by backslash)
%              3 - QR decomposition of H with minimum number of fill-ins
%              4 - Hybrid method (decompose H by QR, but only compute R)

%OUTPUTS
%V           - vector of voltage magnitudes of ALL buses in the system
%theta       - vector of voltage phase angles for ALL buses in the system
%eps_all     - vector, each element of which is equal to max(|delta_x|) at
%              a certain iteration. Shows the convergence pattern
%time        - total computation time for all iterations for solving linear
%              system. Time for computing H and h is not included here
%convergence - if 1, GN has converged to a solution, otherwise 0

%**************************************************************************
%FUNCTION CODE
%**************************************************************************
eps=eps_tol+1; %making sure we'll enter the loop
eps_all=zeros(Max_iter,1); %initializing eps_all to maximum possible size
iter=1;  %set the initial iteration counter to 1
time=0;  %initialize the simulation time
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
eps_all(iter)=eps;
end

%retaining eps_all only for actually performed iterations
eps_all=eps_all(1:iter-1);

%checking convergence
if (eps<=eps_tol)
    convergence=1;
else
    convergence=0;
end

end