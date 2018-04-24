function [ V, theta, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver )
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
if (H_decoupled==0) %if full Jacobian is used
    %Starting GN loop
    while iter<=Max_iter && eps>eps_tol
        %getting measurement Jacobian
        [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        %getting vector of operating point parameters
        [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
        %choosing linear solver
        tic;
        switch linsolver
            case 1  %matrix inverse of G
                G=H'*W*H;
                rhs=H'*W*(z-h);
                Delta_X=inv(G)*rhs;
            case 2  %Cholesky factorization of G
                G=H'*W*H;
                rhs=H'*W*(z-h);
                Delta_X=G\rhs; %backslash in this case will use Cholesky
            case 3  %QR decomposition with minimum number of fill-ins
                Htilde =Wsqrt*H;
                [Q,R,e] = qr(Htilde,0); %e is a permutation vector of columns of Htilde, used to minimize the number of nonzeros in R and Q
                rhs = Q'*Wsqrt*(z-h);
                Delta_X=R\rhs;
                Delta_X(e)=Delta_X; %permute vector of updates back to original
            case 4  %Hybrid method (with preordering the columns of H)
                Htilde =Wsqrt*H;
                p = colamd(Htilde); %getting the permutation vector of columns of Htilde to minimize the number of fill-ins in R
                R=qr(Htilde(:,p),0); %doing QR without computing Q matrix -> makes QR much faster
                rhs = Htilde(:,p)'*Wsqrt*(z-h);
                Delta_X=R\(R'\rhs);
                Delta_X(p)=Delta_X; %permute vector of updates back to original
        end
        %checking stopping criterion
        eps=max(abs(Delta_X));
        %recording maximum change in delta_x at the current iteration
        eps_all(iter,1)=eps;
        %updating variables
        theta(2:topo.nBus)=theta(2:topo.nBus)+Delta_X(1:topo.nBus-1);
        V=V+Delta_X(topo.nBus:end);
        time=time+toc; %recording the time
        %updating iteration counter
        iter=iter+1;
    end
    
else %if decoupled Jacobian is used (NOTE: in this case, we'll only use Cholesky factorization)
    %THIS CODE WILL NOT BE DISTRIBUTED
    eps_all=ones(Max_iter,1);
    iter=Max_iter+1;
    time=inf;
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

