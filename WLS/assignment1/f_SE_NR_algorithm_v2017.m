function [ V, theta, eps_all, time, convergence] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
    ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver )
%**************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function solves the State Estimation problem cast as a WLS problem
%using the Gauss-Newton (GN) method. The possible solution options include
%choosing different algorithms for solving a linear system at each
%iteration as well as using either a full GN method or fast decoupled GN
%method.

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
%eps_tol     - stopping criterion. GN should stop if max(|delta_x|)<eps_tol
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
%theta       - vector of voltage phase angles of ALL buses in the system
%eps_all     - vector, each element of which is equal to max(|delta_x|) at
%              a certain iteration. Shows the convergence pattern
%time        - total computation time for all iterations for solving linear
%              system. Time for computing H and h should NOT be included here
%convergence - if 1, GN has converged to a solution, otherwise 0

%**************************************************************************
%FUNCTION CODE
%**************************************************************************


%STUDENT CODE
%for computing matrix H and vector h, use the following functions:
i=0;
deltax=eps_tol+1;
a=size(topo.busNumbers,1);
tic;
if linsolver==1
    
    while (i<Max_iter) && (max(abs(deltax))>eps_tol)
    
        i=i+1;
        [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
        Htran=H.';
        G=Htran*W*H;
        invG=inv(G);
        second=Htran*W*(z-h);
        deltax=invG*second;
        eps_all(i)=max(abs(deltax));
        deltatheta=deltax(1:a-1);
        theta=theta+[0;deltatheta];
        deltaV=deltax(a:a*2-1);
        V=V+deltaV; 
    
    end
    
elseif linsolver==2
        
    while (i<Max_iter) && (max(abs(deltax))>eps_tol)
    
        i=i+1;
        [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
        Htran=H.';
        G=Htran*W*H;
        
        second=Htran*W*(z-h);
        deltax=G\second;
        eps_all(i)=max(abs(deltax));
        deltatheta=deltax(1:a-1);
        theta=theta+[0;deltatheta];
        deltaV=deltax(a:a*2-1);
        V=V+deltaV;
        
    end

elseif linsolver==3
    while (i<Max_iter) && (max(abs(deltax))>eps_tol)
    
        i=i+1;
        [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
        Hs=Wsqrt*H;
        [Q,R,e]=qr(Hs,0);
        deltax=R\((Q.')*Wsqrt*(z-h));
        deltax(e)=deltax;
        
        deltatheta=deltax(1:a-1);
        theta=theta+[0;deltatheta];
        deltaV=deltax(a:a*2-1);
        V=V+deltaV;
       
        eps_all(i)=max(abs(deltax));
    end
elseif linsolver==4
    while (i<Max_iter) && (max(abs(deltax))>eps_tol)
    
        i=i+1;
        [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse);
        [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
        Hs=Wsqrt*H;
        p=colamd(Hs);
        R=qr(Hs(:,p),0);
        zs=Wsqrt*(z-h);
        zh=(Hs(:,p).')*zs;
        deltax=R\((R.')\zh);
        deltax(p)=deltax;

        deltatheta=deltax(1:a-1);
        theta=theta+[0;deltatheta];
        deltaV=deltax(a:a*2-1);
        V=V+deltaV;
        
        deltaxback(:,i)=deltax;
        eps_all(i)=max(abs(deltax));
    end
end
toc;
time=toc;

if(max(abs(deltax))<eps_tol)
    convergence = 1;
else
    convergence = 0;
end