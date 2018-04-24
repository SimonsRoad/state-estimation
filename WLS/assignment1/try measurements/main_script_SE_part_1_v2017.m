clc;clear;
%LOADING THE DATA
load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-5;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=1;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid, %0 - all of them at once

%CHOOSE WHAT TASKS YOU WANT TO COMPLETE (1 if you want, 0 otherwise)
Task1=1; %solve SE with GN using different linear solvers and sparse/dense H
Task2=0; %compute condition numbers and density factors of matrices
Task3=0; %compare the results of solving SE with full/decoupled GN


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

%constructing the matrix of weights W using standard deviations provided in structure 'meas'
%NOTE: matrix W must be constructed in sparse format!!!
%STUDENT CODE 2
W_diag=[meas.std_BusV(ind_meas.ind_meas_V); meas.std_BusInjP(ind_meas.ind_meas_Pinj); ...
    meas.std_BusInjQ(ind_meas.ind_meas_Qinj); meas.std_LineP_FT(ind_meas.ind_meas_Pij); ...
    meas.std_LineP_TF(ind_meas.ind_meas_Pji); meas.std_LineQ_FT(ind_meas.ind_meas_Qij); ...
    meas.std_LineQ_TF(ind_meas.ind_meas_Qji)];
W=sparse(1:N_total,1:N_total,1./((W_diag).^2),N_total,N_total);
Wsqrt=sparse(1:N_total,1:N_total,1./(W_diag),N_total,N_total);

%choosing the initial point
%STUDENT CODE 3
V=ones(topo.nBus,1);
theta=zeros(topo.nBus,1);


%TASK 1: SOLVING SE PROBLEM WITH GN ALGORITHM
%STUDENT CODE 4
if (Task1==1) %if this task is selected for completion
    fprintf('----------------------------\n');
    fprintf('COMPLETING TASK 1\n')
    fprintf('----------------------------\n');
    %preparing to print information about solution
    H_sparsity={'dense' 'sparse'}; H_form={'Full' 'Decoupled'};
    fprintf('%s matrix H is used and it is treated as a %s matrix.\n\n',H_form{H_decoupled+1},H_sparsity{H_sparse+1})
    if (linsolver~=0) %if just one linear solver has been chosen
        linsolver_name={'Matrix inverse' 'Cholesky' 'QR' 'Hybrid'};
        fprintf('For solving a linear system, the following method is used: %s.\n\n',linsolver_name{linsolver})
        %calling the GN algorithm
        [ V, theta, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
            ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
        if (convergence==1)
            fprintf('GN converged in %d iterations and total time for solving linear system is %.4f sec.\n', numel(eps_all), time)
        else
            fprintf('GN did not converge in %d iterations. At the last iteration, max(|delta_x|)=%f.\n',numel(eps_all), eps_all(end))
        end
    else %if the results for all for four methods are desired
        linsolver_name={'Matrix inverse' 'Cholesky      ' 'QR            ' 'Hybrid        '};
        fprintf('The following results were obtained while using different solvers for a linear system:\n\n')
        for linsolver1=1:4
            %calling the GN algorithm
            [ ~, ~, eps_all, time, convergence ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
                ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver1 );
            if (convergence==1)
                fprintf('%s - GN converged in %d iterations and total time for solving linear system is %.4f sec.\n', linsolver_name{linsolver1}, numel(eps_all), time)
            else
                fprintf('%s - GN did not converge in %d iterations. At the last iteration, max(|delta_x|)=%f.\n', linsolver_name{linsolver1}, numel(eps_all), eps_all(end))
            end
        end
    end
    fprintf('----------------------------\n\n\n\n');
end


%TASK 2: COMPUTING CONDITION NUMBERS AND DENSITY FACTORS OF MATRICES
%STUDENT CODE 5
if (Task2==1) %if this task is selected for completion
    fprintf('----------------------------\n');
    fprintf('COMPLETING TASK 2\n')
    fprintf('----------------------------\n');
    %selecting flat start as staring point
    V=ones(topo.nBus,1);
    theta=zeros(topo.nBus,1);
    %making sure we will work with sparse matrices
    H_sparse1=1;  %otherwise Cholesky and QR will fail to produce factorization with minimum number of fill-ins
    %computing the values of the measurement Jacobian
    [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse1);
    
    %computing the density factor of the measurement Jacobian (in percent)
    density.H=100*nnz(H)/numel(H);

    %computing the density factor of matrix G
    G=H'*W*H;
    density.G=100*nnz(G)/numel(G);
    
    %computing condition number of matrix G
    condition_number.G=cond(full(G));
    
    %computing density factor of matrix inverse
    G_inv=inv(G);
    density.G_inv=100*nnz(G_inv)/numel(G_inv);
    
    %using Cholesky and measuring density factor of L
    [L,~,~] = chol(G);
    density.L=100*nnz(L)/numel(L);
    
    %using QR decomposition and measuring density factor of R and Q, and condition number of Htilde
    Htilde = Wsqrt*H;
    [Q,R,e] = qr(Htilde,0);
    density.QR_R=100*(nnz(R))/(numel(R));
    density.QR_Q=100*(nnz(Q))/(numel(Q));
    condition_number.Htilde=cond(full(R)); %condition numbers are the same for R and Htilde
    
    %printing the results
    fprintf('Density factors of matrices are as follows:\n  H      - %3.2f%% \n  G      - %3.2f%% \n  G_inv  - %3.2f%% \n  L_chol - %3.2f%% \n  R      - %3.2f%% \n  Q      - %3.2f%% \n\n', ...
        density.H, density.G, density.G_inv, density.L, density.QR_R, density.QR_Q)
    fprintf('Condition numbers of matrices are as follows:\n  Htilde - %.2e \n  G      - %.2e \n', condition_number.Htilde, condition_number.G)
    fprintf('----------------------------\n\n\n\n');
end


%TASK 3: COMPARING PERFORMANCE OF FULL/DECOUPLED GN
%STUDENT CODE 6
if (Task3==1) %if this task is selected for completion
    fprintf('----------------------------\n');
    fprintf('COMPLETING TASK 3\n')
    fprintf('----------------------------\n');
    %selecting flat start as staring point
    V=ones(topo.nBus,1);
    theta=zeros(topo.nBus,1);
    %making sure Cholesky factorization is used for full GN (so that the comparison is fair)
    linsolver=2;
    
    %solving SE with full GN
    H_decoupled=0;
    [ ~, ~, Eps.Full_NR, Time.Full_NR, Convergence.Full_NR ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
        ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
    if (Convergence.Full_NR==1)
        fprintf('Full GN converged in %d iterations and total time for solving linear system is %.4f sec.\n\n', numel(Eps.Full_NR), Time.Full_NR)
    else
        fprintf('Full GN did not converge in %d iterations. At the last iteration, max(|delta_x|)=%f.\n\n',numel(Eps.Full_NR), Eps.Full_NR(end))
    end
    
    %solving SE with fast decoupled GN
    H_decoupled=1;
    [ ~, ~, Eps.Decoupled_NR, Time.Decoupled_NR, Convergence.Decoupled_NR ] = f_SE_NR_algorithm_v2017 ( V, theta, topo, Y_bus, z, W, Wsqrt, ...
        ind_meas, N_meas, eps_tol, Max_iter, H_decoupled, H_sparse, linsolver );
    if (Convergence.Decoupled_NR==1)
        fprintf('Decoupled GN converged in %d iterations and total time for solving linear system is %.4f sec.\n', numel(Eps.Decoupled_NR), Time.Decoupled_NR)
    else
        fprintf('Decoupled GN did not converge in %d iterations. At the last iteration, max(|delta_x|)=%f.\n',numel(Eps.Decoupled_NR), Eps.Decoupled_NR(end))
    end
    
    %comparing and printing the results
    if (Convergence.Decoupled_NR==1 && Convergence.Full_NR==1) %comparison only makes sense if both algorithms have converged
        if (Time.Decoupled_NR>Time.Full_NR) %if decoupled GN is slower than full GN
            Time.Decoupled_to_Full=Time.Decoupled_NR/Time.Full_NR; %computing the slow-down factor
            fprintf('\nDecoupled GN is %.1f times slower than full GN.\n\n', Time.Decoupled_to_Full)
        else %if decoupled GN is faster than full GN
            Time.Full_to_Decoupled=Time.Full_NR/Time.Decoupled_NR; %computing the speed-up factor
            fprintf('\nDecoupled GN is %.1f times faster than full GN.\n\n', Time.Full_to_Decoupled)
        end
        %computing the speed-up factor of each iteration of decoupled GN over
        Time.Decoupled_to_Full_one_iteration=(Time.Full_NR/numel(Eps.Full_NR))/(Time.Decoupled_NR/numel(Eps.Decoupled_NR));
        fprintf('Each iteration of decoupled GN is %.1f times faster than of full GN.\n', Time.Decoupled_to_Full_one_iteration)
        
        %plotting the graph of log(eps_all)=f(iterations). Using log is better because otherwise the graph will be hard to read
        plot(1:numel(Eps.Decoupled_NR),log(Eps.Decoupled_NR),1:numel(Eps.Full_NR),log(Eps.Full_NR), ...
            1:max(numel(Eps.Decoupled_NR),numel(Eps.Full_NR)), log(ones(max(numel(Eps.Decoupled_NR),numel(Eps.Full_NR)),1)*eps_tol))
        title('Convergence Pattern')
        legend('Decoupled GN', 'Full GN', 'Set tolerance');
        xlabel('Iterations');
        ylabel('log(max(|\DeltaX|))');
    end
    fprintf('----------------------------\n\n');
end






