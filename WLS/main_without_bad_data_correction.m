clc;clear;
%LOADING THE DATA
%load GoodMeasurement_14_bus.mat;
%load GoodMeasurement_1354_bus.mat;
load Measurement.mat

%HERE COME USER-DEFINED PARAMETERS OF GN algorithm
eps_tol=10^-6;  %stopping criterion for max(abs(delta_x))
Max_iter=100; %maximum number of GN iterations
H_decoupled=0; %if 0, full H and 'normal' GN are used. If 1, decoupled H and fast decoupled GN are used
H_sparse=1; %if 1, H is created as a sparse matrix, if 0 - as a dense matrix
linsolver=2;  %1 - matrix inverse, 2 - Cholesky, 3 - QR, 4 - Hybrid, %0 - all of them at once

%SOLTUION ALGORITHM STARTS HERE
a=size(meas.BusV,2);%getting the total number of measurements
meas_new=meas;
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
    state(:,i)=[theta;V];
    [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas);
    J(i)=(z-h)'*W*(z-h);
end





