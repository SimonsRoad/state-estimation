function [ Y_bus ] = f_Y_bus_compressed_v2017(topo)
%% ************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function calculates admittance matrix in a compressed form (i.e. only
%structurally nonzero elements are computed). The order of the elements in
%the output vector corresponds to the order, in which these elements are 
%used for computing the Jacobian of power balance equations: for node i,
%first, admittances of all adjacent nodes y_ij are recorded and, finally, 
%admittance y_ii is recorded. Storing admittance elements in this order
%helps speed up the calculation of the Jacobian of power balance equations.

%INPUTS
%topo     - structure containing all information about the system
%           (topology, parameters of the elements, etc)

%OUTPUTS
%Y_bus    - structure that contains real and imaginary parts of admittances
%           of structurally nonzero elements in the admittance matrix.
%           NOTE: G and B in this case are vectors, not matrices.


%% ************************************************************************
%INITIALIZATION OF VARIABLES
%**************************************************************************
%extracting data from the structure to speed up subsequent calculations
N=topo.nBus; L=topo.nBranch; branch_R=topo.lineR; branch_X=topo.lineX;
node_B_sh=topo.busShuntB; node_G_sh=topo.busShuntG; branch_Bc=topo.lineShuntB;
branch_Kt_real=topo.internal.lineTapRatioReal; branch_Kt_imag=topo.internal.lineTapRatioImag;
number_nonzeros_Y=topo.internal.number_nonzeros_Y; branch_info_1=topo.internal.branch_info_1;
branch_info_2=topo.internal.branch_info_2; row_start_indexes=topo.internal.row_start_indexes;
number_connected_nodes=topo.internal.number_connected_nodes; 
node_all_connections=topo.internal.node_all_connections;
%initializing the vectors of active and reactive parts of admittance elements
G=zeros(number_nonzeros_Y,1);
B=zeros(number_nonzeros_Y,1);


%% ************************************************************************
%COMPUTATION
%**************************************************************************
for i=1:L
    R = branch_R(i);
    X = branch_X(i);
    K1 = branch_Kt_real(i);
    K2 = branch_Kt_imag(i);
    %getting the positions of the first and second nodes of this branch in the vector of nodes;
    i_begin = branch_info_1(i);
    i_end = branch_info_2(i);
    
    %getting position of second node in a list of all nodes connected to first node
    i_start=row_start_indexes(i_begin);
    i_finish=i_start+number_connected_nodes(i_begin)-1;
    for ii=i_start:i_finish
        if (node_all_connections(ii)==i_end)
            index_end=ii;
            break;
        end
    end
    
    %getting position of first node in a list of all nodes connected to second node
    i_start=row_start_indexes(i_end);
    i_finish=i_start+number_connected_nodes(i_end)-1;
    for ii=i_start:i_finish
        if (node_all_connections(ii)==i_begin)
            index_begin=ii;
            break;
        end
    end
    
    %precomputing some values for simplicity
    buf_R = R / (R*R + X*X);
    buf_X = -X / (R*R + X*X);
    
    if (K1 == 0 && K2 == 0)   %if this branch is not a transformer
        %for Y_ij
        index_buf = index_begin;
        G(index_buf) = G(index_buf) - buf_R;
        B(index_buf) = B(index_buf) - buf_X;
        
        %for Y_ji
        index_buf = index_end;
        G(index_buf) = G(index_buf) - buf_R;
        B(index_buf) = B(index_buf) - buf_X;
          
        %for Y_ii
        index_buf = row_start_indexes(i_begin) + number_connected_nodes(i_begin);
        G(index_buf) = G(index_buf) + buf_R;
        B(index_buf) = B(index_buf) + buf_X - branch_Bc(i) * 0.5;
        
        %for Y_jj
        index_buf = row_start_indexes(i_end) + number_connected_nodes(i_end);
        G(index_buf) = G(index_buf) + buf_R;
        B(index_buf) = B(index_buf) + buf_X - branch_Bc(i) * 0.5;
    else
        %for Y_ij
        index_buf = index_begin;
        G(index_buf) = G(index_buf) - (R*K1 + X*K2) / ((R*R + X*X)*(K1*K1 + K2*K2));
        B(index_buf) = B(index_buf) + (-R*K2 + X*K1) / ((R*R + X*X)*(K1*K1 + K2*K2));
        
        %for Y_ji
        index_buf = index_end;
        G(index_buf) = G(index_buf) - (R*K1 - X*K2) / ((R*R + X*X)*(K1*K1 + K2*K2));
        B(index_buf) = B(index_buf) + (R*K2 + X*K1) / ((R*R + X*X)*(K1*K1 + K2*K2));
        
        %for Y_ii
        index_buf = row_start_indexes(i_begin) + number_connected_nodes(i_begin);
        G(index_buf) = G(index_buf) + buf_R / (K1*K1 + K2*K2);
        B(index_buf) = B(index_buf) + (buf_X - branch_Bc(i)*0.5) / (K1*K1 + K2*K2);
        
        %for Y_jj
        index_buf = row_start_indexes(i_end) + number_connected_nodes(i_end);
        G(index_buf) = G(index_buf) + buf_R;
        B(index_buf) = B(index_buf) + buf_X - branch_Bc(i) * 0.5;
    end
end

%including node shunts
for i=1:N
    index_buf = row_start_indexes(i) + number_connected_nodes(i);
    G(index_buf) = G(index_buf) + node_G_sh(i);
    B(index_buf) = B(index_buf) - node_B_sh(i);
end

%saving the result into a structure
Y_bus.G=G;
Y_bus.B=B;
end