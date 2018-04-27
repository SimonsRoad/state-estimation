function [ H ] = f_measJac_H_v2017( V, theta, Y_bus, topo, ind_meas, N_meas, H_decoupled, H_sparse)
%% ************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function calculates the Jacobian of measurements (H). The measured
%quantities in the system are voltage magnitudes (V), active (Pinj) and 
%reactive (Qinj) power injections at buses, active and reactive line flows
%at both ends of the lines (Pij, Pji, Qij, Qji). The unknown variables are
%voltage magnitudes (V) and phase angles (theta). The Jacobian of
%measurements consists of a number of submatrices, each of which represents
%a partial derivative of a certain type of measured quantity w.r.t. 
%a certain type of unknown variable (e.g. dPinj/dtheta, dQij/dV etc). This
%function computes each submatrix in a separate loop and stores all of
%submtarices in a sparse format. In the end, all submatrices are combined
%into one matrix H. Note: the formulas for the derivatives of line flows
%work only when there are no phase-shifting transformers in the system.

%INPUTS
%V           - vector of voltage magnitudes of ALL buses in the system
%theta       - vector of voltage phase angles for ALL buses in the system
%Y_bus       - structure that contains real and imag parts of admittances
%              of structurally nonzero elements in the admittance matrix
%internal    - structure that contains information on topology of the 
%              system and parameters of system elements
%ind_meas    - structure that for each type of measurement contains info
%              on what are the indexes of available measurements in a list
%              of all possible measurements for this type
%N_meas      - structure that contains the number of available measurements
%              for each type of measurements 
%H_decoupled - control variable. If 0, full Jacobian is constructed. If 1,
%              decoupled Jacobian is constructed
%H_sparse    - control variable. If 1, matrix H will be in a sparse form. 
%              If 0, matrix H will be in a dense form

%OUTPUTS
%H           - Jacobian of measurements, constructed as sparse/dense matrix


%% ************************************************************************
%INITIALIZATION OF VARIABLES
%**************************************************************************
%extracting data from the structure to speed up subsequent calculations
%getting graph connectivity pattern
row_start_indexes=topo.internal.row_start_indexes; 
node_all_connections=topo.internal.node_all_connections;
number_connected_nodes=topo.internal.number_connected_nodes;
branch_index_of_1_node=topo.internal.branch_info_1;
branch_index_of_2_node=topo.internal.branch_info_2;
%getting the Jacobian sparsity pattern
row_indexes=topo.internal.row_indexes_H; column_indexes=topo.internal.column_indexes_H;
number_nonzeros=topo.internal.number_nonzeros_H;
%getting info about which measurements are avaiable
ind_meas_V=ind_meas.ind_meas_V;        N_meas_V=N_meas.N_meas_V;
ind_meas_Pinj=ind_meas.ind_meas_Pinj;  N_meas_Pinj=N_meas.N_meas_Pinj;
ind_meas_Qinj=ind_meas.ind_meas_Qinj;  N_meas_Qinj=N_meas.N_meas_Qinj;
ind_meas_Pij=ind_meas.ind_meas_Pij;    N_meas_Pij=N_meas.N_meas_Pij;
ind_meas_Pji=ind_meas.ind_meas_Pji;    N_meas_Pji=N_meas.N_meas_Pji;
ind_meas_Qij=ind_meas.ind_meas_Qij;    N_meas_Qij=N_meas.N_meas_Qij;
ind_meas_Qji=ind_meas.ind_meas_Qji;    N_meas_Qji=N_meas.N_meas_Qji;
%getting parameters of system elements (admittances etc.)
G_values=Y_bus.G; B_values=Y_bus.B; IsTransformer=topo.internal.IsTransformer;
branch_Kt_real=topo.internal.lineTapRatioReal; branch_g=topo.internal.branch_g; 
branch_b=topo.internal.branch_b; branch_b_sh=topo.lineShuntB; 
N=topo.nBus; N_bus=N-1;




%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dV/dtheta
%**************************************************************************
dV_dtheta = sparse(N_meas_V,N_bus); 


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dV/dV
%**************************************************************************
dV_dV = sparse(1:N_meas_V,ind_meas_V, ones(N_meas_V,1), N_meas_V, N); 


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPinj/dtheta
%**************************************************************************
J_v=zeros(number_nonzeros.dPinj_dtheta,1);
nonzeros=0;
for i0=1:N_meas_Pinj
    i=ind_meas_Pinj(i0);
    sum1=0;
    N1=number_connected_nodes(i);
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        sum0=-G_values(index_y)*sin(theta(i)-theta(j))+B_values(index_y)*cos(theta(i)-theta(j));
        sum1=sum1+V(j)*sum0;
        if (j~=1)
            nonzeros=nonzeros+1;
            J_v(nonzeros)=-V(j)*V(i)*sum0; %dWp_i/dtheta_j
        end
    end
    if (i~=1)
        nonzeros=nonzeros+1;
        J_v(nonzeros)=sum1*V(i);  %dWp_i/dtheta_i
    end
end
dPinj_dtheta = sparse(row_indexes.dPinj_dtheta,column_indexes.dPinj_dtheta, J_v, N_meas_Pinj, N_bus); %derivatives of P_inj w.r.t theta


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPinj/dV
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dPinj_dV,1);
    nonzeros=0;
    for i0=1:N_meas_Pinj
        i=ind_meas_Pinj(i0);
        sum3=0;
        N1=number_connected_nodes(i);
        index_y=row_start_indexes(i)-1;
        for i1=1:N1
            index_y=index_y+1;
            j=node_all_connections(index_y);
            sum2=G_values(index_y)*cos(theta(i)-theta(j))+B_values(index_y)*sin(theta(i)-theta(j));
            sum3=sum3+V(j)*sum2;
            nonzeros=nonzeros+1;
            J_v(nonzeros)=V(i)*sum2; %dWp_i/dV_j
        end
        index_y=index_y+1;
        nonzeros=nonzeros+1;
        J_v(nonzeros)=2*G_values(index_y)*V(i)+sum3; %dWp_i/dV_i
    end
    dPinj_dV = sparse(row_indexes.dPinj_dV,column_indexes.dPinj_dV, J_v, N_meas_Pinj, N); %derivatives of P_inj w.r.t V
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQinj/dtheta
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dQinj_dtheta,1);
    nonzeros=0;
    for i0=1:N_meas_Qinj
        i=ind_meas_Qinj(i0);
        sum3=0;
        N1=number_connected_nodes(i);
        index_y=row_start_indexes(i)-1;
        for i1=1:N1
            index_y=index_y+1;
            j=node_all_connections(index_y);
            sum2=G_values(index_y)*cos(theta(i)-theta(j))+B_values(index_y)*sin(theta(i)-theta(j));
            sum3=sum3+V(j)*sum2;
            if (j~=1)
                nonzeros=nonzeros+1;
                J_v(nonzeros)=-sum2*V(j)*V(i); %dWq_i/dtheta_j
            end
        end
        if (i~=1)
            nonzeros=nonzeros+1;
            J_v(nonzeros)=V(i)*sum3; %dWq_i/dtheta_i
        end
    end
    dQinj_dtheta = sparse(row_indexes.dQinj_dtheta,column_indexes.dQinj_dtheta, J_v, N_meas_Qinj, N_bus); %derivatives of Q_inj w.r.t theta
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQinj/dV
%**************************************************************************
J_v=zeros(number_nonzeros.dQinj_dV,1);
nonzeros=0;
for i0=1:N_meas_Qinj
    i=ind_meas_Qinj(i0);
    sum1=0;
    N1=number_connected_nodes(i);
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        sum0=G_values(index_y)*sin(theta(i)-theta(j))-B_values(index_y)*cos(theta(i)-theta(j));
        sum1=sum1+V(j)*sum0;
        nonzeros=nonzeros+1;
        J_v(nonzeros)=V(i)*sum0; %dWq_i/dV_j
    end
    index_y=index_y+1;
    nonzeros=nonzeros+1;
    J_v(nonzeros)=-2*B_values(index_y)*V(i)+sum1; %dWq_i/dV_i
end
dQinj_dV = sparse(row_indexes.dQinj_dV,column_indexes.dQinj_dV, J_v, N_meas_Qinj, N); %derivatives of Q_inj w.r.t V


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPij/dtheta
%**************************************************************************
J_v=zeros(number_nonzeros.dPij_dtheta,1);
nonzeros=0;
for i0=1:N_meas_Pij
    l=ind_meas_Pij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        if (i~=1)
            nonzeros=nonzeros+1;%dPij/dtheta_i
            J_v(nonzeros)=-V(i) * V(j)*(-branch_g(l)*sin(theta(i)-theta(j))+branch_b(l)*cos(theta(i)-theta(j)));
        end
        if (j~=1)
            nonzeros=nonzeros+1;%dPij/dtheta_j
            J_v(nonzeros)=-V(i) * V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)));
        end
    else
        if (i~=1)
            nonzeros=nonzeros+1;%dPij/dtheta_i
            J_v(nonzeros)=-V(i) * V(j)*(-branch_g(l)*sin(theta(i)-theta(j))+branch_b(l)*cos(theta(i)-theta(j)))/branch_Kt_real(l);
        end
        if (j~=1)
            nonzeros=nonzeros+1;%dPij/dtheta_j
            J_v(nonzeros)=-V(i) * V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)))/branch_Kt_real(l);
        end
    end
end
dPij_dtheta = sparse(row_indexes.dPij_dtheta,column_indexes.dPij_dtheta, J_v, N_meas_Pij, N_bus); %derivatives of Pij w.r.t theta


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPij/dV
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dPij_dV,1);
    nonzeros=0;
    for i0=1:N_meas_Pij
        l=ind_meas_Pij(i0);
        i=branch_index_of_1_node(l);
        j=branch_index_of_2_node(l);
        if (IsTransformer(l)==0)
            nonzeros=nonzeros+1;%dPij/dV_i
            J_v(nonzeros)=2*V(i)*branch_g(l)-V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)));
            nonzeros=nonzeros+1;%dPij/dV_j
            J_v(nonzeros)=-V(i) *(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)));
        else
            nonzeros=nonzeros+1;%dPij/dV_i
            J_v(nonzeros)=2*V(i)*branch_g(l)/(branch_Kt_real(l)^2)- V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)))/branch_Kt_real(l);
            nonzeros=nonzeros+1;%dPij/dV_j
            J_v(nonzeros)=-V(i) *(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)))/branch_Kt_real(l);
        end
    end
    dPij_dV = sparse(row_indexes.dPij_dV,column_indexes.dPij_dV, J_v, N_meas_Pij, N); %derivatives of Pij w.r.t V
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPji/dtheta
%**************************************************************************
J_v=zeros(number_nonzeros.dPji_dtheta,1);
nonzeros=0;
for i0=1:N_meas_Pji
    l=ind_meas_Pji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        if (i~=1)
            nonzeros=nonzeros+1;%dPji/dtheta_i
            J_v(nonzeros)=-V(i)*V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)));
        end
        if (j~=1)
            nonzeros=nonzeros+1;%dPji/dtheta_j
            J_v(nonzeros)=-V(i)*V(j)*(-branch_g(l)*sin(theta(j)-theta(i))+branch_b(l)*cos(theta(j)-theta(i)));
        end
    else
        if (i~=1)
            nonzeros=nonzeros+1;%dPji/dtheta_i
            J_v(nonzeros)=-V(i)*V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)))/branch_Kt_real(l);
        end
        if (j~=1)
            nonzeros=nonzeros+1;%dPji/dtheta_j
            J_v(nonzeros)=-V(i)*V(j)*(-branch_g(l)*sin(theta(j)-theta(i))+branch_b(l)*cos(theta(j)-theta(i)))/branch_Kt_real(l);
        end
    end
end
dPji_dtheta = sparse(row_indexes.dPji_dtheta,column_indexes.dPji_dtheta, J_v, N_meas_Pji, N_bus); %derivatives of Pji w.r.t theta


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dPji/dV
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dPji_dV,1);
    nonzeros=0;
    for i0=1:N_meas_Pji
        l=ind_meas_Pji(i0);
        i=branch_index_of_1_node(l);
        j=branch_index_of_2_node(l);
        if (IsTransformer(l)==0)
            nonzeros=nonzeros+1;%dPji/dV_i
            J_v(nonzeros)=-V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)));
            nonzeros=nonzeros+1;%dPji/dV_j
            J_v(nonzeros)=2*V(j)*branch_g(l)-V(i)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)));
        else
            nonzeros=nonzeros+1;%dPji/dV_i
            J_v(nonzeros)=-V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)))/branch_Kt_real(l);
            nonzeros=nonzeros+1;%dPji/dV_j
            J_v(nonzeros)=2*V(j)*branch_g(l)-V(i)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)))/branch_Kt_real(l);
        end
    end
    dPji_dV = sparse(row_indexes.dPji_dV,column_indexes.dPji_dV, J_v, N_meas_Pji, N); %derivatives of Pji w.r.t V
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQij/dtheta
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dQij_dtheta,1);
    nonzeros=0;
    for i0=1:N_meas_Qij
        l=ind_meas_Qij(i0);
        i=branch_index_of_1_node(l);
        j=branch_index_of_2_node(l);
        if (IsTransformer(l)==0)
            if (i~=1)
                nonzeros=nonzeros+1;%dQij/dtheta_i
                J_v(nonzeros)=-V(i) * V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)));
            end
            if (j~=1)
                nonzeros=nonzeros+1;%dQij/dtheta_j
                J_v(nonzeros)=-V(i) * V(j)*(-branch_g(l)*cos(theta(i)-theta(j))-branch_b(l)*sin(theta(i)-theta(j)));
            end
        else
            if (i~=1)
                nonzeros=nonzeros+1;%dQij/dtheta_i
                J_v(nonzeros)=-V(i) * V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)))/branch_Kt_real(l);
            end
            if (j~=1)
                nonzeros=nonzeros+1;%dQij/dtheta_j
                J_v(nonzeros)=-V(i) * V(j)*(-branch_g(l)*cos(theta(i)-theta(j))-branch_b(l)*sin(theta(i)-theta(j)))/branch_Kt_real(l);
            end
        end
    end
    dQij_dtheta = sparse(row_indexes.dQij_dtheta,column_indexes.dQij_dtheta, J_v, N_meas_Qij, N_bus); %derivatives of Qij w.r.t theta
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQij/dV
%**************************************************************************
J_v=zeros(number_nonzeros.dQij_dV,1);
nonzeros=0;
for i0=1:N_meas_Qij
    l=ind_meas_Qij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        nonzeros=nonzeros+1;%dQij/dV_i
        J_v(nonzeros)=-2*V(i)*(branch_b(l)-branch_b_sh(l)/2)-V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)));
        nonzeros=nonzeros+1;%dQij/dV_j
        J_v(nonzeros)=-V(i)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)));
    else
        nonzeros=nonzeros+1;%dQij/dV_i
        J_v(nonzeros)=-2*V(i)*(branch_b(l)-branch_b_sh(l)/2)/(branch_Kt_real(l)^2)-V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)))/branch_Kt_real(l);
        nonzeros=nonzeros+1;%dQij/dV_j
        J_v(nonzeros)=-V(i)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)))/branch_Kt_real(l);
    end
end
dQij_dV = sparse(row_indexes.dQij_dV,column_indexes.dQij_dV, J_v, N_meas_Qij, N); %derivatives of Qij w.r.t V


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQji/dtheta
%**************************************************************************
if (H_decoupled==0) %only for full Jacobian
    J_v=zeros(number_nonzeros.dQji_dtheta,1);
    nonzeros=0;
    for i0=1:N_meas_Qji
        l=ind_meas_Qji(i0);
        i=branch_index_of_1_node(l);
        j=branch_index_of_2_node(l);
        if (IsTransformer(l)==0)
            if (i~=1)
                nonzeros=nonzeros+1;%dQji/dtheta_i
                J_v(nonzeros)=-V(i)*V(j)*(-branch_g(l)*cos(theta(j)-theta(i))-branch_b(l)*sin(theta(j)-theta(i)));
            end
            if (j~=1)
                nonzeros=nonzeros+1;%dQji/dtheta_j
                J_v(nonzeros)=-V(i)*V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)));
            end
        else
            if (i~=1)
                nonzeros=nonzeros+1;%dQji/dtheta_i
                J_v(nonzeros)=-V(i)*V(j)*(-branch_g(l)*cos(theta(j)-theta(i))-branch_b(l)*sin(theta(j)-theta(i)))/branch_Kt_real(l);
            end
            if (j~=1)
                nonzeros=nonzeros+1;%dQji/dtheta_j
                J_v(nonzeros)=-V(i)*V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)))/branch_Kt_real(l);
            end
        end
    end
    dQji_dtheta = sparse(row_indexes.dQji_dtheta,column_indexes.dQji_dtheta, J_v, N_meas_Qji, N_bus); %derivatives of Qji w.r.t theta
end


%% ************************************************************************
%COMPUTING MATRIX OF DERIVATIVES dQji/dV
%**************************************************************************
J_v=zeros(number_nonzeros.dQji_dV,1);
nonzeros=0;
for i0=1:N_meas_Qji
    l=ind_meas_Qji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        nonzeros=nonzeros+1;%dQji/dV_i
        J_v(nonzeros)=-V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)));
        nonzeros=nonzeros+1;%dQji/dV_j
        J_v(nonzeros)=-2*V(j)*(branch_b(l)-branch_b_sh(l)/2)-V(i)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)));
    else
        nonzeros=nonzeros+1;%dQji/dV_i
        J_v(nonzeros)=-V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)))/branch_Kt_real(l);
        nonzeros=nonzeros+1;%dQji/dV_j
        J_v(nonzeros)=-2*V(j)*(branch_b(l)-branch_b_sh(l)/2)-V(i)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)))/branch_Kt_real(l);
    end
end
dQji_dV = sparse(row_indexes.dQji_dV,column_indexes.dQji_dV, J_v, N_meas_Qji, N); %derivatives of Qji w.r.t V


%% ************************************************************************
% CONSTRUCTING MEASUREMENT JACOBIAN
%**************************************************************************
if (H_decoupled==0)
    % full Jacobian
    H= ...
        [dV_dtheta     dV_dV;
        dPinj_dtheta   dPinj_dV;
        dQinj_dtheta   dQinj_dV;
        dPij_dtheta    dPij_dV;
        dPji_dtheta    dPji_dV;
        dQij_dtheta    dQij_dV;
        dQji_dtheta    dQji_dV];
else
    % decoupled Jacobian
    H= ...
        [dV_dtheta                  dV_dV;
        dPinj_dtheta                sparse(N_meas_Pinj, N);
        sparse(N_meas_Qinj, N_bus)  dQinj_dV;
        dPij_dtheta                 sparse(N_meas_Pij, N);
        dPji_dtheta                 sparse(N_meas_Pji, N);
        sparse(N_meas_Qij, N_bus)   dQij_dV;
        sparse(N_meas_Qji, N_bus)   dQji_dV];
end

if (H_sparse~=1)  %converting to dense form if required
    H=full(H);
end

end

 
 
 
 
 