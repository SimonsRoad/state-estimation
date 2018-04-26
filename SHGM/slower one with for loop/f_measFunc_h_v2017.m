function [ h ] = f_measFunc_h_v2017( V, theta, Y_bus, topo, ind_meas, N_meas)
%% ************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function computes the vector of operation parameters that correspond 
%to collected measurements (e.g. line flows, bus injections, voltage 
%magnitudes). Each subvector corresponding to a particular type of 
%measurement is computed in a separate loop. NOTE: formulas for line flows
%work only when there are no phase-shifting transformers in the system.

%INPUTS
%V        - vector of voltage magnitudes of ALL buses in the system
%theta    - vector of voltage phase angles for ALL buses in the system
%Y_bus    - structure that contains real and imaginary parts of admittances
%           of structurally nonzero elements in the admittance matrix
%internal - structure that contains information on topology of the 
%           system and parameters of system elements
%ind_meas - structure that for each type of measurement contains info
%           on what are the indexes of available measurements in a list of 
%           all possible measurements for this type
%N_meas   - structure that contains the number of available measurements
%           for each type of measurements 

%OUTPUTS
%h        - vector of operation parameters that correspond to given 
%           measurements (e.g. line flows, bus injections, etc)


%% ************************************************************************
%THEORY
%**************************************************************************
%node injections are computed according to the following formulas:
% Pinj = V .* ((G .* cos(theta) + B .* sin(theta)) * V);
% Qinj = V .* ((G .* sin(theta) - B .* cos(theta)) * V);

%line flows are computed according to the following formulas:
% Q_ij = - a_ij^2 * V(i)^2 * (b_ij + b_ij_sh/2) ...
%        + a_ij * V(i) * V(j) * b_ij * cos(theta(i) - theta(j)) ...
%        - a_ij * V(i) * V(j) * g_ij * sin(theta(i) - theta(j));
% Q_ji = - V(j)^2 * (b_ij + b_ij_sh/2) ...
%        + a_ij * V(i) * V(j) * b_ij * cos(theta(j) - theta(i)) ...
%        - a_ij * V(i) * V(j) * g_ij * sin(theta(j) - theta(i));
% P_ij =   a_ij^2 * V(i)^2 * g_ij ...
%        - a_ij * V(i) * V(j) * g_ij * cos(theta(i) - theta(j)) ...
%        - a_ij * V(i) * V(j) * b_ij * sin(theta(i) - theta(j));
% P_ji =   V(j)^2 * g_ij ...
%        - a_ij * V(i) * V(j) * g_ij * cos(theta(j) - theta(i)) ...
%        - a_ij * V(i) * V(j) * b_ij * sin(theta(j) - theta(i));
% where a_ij - 1/TapRatio_ij

   
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

%initializing the output vector
h=zeros(N_meas_V+N_meas_Pinj+N_meas_Qinj+N_meas_Pij+N_meas_Pji+N_meas_Qij+N_meas_Qji,1);


%% ************************************************************************
%COMPUTING VECTOR V
%**************************************************************************
h(1:N_meas_V)=V(ind_meas_V);


%% ************************************************************************
%COMPUTING VECTOR Pinj
%**************************************************************************
start_index=N_meas_V;
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
    end
    index_y=index_y+1;
    h(start_index+i0)=V(i)*(V(i)*G_values(index_y)+sum3);
end


%% ************************************************************************
%COMPUTING VECTOR Qinj
%**************************************************************************
start_index=start_index+N_meas_Pinj;
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
    end
    index_y=index_y+1;
    h(start_index+i0)=V(i)*(-V(i)*B_values(index_y)+sum1);
end


%% ************************************************************************
%COMPUTING VECTOR Pij
%**************************************************************************
start_index=start_index+N_meas_Qinj;
for i0=1:N_meas_Pij
    l=ind_meas_Pij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        h(start_index+i0)=V(i)^2*branch_g(l)-V(i) * V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)));
    else
        h(start_index+i0)=V(i)^2*branch_g(l)/(branch_Kt_real(l)^2)-V(i) * V(j)*(branch_g(l)*cos(theta(i)-theta(j))+branch_b(l)*sin(theta(i)-theta(j)))/branch_Kt_real(l);
    end
end


%% ************************************************************************
%COMPUTING VECTOR Pji
%**************************************************************************
start_index=start_index+N_meas_Pij;
for i0=1:N_meas_Pji
    l=ind_meas_Pji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        h(start_index+i0)=V(j)^2*branch_g(l)-V(i)*V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)));
    else
        h(start_index+i0)=V(j)^2*branch_g(l)-V(i)*V(j)*(branch_g(l)*cos(theta(j)-theta(i))+branch_b(l)*sin(theta(j)-theta(i)))/branch_Kt_real(l);
    end
end


%% ************************************************************************
%COMPUTING VECTOR Qij
%**************************************************************************
start_index=start_index+N_meas_Pji;
for i0=1:N_meas_Qij
    l=ind_meas_Qij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        h(start_index+i0)=-V(i)^2*(branch_b(l)-branch_b_sh(l)/2)-V(i) * V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)));
    else
        h(start_index+i0)=-V(i)^2*(branch_b(l)-branch_b_sh(l)/2)/(branch_Kt_real(l)^2)-V(i) * V(j)*(branch_g(l)*sin(theta(i)-theta(j))-branch_b(l)*cos(theta(i)-theta(j)))/branch_Kt_real(l);
    end
end


%% ************************************************************************
%COMPUTING VECTOR Qji
%**************************************************************************
start_index=start_index+N_meas_Qij;
for i0=1:N_meas_Qji
    l=ind_meas_Qji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (IsTransformer(l)==0)
        h(start_index+i0)=-V(j)^2*(branch_b(l)-branch_b_sh(l)/2)-V(i)*V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)));
    else
        h(start_index+i0)=-V(j)^2*(branch_b(l)-branch_b_sh(l)/2)-V(i)*V(j)*(branch_g(l)*sin(theta(j)-theta(i))-branch_b(l)*cos(theta(j)-theta(i)))/branch_Kt_real(l);
    end
end


end

 
 
 
 
 