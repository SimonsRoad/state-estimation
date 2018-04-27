function [ topo, ind_meas, N_meas ] =f_meas_indices_and_H_sparsity_pattern_v2017( topo, meas)
%% ************************************************************************
%DESCRIPTION OF FUNCTION (VERSION 2017)
%**************************************************************************
%This function has two purposes. First, for each type of measurements,
%it calculates the number of available measurements and indices of these
%measurements in the list of all possible measurements of this type. E.g.
%if we have three buses in the system and are interested in active power 
%injections at buses, in total we can have three measurements (one for each
%bus). Suppose that we only have the measuring devices at buses 1 and 3. 
%Then the number of available measurements will be 2 and their indices in 
%the list of all possible measurements will be [1, 3].
%Second, this function calculates the sparsity pattern of matrix H. The 
%sparsity pattern of each submatrix of H is computed separately and stored
%as a field of a MATLAB structure. This is later used when matrix H is
%constructed at each step of GN method.

%INPUTS
%topo              - structure that contains information on topology of the 
%                    system and parameters of system elements
%meas              - structure that contains values and standard deviations
%                    of all measurements in the system. If a particular 
%                    measurement is unavailable, its value is equal to NaN.

%OUTPUTS
%row_indexes_H     - structure that for each submatrix of matrix H contains 
%                    information on what are the row indices of nonzero
%                    elements in this submatrix (stored as a field of
%                    structure topo.internal)
%column_indexes_H  - structure that for each submatrix of matrix H contains 
%                    information on what are the column indices of nonzero
%                    elements in this submatrix (stored as a field of
%                    structure topo.internal)
%number_nonzeros_H - structure that for each submatrix of matrix H contains 
%                    the number of nonzero elements in this submatrix
%                    (stored as a field of structure topo.internal)
%ind_meas          - structure that for each type of measurement contains 
%                    info on what are the indexes of available measurements
%                    in a list of all possible measurements for this type
%N_meas            - structure that contains the number of available 
%                    measurements for each type of measurements


%% ************************************************************************
%INITIALIZATION OF VARIABLES
%**************************************************************************
%extracting data from the structure to speed up subsequent calculations
row_start_indexes=topo.internal.row_start_indexes;
node_all_connections=topo.internal.node_all_connections;
number_connected_nodes=topo.internal.number_connected_nodes;
branch_index_of_1_node=topo.internal.branch_info_1;
branch_index_of_2_node=topo.internal.branch_info_2;
N_bus=topo.nBus-1;
%inititializing internal and output variables
buf=4*(N_bus+sum(number_connected_nodes)+1); %maximum possible number of nonzeros
row_indexes=zeros(buf,1);
column_indexes=zeros(buf,1);


%% ************************************************************************
%GETTING INDICES OF AVAILABLE MEASUREMENTS OF EACH TYPE AND THEIR NUMBER
%**************************************************************************
%Obtaining the indices of available measurements in the list of all possible
%measurements for each type of measurement (e.g. what are the indices of
%available measurements of Pinj in the list of all possible measurements of Pinj)
ind_meas_V=find(~isnan(meas.BusV));
ind_meas_Pinj=find(~isnan(meas.BusInjP));
ind_meas_Qinj=find(~isnan(meas.BusInjQ));
ind_meas_Pij=find(~isnan(meas.LineP_FT));
ind_meas_Pji=find(~isnan(meas.LineP_TF));
ind_meas_Qij=find(~isnan(meas.LineQ_FT));
ind_meas_Qji=find(~isnan(meas.LineQ_TF));

%Obtaining the number of available measurements of each type
N_meas_V=length(ind_meas_V);
N_meas_Pinj=length(ind_meas_Pinj);
N_meas_Qinj=length(ind_meas_Qinj);
N_meas_Pij=length(ind_meas_Pij);
N_meas_Pji=length(ind_meas_Pji);
N_meas_Qij=length(ind_meas_Qij);
N_meas_Qji=length(ind_meas_Qji);

%saving the obtained vectors of indices to structure ind_meas
ind_meas=struct('ind_meas_V',ind_meas_V,'ind_meas_Pinj',ind_meas_Pinj,...
    'ind_meas_Qinj',ind_meas_Qinj,'ind_meas_Pij',ind_meas_Pij,...
    'ind_meas_Pji',ind_meas_Pji,'ind_meas_Qij',ind_meas_Qij,...
    'ind_meas_Qji',ind_meas_Qji);
%saving the obtained numbers of available measurements to structure N_meas
N_meas=struct('N_meas_V',N_meas_V,'N_meas_Pinj',N_meas_Pinj,...
    'N_meas_Qinj',N_meas_Qinj,'N_meas_Pij',N_meas_Pij,'N_meas_Pji',N_meas_Pji,...
    'N_meas_Qij',N_meas_Qij,'N_meas_Qji',N_meas_Qji);


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPinj/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pinj
    i=ind_meas_Pinj(i0);%%looking for the bus index i of the i0th meas_Pinj
    N1=number_connected_nodes(i);%%N1 is the number of buses connected to bus i
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        if (j~=1)
            nonzeros=nonzeros+1; %dPinj_i/dtheta_j
            row_indexes(nonzeros)=i0;
            column_indexes(nonzeros)=j-1;
        end
    end
    if (i~=1)
        nonzeros=nonzeros+1;  %dPinj_i/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
end
row_indexes_H.dPinj_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dPinj_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dPinj_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPinj/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pinj
    i=ind_meas_Pinj(i0);
    N1=number_connected_nodes(i);
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        nonzeros=nonzeros+1;  %dPinj_i/dV_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j;
    end
    nonzeros=nonzeros+1; %dPinj_i/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
end
row_indexes_H.dPinj_dV=row_indexes(1:nonzeros);
column_indexes_H.dPinj_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dPinj_dV=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQinj/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qinj
    i=ind_meas_Qinj(i0);
    N1=number_connected_nodes(i);
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        if (j~=1)
            nonzeros=nonzeros+1; %dQinj_i/dtheta_j
            row_indexes(nonzeros)=i0;
            column_indexes(nonzeros)=j-1;
        end
    end
    if (i~=1)
        nonzeros=nonzeros+1; %dQinj_i/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
end
row_indexes_H.dQinj_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dQinj_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dQinj_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQinj/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qinj
    i=ind_meas_Qinj(i0);
    N1=number_connected_nodes(i);
    index_y=row_start_indexes(i)-1;
    for i1=1:N1
        index_y=index_y+1;
        j=node_all_connections(index_y);
        nonzeros=nonzeros+1; %dQinj_i/dV_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j;
    end
    nonzeros=nonzeros+1; %dQinj_i/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
end
row_indexes_H.dQinj_dV=row_indexes(1:nonzeros);
column_indexes_H.dQinj_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dQinj_dV=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPij/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pij
    l=ind_meas_Pij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (i~=1)
        nonzeros=nonzeros+1; %dPij/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
    if (j~=1)
        nonzeros=nonzeros+1; %dPij/dtheta_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j-1;
    end
end
row_indexes_H.dPij_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dPij_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dPij_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPij/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pij
    l=ind_meas_Pij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    nonzeros=nonzeros+1; %dPij/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
    nonzeros=nonzeros+1; %dPij/dV_j
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=j;
end
row_indexes_H.dPij_dV=row_indexes(1:nonzeros);
column_indexes_H.dPij_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dPij_dV=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPji/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pji
    l=ind_meas_Pji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (i~=1)
        nonzeros=nonzeros+1; %dPji/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
    if (j~=1)
        nonzeros=nonzeros+1; %dPji/dtheta_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j-1;
    end
end
row_indexes_H.dPji_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dPji_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dPji_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dPji/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Pji
    l=ind_meas_Pji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    nonzeros=nonzeros+1; %dPji/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
    nonzeros=nonzeros+1; %dPji/dV_j
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=j;
end
row_indexes_H.dPji_dV=row_indexes(1:nonzeros);
column_indexes_H.dPji_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dPji_dV=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQij/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qij
    l=ind_meas_Qij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (i~=1)
        nonzeros=nonzeros+1; %dQij/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
    if (j~=1)
        nonzeros=nonzeros+1; %dQij/dtheta_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j-1;
    end
end
row_indexes_H.dQij_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dQij_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dQij_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQij/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qij
    l=ind_meas_Qij(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    nonzeros=nonzeros+1; %dQij/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
    nonzeros=nonzeros+1; %dQij/dV_j
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=j;
end
row_indexes_H.dQij_dV=row_indexes(1:nonzeros);
column_indexes_H.dQij_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dQij_dV=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQji/dtheta
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qji
    l=ind_meas_Qji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    if (i~=1)
        nonzeros=nonzeros+1; %dQji/dtheta_i
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=i-1;
    end
    if (j~=1)
        nonzeros=nonzeros+1; %dQji/dtheta_j
        row_indexes(nonzeros)=i0;
        column_indexes(nonzeros)=j-1;
    end
end
row_indexes_H.dQji_dtheta=row_indexes(1:nonzeros);
column_indexes_H.dQji_dtheta=column_indexes(1:nonzeros);
number_nonzeros_H.dQji_dtheta=nonzeros;


%% ************************************************************************
%GETTING SPARSITY PATTERN OF MATRIX OF DERIVATIVES dQji/dV
%**************************************************************************
nonzeros=0;
for i0=1:N_meas_Qji
    l=ind_meas_Qji(i0);
    i=branch_index_of_1_node(l);
    j=branch_index_of_2_node(l);
    nonzeros=nonzeros+1; %dQji/dV_i
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=i;
    nonzeros=nonzeros+1; %dQji/dV_j
    row_indexes(nonzeros)=i0;
    column_indexes(nonzeros)=j;
end
row_indexes_H.dQji_dV=row_indexes(1:nonzeros);
column_indexes_H.dQji_dV=column_indexes(1:nonzeros);
number_nonzeros_H.dQji_dV=nonzeros;


%% ************************************************************************
%saving sparsity patterns of all submatrices to internal topology structure
topo.internal.row_indexes_H=row_indexes_H;
topo.internal.column_indexes_H=column_indexes_H;
topo.internal.number_nonzeros_H=number_nonzeros_H;
end