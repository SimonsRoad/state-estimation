clear all;
%% import datas
buses=importdata('buses.txt');
branches=importdata('branches.txt');
shunts=importdata('shunts.txt');
active_flow=importdata('Conventional Active Flow.txt');
reactive_flow=importdata('Conventional Reactive Flow.txt');
Pinj=importdata('Conventional Pinj.txt');
Qinj=importdata('Conventional Qinj.txt');
Voltage=importdata('Conventional Voltage.txt');
load('GoodMeasurement_14_bus.mat');

%% for topo

%combine brances and conventional active flow together

%revise 'busNumbers'
a=size(buses,1);
for m=1:a
    topo.busNumbers(m)=m;
end

%revise 'FromBusRef' and 'ToBusRef'
%FromBusRef must smaller than ToBusRef
for m=1:size(branches,1)
    if branches(m,1)>branches(m,2)
        c=branches(m,2);
        branches(m,2)=branches(m,1);
        branches(m,1)=c;
    end
end    
branches=sortrows(branches);
branches(:,1)=branches(:,1)-1;%make bus number start from 1
branches(:,2)=branches(:,2)-1;%make bus number start from 1
topo.FromBusRef=branches(:,1);
topo.ToBusRef=branches(:,2);

%revise nBus
topo.nBus=a;

%revise nBranch
topo.nBranch=size(branches,1);

%revise S_baseMVA
topo.S_baseMVA=branches(1,7);

%revise lineR , lineX...
topo.lineR=branches(:,4);%lineR
topo.lineX=branches(:,5);%lineX
topo.lineShuntB=branches(:,6);%lineX
topo.lineTapRatio=zeros(size(branches,1),1);%lineTapRatio are 0
topo.linePhaseShift=zeros(size(branches,1),1);%linePhaseShhift are 0
%bus shunt G and B
topo.busShuntG=zeros(size(branches,1),1);
topo.busShuntB=zeros(size(branches,1),1);
shunts(:,1)=shunts(:,1)-1;
for m=1:size(shunts,1)
    topo.busShuntG(shunts(m,1))=shunts(m,3);
    topo.busShuntB(shunts(m,1))=shunts(m,4);
end

%revise topo.internal
topo.internal.branch_info_1=sparse(topo.FromBusRef);
topo.internal.branch_info_2=sparse(topo.ToBusRef);
d=1;
j=0;
c=0;
%%for row_start_indexes
topo.internal.row_start_indexes=[];
for m=1:size(buses,1)
    if any(branches(:,1)==m)
        j=j+1;
        topo.internal.row_start_indexes(j,1)=c+d;
        c=c+d;
        b=find(any(branches(:,1)==m,2));%finding the rows start from bus m
        d=size(b,1);
    end
end
%%for number_connected_nodes
topo.internal.number_connected_nodes=[];
for m=1:size(buses,1)
    e=size(find(any(branches(:,1)==m,2)),1);
    f=size(find(any(branches(:,2)==m,2)),1);
    topo.internal.number_connected_nodes(m,1)=e+f;
end   
%%for node_all_connections
topo.internal.node_all_connections=[];
e=1;
for m=1:size(buses,1)
    if any(branches(:,1)==m)
        b=find(any(branches(:,1)==m,2));%index of rows whose 'from_buses' are equal to m
        c=branches(b,2);%buses index of b
        d=size(b,1);%size of b
        %insert bus index into node_all_connections        
        topo.internal.node_all_connections(e:e+d-1,1)=c;
        topo.internal.node_all_connections(e+d,1)=m;
        e=e+d+1;
    end
end

%topo.internal.branch_g and topo.internal.branch_b
A=1./(topo.lineR+1i*topo.lineX);
topo.internal.branch_g=real(A);
topo.internal.branch_b=imag(A);

%other parameters in topo.internal
topo.internal.number_nonzeros_Y=size(topo.internal.node_all_connections,1);
topo.internal.IsTransformer=zeros(size(branches,1),1);
topo.internal.lineTapRatioImag=zeros(size(branches,1),1);
topo.internal.lineTapRatioReal=zeros(size(branches,1),1);

%% for meas

%revise 'active flow' and 'reactive flow'
%actice flow
for m=1:size(active_flow,1)
    if active_flow(m,1)>active_flow(m,2)
        c=active_flow(m,2);
        active_flow(m,2)=active_flow(m,1);
        active_flow(m,1)=c;
        active_flow(m,4:size(active_flow,2)-1)=-active_flow(m,4:size(active_flow,2)-1);%from bus to bus change, sign of flow also change
    end
end    
active_flow=sortrows(active_flow);
active_flow(:,1)=active_flow(:,1)-1;%make bus number start from 1
active_flow(:,2)=active_flow(:,2)-1;%make bus number start from 1
meas.LineP_FT=active_flow(:,4:size(active_flow,2)-1);
meas.LineP_TF=-active_flow(:,4:size(active_flow,2)-1);
meas.std_LineP_FT=repmat(active_flow(:,size(active_flow,2)),1,size(active_flow,2)-4);
meas.std_LineP_TF=repmat(active_flow(:,size(active_flow,2)),1,size(active_flow,2)-4);
%reactice flow
for m=1:size(reactive_flow,1)
    if reactive_flow(m,1)>reactive_flow(m,2)
        c=reactive_flow(m,2);
        reactive_flow(m,2)=reactive_flow(m,1);
        reactive_flow(m,1)=c;
        reactive_flow(m,4:size(reactive_flow,2)-1)=-reactive_flow(m,4:size(reactive_flow,2)-1);
    end
end    
reactive_flow=sortrows(reactive_flow);
reactive_flow(:,1)=reactive_flow(:,1)-1;%make bus number start from 1
reactive_flow(:,2)=reactive_flow(:,2)-1;%make bus number start from 1
meas.LineQ_FT=active_flow(:,4:size(reactive_flow,2)-1);
meas.LineQ_TF=-active_flow(:,4:size(reactive_flow,2)-1);
meas.std_LineQ_FT=repmat(active_flow(:,size(reactive_flow,2)),1,size(active_flow,2)-4);
meas.std_LineQ_TF=repmat(active_flow(:,size(reactive_flow,2)),1,size(active_flow,2)-4);

%%%revise Voltage
meas.BusV=Voltage(:,3:size(Voltage,2)-1);
meas.std_BusV=repmat(Voltage(:,size(Voltage,2)),1,size(active_flow,2)-4);

%%%revise Pinj and Qinj
Pinj(:,1)=Pinj(:,1)-1;
Qinj(:,1)=Qinj(:,1)-1;
%%%Pinj
for m=1:size(buses,1)
    if any(Pinj(:,1)==m)
        b=find(any(Pinj(:,1)==m,2));
        meas.BusInjP(m,1:size(Pinj,2)-2)=Pinj(b,2:size(Pinj,2)-1);
    else
        meas.BusInjP(m,1:size(Pinj,2)-2)=NaN;
    end
end    
meas.std_BusInjP=repmat(Pinj(:,size(Pinj,2)),1,size(active_flow,2)-4);
%%%Qinj
for m=1:size(buses,1)
    if any(Qinj(:,1)==m)
        b=find(any(Qinj(:,1)==m,2));
        meas.BusInjQ(m,1:size(Qinj,2)-2)=Pinj(b,2:size(Qinj,2)-1);
    else
        meas.BusInjQ(m,1:size(Qinj,2)-2)=NaN;
    end
end    
meas.std_BusInjQ=repmat(Qinj(:,size(Qinj,2)),1,size(active_flow,2)-4);

%% save data
save('Measurement.mat', 'meas','topo');
