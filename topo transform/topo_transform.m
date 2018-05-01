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
load('theta')

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
topo.lineShuntB=-branches(:,6);%lineB
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


%%for number_connected_nodes
topo.internal.number_connected_nodes=[];
for m=1:size(buses,1)
    e=size(find(any(branches(:,1)==m,2)),1);
    f=size(find(any(branches(:,2)==m,2)),1);
    topo.internal.number_connected_nodes(m,1)=e+f;
end   
%%for node_all_connections structure:x x x x 1 x x x 2 x x 3... x x N_buses
topo.internal.node_all_connections=[];

e=1;
for m=1:size(buses,1)
    b=find(any(branches(:,1)==m,2));%index of rows whose 'from_buses' are equal to bus m
    c=branches(b,2);%buses index of b
    g=find(any(branches(:,2)==m,2));%index of rows whose 'to_buses' are equal to bus m
    h=branches(g,1);%buses index of g
    k=sort([h;c]);
    d=size(k,1);%size of k
    %insert bus index into node_all_connections        
    topo.internal.node_all_connections(e:e+d-1,1)=k;
    topo.internal.node_all_connections(e+d,1)=m;
    e=e+d+1;
    
end

%%for row_start_indexes
topo.internal.row_start_indexes=[];
c=1;
d=0;
for m=1:size(buses,1)
    j=j+1;
    topo.internal.row_start_indexes(j,1)=c+d;
    c=c+d;
    b=find(any(branches(:,1)==m,2));%finding the rows start from bus m
    g=find(any(branches(:,2)==m,2));%finding rows whose 'from_buses' are equal to bus m
    l=[b;g];
    d=size(l,1)+1;
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

%%%revise Voltage
meas.BusV=Voltage(:,3:size(Voltage,2)-1);
meas.std_BusV=repmat(Voltage(:,size(Voltage,2)),1,size(active_flow,2)-4);
Voltage(:,1)=Voltage(:,1)-1;%make bus number start from 1

%revise 'active flow' and 'reactive flow'
%actice flow
active_flow(:,1)=active_flow(:,1)-1;%make bus number start from 1
active_flow(:,2)=active_flow(:,2)-1;%make bus number start from 1
for m=1:size(active_flow,1)
    if active_flow(m,1)>active_flow(m,2)%the bus index in frombus must smaller than tobus
        V1=Voltage(active_flow(m,1),3:size(Voltage,2)-1);
        V2=Voltage(active_flow(m,2),3:size(Voltage,2)-1);
        theta1=theta(active_flow(m,1),:)/180*pi;%from degree to rad
        theta2=theta(active_flow(m,2),:)/180*pi;%from degree to rad
        V1=V1.*cos(theta1)+1i*V1.*sin(theta1);%V1 and V2 in complex form
        V2=V2.*cos(theta2)+1i*V2.*sin(theta2);
        P1=active_flow(m,4:size(active_flow,2)-1);
        %obtain the R between bus V1 and V2
        g_index=find(branches(:,1)==active_flow(m,2) & branches(:,2)==active_flow(m,1));%because the bus:active_flow(m,1) and bus: active_flow(m,2) have to switch place, so branches(:,1)==active_flow(m,2) and branches(:,2)==active_flow(m,1)) 
        g=topo.internal.branch_g(g_index);
        active_flow(m,4:size(active_flow,2)-1)=-(P1-(abs(V1-V2)).^2*g);%Pij changes, Pij+Pji=g*|V1-V2|^2,power system analysis P22,equation (3.7)
        c=active_flow(m,2);
        active_flow(m,2)=active_flow(m,1);
        active_flow(m,1)=c;       
    end
end    
active_flow=sortrows(active_flow);
meas.LineP_FT=active_flow(:,4:size(active_flow,2)-1);
%calculate P_TF based on P_FT
meas.LineP_TF=zeros(size(meas.LineP_FT));
for m=1:size(active_flow,1)
	V1=Voltage(active_flow(m,1),3:size(Voltage,2)-1);
	V2=Voltage(active_flow(m,2),3:size(Voltage,2)-1);
	theta1=theta(active_flow(m,1),:)/180*pi;%from degree to rad
	theta2=theta(active_flow(m,2),:)/180*pi;%from degree to rad
	V1=V1.*cos(theta1)+1i*V1.*sin(theta1);%V1 and V2 in complex form
	V2=V2.*cos(theta2)+1i*V2.*sin(theta2);
    P1=active_flow(m,4:size(active_flow,2)-1);
    g_index=find(branches(:,1)==active_flow(m,1) & branches(:,2)==active_flow(m,2));
    g=topo.internal.branch_g(g_index);
    meas.LineP_TF(m,:)=-(P1-(abs(V1-V2)).^2*g);
end    
meas.std_LineP_FT=repmat(active_flow(:,size(active_flow,2)),1,size(active_flow,2)-4);
meas.std_LineP_TF=repmat(active_flow(:,size(active_flow,2)),1,size(active_flow,2)-4);

%reactice flow
reactive_flow(:,1)=reactive_flow(:,1)-1;%make bus number start from 1
reactive_flow(:,2)=reactive_flow(:,2)-1;%make bus number start from 1
for m=1:size(reactive_flow,1)
    if reactive_flow(m,1)>reactive_flow(m,2)
        U1=Voltage(reactive_flow(m,1),3:size(Voltage,2)-1);
        U2=Voltage(reactive_flow(m,2),3:size(Voltage,2)-1);
        theta1=theta(reactive_flow(m,1),:)/180*pi;%from degree to rad
        theta2=theta(reactive_flow(m,2),:)/180*pi;%from degree to rad
        V1=U1.*cos(theta1)+1i*U1.*sin(theta1);%V1 and V2 in complex form
        V2=U2.*cos(theta2)+1i*U2.*sin(theta2);
        Q1=reactive_flow(m,4:size(reactive_flow,2)-1);
        b_line_index=find(branches(:,1)==reactive_flow(m,2) & branches(:,2)==reactive_flow(m,1));
        b_line=topo.internal.branch_b(b_line_index);
        b_sh_index=find(branches(:,1)==reactive_flow(m,2) & branches(:,2)==reactive_flow(m,1));
        b_sh=-topo.lineShuntB(b_sh_index)/2;%should divided by two and add a minus sign ahead
        Q_loss=-b_sh*(U1.^2+U2.^2)-b_line*(abs(V1-V2).^2);
        reactive_flow(m,4:size(reactive_flow,2)-1)=-(Q1-Q_loss);       
        c=reactive_flow(m,2);
        reactive_flow(m,2)=reactive_flow(m,1);
        reactive_flow(m,1)=c;
    end
end    
reactive_flow=sortrows(reactive_flow);
meas.LineQ_FT=reactive_flow(:,4:size(reactive_flow,2)-1);
%calculate Q_TF based on Q_FT
meas.LineQ_TF=zeros(size(meas.LineP_FT));
for m=1:size(reactive_flow,1)
	U1=Voltage(reactive_flow(m,1),3:size(Voltage,2)-1);
	U2=Voltage(reactive_flow(m,2),3:size(Voltage,2)-1);
	theta1=theta(reactive_flow(m,1),:)/180*pi;%from degree to rad
	theta2=theta(reactive_flow(m,2),:)/180*pi;%from degree to rad
	V1=U1.*cos(theta1)+1i*U1.*sin(theta1);%V1 and V2 in complex form
	V2=U2.*cos(theta2)+1i*U2.*sin(theta2);
	Q1=reactive_flow(m,4:size(reactive_flow,2)-1);
	b_line_index=find(branches(:,1)==reactive_flow(m,1) & branches(:,2)==reactive_flow(m,2));
	b_line=topo.internal.branch_b(b_line_index);
	b_sh_index=find(branches(:,1)==reactive_flow(m,1) & branches(:,2)==reactive_flow(m,2));
	b_sh=-topo.busShuntB(b_sh_index)/2;%devided by 2??
	Q_loss=-b_sh*(U1.^2+U2.^2)-b_line*(abs(V1-V2).^2);
	meas.LineQ_TF(m,:)=-(Q1-Q_loss);
end 
meas.std_LineQ_FT=repmat(reactive_flow(:,size(reactive_flow,2)),1,size(reactive_flow,2)-4);
meas.std_LineQ_TF=repmat(reactive_flow(:,size(reactive_flow,2)),1,size(reactive_flow,2)-4);



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
meas.std_BusInjP=repmat(Pinj(1,size(Pinj,2)),a,size(active_flow,2)-4);
%%%Qinj
for m=1:size(buses,1)
    if any(Qinj(:,1)==m)
        b=find(any(Qinj(:,1)==m,2));
        meas.BusInjQ(m,1:size(Qinj,2)-2)=Qinj(b,2:size(Qinj,2)-1);
    else
        meas.BusInjQ(m,1:size(Qinj,2)-2)=NaN;
    end
end    
meas.std_BusInjQ=repmat(Qinj(1,size(Qinj,2)),a,size(active_flow,2)-4);

%% save data
save('Measurement.mat', 'meas','topo');
