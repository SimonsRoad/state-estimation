clear  i pi
global  zpseudo  IPMUUT

% Select Buses with PMU
BusPMU=[24;40;59;69;75;80;100;103;113;114]; % Location of PMUs

% Network data input
% ===========================
A=load('branches.txt');
Bus=load('buses.txt');
Shunt=load('shunts.txt');

[m,n]=size(A);
ci=A(1:m,1); 
cj=A(1:m,2);
circ=A(1:m,3);
rgr=A(1:m,4);
xgr=A(1:m,5);
bgro=A(1:m,6);
a=A(1:m,8);

zgr=rgr+xgr*i;  % Branch impedance

for k=1:m;
    ygr(k,1)=1/zgr(k);    % Branch admitance
    ggr(k,1)=real(ygr(k));    
    bgr(k,1)=imag(ygr(k));
end

nbus=length(Bus(:,1));% Number of network nodes
B=[ci,cj,circ,ggr,bgr,bgro,a];

% V phasor- PMU
%========================
F=load('V phasor PMU.txt');
nupmu=length(F(:,1));
up=F(:,3);
c=F(:,4);

% I Phasor PMU
%=========================
K=load('I phasor PMU.txt');
nI=length(K(:,1));
nIc=length(K(:,1));
im=K(:,4);
ic=K(:,5);

i1=1;
for i=1:length(F(:,1))
    for j=1:length(BusPMU)
        if BusPMU(j)==F(i,1)
           vpmu(i1,:)=F(i,:);
           i1=i1+1;
        end
    end
end

i1=1;
for i=1:length(K(:,1))
    for j=1:length(BusPMU)
         if BusPMU(j)==K(i,1)
           ipmu(i1,:)=K(i,:);
           i1=i1+1;
        end
    end
end

% original values
vpmuog_=vpmu;
vog=vpmuog_(:,3);
ipmuog_=ipmu;
imog=ipmuog_(:,4);
icog=ipmuog_(:,5);

% Include Measurements Noise and Standard Deviation in PMU Measurements

for i=1:length(vpmu(:,1))
    nvm(i)=abs(vpmu(i,3)*0.0002/sqrt(3));  % 0.0002    % sqrt(3) is used to obtain the standard uncertainties
    nvc(i)=abs(1*0.01*pi/180/sqrt(3));  
end

for i=1:length(ipmu(:,1))
    nim(i)=abs(ipmu(i,4)*0.0003/sqrt(3)); %  0.0003
    nic(i)=abs(1*0.01*pi/180/sqrt(3));    % 0.01*pi
end

for i=1:length(vpmu(:,1))
    vpmu(i,3)=vpmu(i,3)-nvm(i)*sqrt(3)+2*rand(1)*nvm(i)*sqrt(3); %+randn(1)*nvm(i);
    vpmu(i,4)=vpmu(i,4)-nvc(i)*sqrt(3)+2*rand(1)*nvc(i)*sqrt(3); %+randn(1)*nvc(i);
end

for i=1:length(ipmu(:,1))
    ipmu(i,4)=ipmu(i,4)-nim(i)*sqrt(3)+2*rand(1)*nim(i)*sqrt(3); %+randn(1)*nim(i);
    ipmu(i,5)=ipmu(i,5)-nic(i)*sqrt(3)+2*rand(1)*nic(i)*sqrt(3); %+randn(1)*nic(i);
end


VPMU=zeros(length(vpmu(:,1)),6);
VPMU(:,1:4)=vpmu;
VPMU(:,5)=nvm;
VPMU(:,6)=nvc;
IPMU=zeros(length(ipmu(:,1)),7);
IPMU(:,1:5)=ipmu;
IPMU(:,6:7)=[nim',nic'];

UTtransvolt(VPMU,IPMU,B,Bus)
% pseudovolt(VPMU,IPMU,B,Bus)

% original values
Sog=zeros(length(zpseudo(:,1)),1);
for i=1:length(zpseudo(:,1))
    for k=1:length(F(:,1))
        if zpseudo(i,1)==F(k,1)
           if zpseudo(i,3)==1
              Sog(i)=F(k,3);
           else
              Sog(i)=F(k,4);
           end
        end
    end
end

% Transform Polar form of currents into Rectangular form
UTtranspolar(ipmu,nim,nic)
IPMUr=IPMUUT;

% transpolar(ipmu,nim,nic)
% IPMUr=IPMUn;

% original values

for i=1:length(IPMU(:,1))
    IRog(i,4)=ipmuog_(i,4)*cos(ipmuog_(i,5));
    IRog(i,5)=ipmuog_(i,4)*sin(ipmuog_(i,5));
end
irog=IRog(:,4);
iiog=IRog(:,5);

 i1=1;
 for i=1:length(irog(:,1))
     Iog(i1,:)=irog(i);
     Iog(i1+1,:)=iiog(i);
     i1=i1+2;
 end


%% Conventional Measurements

pis=[3;4;8;12;13;16;15;19;20;22;24;25;27;31;33;35;36;42;44;...
    46;47;49;52;53;54;55;61;66;70;77;79;83;85;86;89;90;...
    92;96;97;98;99;102;104;105;110;111;112;116;117;118];% Selects buses with measurements of Pinj and Qinj
pt=[1;13;4;5;6;9;10;11;15;17;19;22;23;28;30;35;36;39;41;49;53;55;56;...
    59;60;64;68;72;73;77;80;81;83;88;91;97;101;290;105;112;113;116;...
    117;118;122;127;128;132;134;135;143;144;146;147;149;154;156;158;...
    159;160;162;167;170;171;173;174;175;177;178;181;182;184;185;186]; % Selects branches with measurements of Pflow and Qflow
vi=[4;10;12;18;25;27;36;40;59;73;76;82;86;92;107;111;112]; % Selects buses with measurements of voltage 


% Read Injected Power - clean measurements
% ==========================
Sinjp=load('Pinj Measured.txt');
Sinjq=load('Qinj Measured.txt');

% Read Power Flow - clean measurements
% =====================
Flowsp=load('Active Flow Measured.txt');
Flowsq=load('Reactive Flow Measured.txt');

% Bus voltages - clean measurements
% ========================
Bus1=load('Bus Voltage Measured.txt');

%original values
uog_=Bus1(:,3);
Piog_=Sinjp;
Qiog_=Sinjq;
Ptog_=Flowsp;
Qtog_=Flowsq;

% sigma definition
sigmaip=Sinjp(:,3);
sigmaiq=Sinjq(:,3);
sigmatp=Flowsp(:,5);
sigmatq=Flowsq(:,5);
sigmau=Bus1(:,4);

%%RANDOM ERROR CALCULATION
% erru=sigmau.*randn(length(sigmau),1);
% errip=sigmaip.*randn(length(sigmaip),1);
% erriq=sigmaiq.*randn(length(sigmaiq),1);
% errtp=sigmatp.*randn(length(sigmatp),1);
% errtq=sigmatq.*randn(length(sigmatq),1);

erru=-sigmau+2*sigmau.*rand(length(sigmau),1);
errip=-sigmaip+2*sigmaip.*rand(length(sigmaip),1);
erriq=-sigmaiq+2*sigmaiq.*rand(length(sigmaiq),1);
errtp=-sigmatp+2*sigmatp.*rand(length(sigmatp),1);
errtq=-sigmatq+2*sigmatq.*rand(length(sigmatq),1);
    

%%RANDOM ERROR INCLUSION
Bus1(:,3)=Bus1(:,3)+erru;
Flowsp(:,4)=Flowsp(:,4)+errtp;
Flowsq(:,4)=Flowsq(:,4)+errtq;
Sinjp(:,2)=Sinjp(:,2)+errip;
Sinjq(:,2)=Sinjq(:,2)+erriq;


for i=1:length(Sinjp(:,1))
    for k=1:length(pis) 
        if i==pis(k)
           Piout(k,:)=Sinjp(i,:);
           Qiout(k,:)=Sinjq(i,:);
           Piogaux(k,:)=Piog_(i,:);
           Qiogaux(k,:)=Qiog_(i,:);
        end
    end
end

for i=1:length(Flowsp(:,1))
    for k=1:length(pt) 
        if i==pt(k)
           Ptout(k,:)=Flowsp(i,:);
           Qtout(k,:)=Flowsq(i,:);
           Ptogaux(k,:)=Ptog_(i,:);
           Qtogaux(k,:)=Qtog_(i,:);
        end
   end
end

for i=1:length(Bus1(:,1))
    for k=1:length(vi) 
        if i==vi(k)
           VolOut(k,:)=Bus1(i,:);
           uogaux(k,:)=uog_(i,:);
        end
    end
end

Piog=Piogaux(:,2);
Qiog=Qiogaux(:,2);
Ptog=Ptogaux(:,4);
Qtog=Qtogaux(:,4);
uog=uogaux;

%%
d1=fopen('V phasor PMU Measured.txt','w');
for k=1:length(VPMU(:,1))
   fprintf(d1,'%3i %3i %8.6f %8.6f %8.8f %8.8f \n',VPMU(k,1:6));
end   
fclose(d1);

d2=fopen('I phasor PMU Measured.txt','w');
for k=1:length(IPMU(:,1))
   fprintf(d2,'%3i %3i %3i %8.6f %8.6f %8.8f %8.8f \n',IPMU(k,1:7));
end   
fclose(d2);

d3=fopen('Pseudo Voltage PMU Measured.txt','w');
for k=1:length(zpseudo(:,1))
   fprintf(d3,'%3i %3i %3i %8.8f %8.14f %8.14f \n',zpseudo(k,1:6));
end   
fclose(d3);

d4=fopen('Conventional Voltage.txt','w');
for k=1:length(VolOut(:,1))
   fprintf(d4,'%3i %3i %8.8f %8.8f \n',VolOut(k,1:4));
end   
fclose(d4);

d5=fopen('Conventional Active Flow.txt','w');
for k=1:length(Ptout(:,1))
   fprintf(d5,'%3i %3i %3i %8.8f %8.8f \n',Ptout(k,1:5));
end   
fclose(d5);


d6=fopen('Conventional Reactive Flow.txt','w');
for k=1:length(Qtout(:,1))
   fprintf(d6,'%3i %3i %3i %8.8f %8.8f \n',Qtout(k,1:5));
end   
fclose(d6);

d7=fopen('Conventional Pinj.txt','w');
for k=1:length(Piout(:,1))
   fprintf(d7,'%3i %8.8f %8.8f \n',Piout(k,1:3));
end   
fclose(d7);

d8=fopen('Conventional Qinj.txt','w');
for k=1:length(Qiout(:,1))
   fprintf(d8,'%3i %8.8f %8.8f \n',Qiout(k,1:3));
end   
fclose(d8);

d9=fopen('I phasor PMU Measured Rectangular.txt','w');
for k=1:length(IPMUr(:,1))
   fprintf(d9,'%3i %3i %3i %8.6f %8.14f %8.14f \n',IPMUr(k,1:6));
end   
fclose(d9);


disp('                                     MEASUREMENTS READY!!!');