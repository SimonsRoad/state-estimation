function reportLF(LF,Flows,Currents,GentorLF,Load,Motor,Shunt,Slack)
global Report

Report=zeros(length(LF(:,1)),10);
P=zeros(length(LF(:,1)),2);
Q=zeros(length(LF(:,1)),2);
B=zeros(length(LF(:,1)),2);
G=zeros(length(LF(:,1)),2);
%
for i=1:length(LF(:,1))
    LF(i,6)=LF(i,3)*LF(i,5);
    P(i,1)=LF(i,1);
    Q(i,1)=LF(i,1);
    B(i,1)=LF(i,1);
    G(i,1)=LF(i,1);
end
%
for i=1:length(P(:,1))
for g=1:length(GentorLF(:,1))
    if GentorLF(g,1)==P(i,1)
       P(i,2)=P(i,2)+GentorLF(g,4);
       Q(i,2)=Q(i,2)+GentorLF(g,5);
    end
end
%
for l=1:length(Load(:,1))
    if Load(l,1)==P(i,1)
       P(i,2)=P(i,2)-Load(l,4);
       Q(i,2)=Q(i,2)-Load(l,5);
    end
end
%
if length(Motor)~=0
    for m=1:length(Motor(:,1))
        if Motor(m,1)==P(i,1)
        P(i,2)=P(i,2)-Motor(m,5);
        Q(i,2)=Q(i,2)-Motor(m,6);
        end
   end
end
%
if Slack(1)==P(i,1) && Slack(4)==1
   P(i,2)=P(i,2)+Slack(2);
   Q(i,2)=Q(i,2)+Slack(3);
end
end
%
for i=1:length(G(:,1))
    if isempty(Shunt)==0
        for s=1:length(Shunt(:,1))
           if Shunt(s,1)==G(i,1)
           G(i,2)=Shunt(s,3);
           B(i,2)=Shunt(s,4);
           end
       end
    end
end

sigmau=0.002*abs(LF(:,3));
sigmaip=0.02*abs(P(:,2));
sigmaiq=0.02*abs(Q(:,2));
sigmatp=0.02*abs(Flows(:,4));
sigmatq=0.02*abs(Flows(:,5));


for i=1:length(sigmaip)
    if sigmaip(i)==0 
       sigmaip(i)=0.001;
    end
    if sigmaiq(i)==0
       sigmaiq(i)=0.001; 
    end
end

for i=1:length(sigmatp)
    if sigmatp(i)==0 
       sigmatp(i)=0.001;
    end
    if sigmatq(i)==0
       sigmatq(i)=0.001; 
    end
end

Report=[LF,P(:,2),Q(:,2),G(:,2),B(:,2)]
Bus1=[LF(:,1:3),sigmau];
Bus2=[LF(:,1:2),LF(:,3:4)];
Flowsp=[Flows(:,1:4),sigmatp];
Flowsq=[Flows(:,1:3),Flows(:,5),sigmatq];
Sinjp=[Report(:,1),Report(:,7),sigmaip];
Sinjq=[Report(:,1),Report(:,8),sigmaiq];
Curre=[Currents];


d1=fopen('Bus Voltage Measured.txt','w');
for k=1:length(Bus1(:,1))
   fprintf(d1,'%3i %3i %8.8f %8.8f \n',Bus1(k,1:4));
end   
fclose(d1);

d2=fopen('Active Flow Measured.txt','w');
for k=1:length(Flowsp(:,1))
   fprintf(d2,'%3i %3i %3i %8.8f %8.8f \n',Flowsp(k,1:5));
end   
fclose(d2);

d3=fopen('Reactive Flow Measured.txt','w');
for k=1:length(Flowsq(:,1))
   fprintf(d3,'%3i %3i %3i %8.8f %8.8f \n',Flowsq(k,1:5));
end   
fclose(d3);

d4=fopen('Pinj Measured.txt','w');
for k=1:length(Sinjp(:,1))
   fprintf(d4,'%3i %8.8f %8.8f \n',Sinjp(k,1:3));
end   
fclose(d4);

d5=fopen('Qinj Measured.txt','w');
for k=1:length(Sinjp(:,1))
   fprintf(d5,'%3i %8.8f %8.8f \n',Sinjq(k,1:3));
end   
fclose(d5);

d6=fopen('I phasor PMU.txt','w');
for k=1:length(Curre(:,1))
   fprintf(d6,'%3i %3i %3i %8.8f %8.8f \n',Curre(k,1:5));
end   
fclose(d6);

d7=fopen('V phasor PMU.txt','w');
for k=1:length(Bus2(:,1))
   fprintf(d7,'%3i %3i %8.8f %8.8f \n',Bus2(k,1:4));
end   
fclose(d7);

