
function mismatchlf(Bus,Ym,Gentor,Load,Motor)
global Pmist Qmist Pmistmax Qmistmax Pmot Qmot;
%Calculation of mismatches
%
%Size of vector of mismatches in load flow formulation
cp=0;cq=0;% auxiliar used to exclude swing bus eqs and Q eq in PV buses
for i=1:length(Bus(:,1))
    if Bus(i,2)==0 || Bus(i,2)==4   %This is for PQ buses
        cp=cp+1;
        cq=cq+1;
    end
    if Bus(i,2)==1% This is for PV buses
    cp=cp+1;
    end
end
%
Pmist=zeros(cp,1);Qmist=zeros(cq,1);Pi=zeros(cp,1);Qi=zeros(cq,1); % Mismatch vector definition
%
i1=1;i2=1;% auxiliar for mismatch vector
for i=1:length(Bus(:,1))
   if Bus(i,2)==0 || Bus(i,2)==4 % This is for PQ buses
      for j=1:length(Bus(:,1)) 
      Pi(i1)=Pi(i1)+Bus(i,3)*Bus(j,3)*(real(Ym(i,j))*cos(Bus(i,4)-Bus(j,4))+imag(Ym(i,j))*sin(Bus(i,4)-Bus(j,4)));
      Qi(i2)=Qi(i2)+Bus(i,3)*Bus(j,3)*(real(Ym(i,j))*sin(Bus(i,4)-Bus(j,4))-imag(Ym(i,j))*cos(Bus(i,4)-Bus(j,4)));
      end
      %
      if length(Motor)~=0   % Inclusion of motor powers
        for m=1:length(Motor(:,1))  
           if Motor(m,1)==Bus(i,1)
           Pi(i1)=Pi(i1)+Motor(m,5);  
           Qi(i2)=Qi(i2)+Motor(m,6);  
           end
        end
     end
      %     
      for l=1:length(Load(:,1))
          if Load(l,1)==Bus(i,1)
          Pi(i1)=Pi(i1)+Load(l,4);
          Qi(i2)=Qi(i2)+Load(l,5);
          end
      end
      %
   Pmist(i1)=-Pi(i1);  
   Qmist(i2)=-Qi(i2); 
   i1=i1+1;
   i2=i2+1;
  end
   %
  if Bus(i,2)==1% This is for PV buses
     for j=1:length(Bus(:,1))
     Pi(i1)=Pi(i1)+Bus(i,3)*Bus(j,3)*(real(Ym(i,j))*cos(Bus(i,4)-Bus(j,4))+imag(Ym(i,j))*sin(Bus(i,4)-Bus(j,4)));
     end
     %
     if length(Motor)~=0   % Inclusion of motor powers
        for m=1:length(Motor(:,1))  
           if Motor(m,1)==Bus(i,1)
           Pi(i1)=Pi(i1)+Motor(m,5);   
           end
        end
     end
      % 
      for g=1:length(Gentor(:,1)) % Inclusion of generator powers
          if Gentor(g,1)==Bus(i,1)
          Pi(i1)=Pi(i1)-Gentor(g,4);
          end
      end
      %
      for l=1:length(Load(:,1)) % Inclusion of load powers
          if Load(l,1)==Bus(i,1)
          Pi(i1)=Pi(i1)+Load(l,4);
          end
      end
      % 
  Pmist(i1)=-Pi(i1);   
  i1=i1+1;
  end
end
%
Pmistmax=max(abs(Pmist));% Max mismatch
Qmistmax=max(abs(Qmist));%Max mismatch


