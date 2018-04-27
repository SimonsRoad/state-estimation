function normalization(Branch_,Gentor_,Motor_,Load_,Sbase)
global  Branch Gentor Motor Load
%Change of bases for impedances (Vbase not included)
%BRANCHES
Branch=Branch_;
Gentor=Gentor_;
Motor=Motor_;
Load=Load_;
% 
% BRANCHES
for i=1:length(Branch(:,1))
    Branch(i,4)=Branch(i,4)*(Sbase/Branch(i,7));
    Branch(i,5)=Branch(i,5)*(Sbase/Branch(i,7));
    Branch(i,6)=Branch(i,6)/(Sbase/Branch(i,7));
end
%
%GENERATORS
for i=1:length(Gentor(:,1))
    Gentor(i,4)=Gentor(i,4)/Sbase;
end
%
%LOADS
for i=1:length(Load(:,1))
    Load(i,4)=Load(i,4)/Sbase;
    Load(i,5)=Load(i,5)/Sbase;
end
%
%MOTORS
if length(Motor)==0
    Motor=[];
else
   for i=1:length(Motor(:,1))
      Motor(i,4)=Motor(i,4)/Sbase;
      Motor(i,5)=Motor(i,5)/Sbase;
      Motor(i,7)=0.025; % initiated at value different from 0
      Motor(i,8)=Motor(i,8)*Sbase/Motor(i,3);
      Motor(i,9)=Motor(i,9)*Sbase/Motor(i,3);
      Motor(i,10)=Motor(i,10)*Sbase/Motor(i,3);
      Motor(i,11)=Motor(i,11)*Sbase/Motor(i,3);
      Motor(i,12)=Motor(i,12)*Sbase/Motor(i,3);
  end
end