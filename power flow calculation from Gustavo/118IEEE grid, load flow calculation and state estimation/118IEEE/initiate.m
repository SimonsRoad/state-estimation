clear
%
global Bus Branch Gentor Motor Load LF GentorLF MotorLF Ym;
%
Bus=load('buses.txt');
Branch_=load('branches.txt'); 
Gentor_=load('generators.txt');
Load_=load('loads.txt');
Motor_=load('motors.txt');
Shunt=load('shunts.txt');
Sbase=100; % System Base in MVA


normalization(Branch_,Gentor_,Motor_,Load_,Sbase)
% %
loadflow(Bus,Branch,Gentor,Load,Motor,Shunt)%the measurements are printed in files:'Bus Voltage Measured.txt'... by 'reportLF' function 
Bus=LF;Gentor=GentorLF;Motor=MotorLF;% Gives the result of the load flow
% %
