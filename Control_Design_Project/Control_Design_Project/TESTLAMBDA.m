clc
clear all
close all

Psi=[];
lambdadot1=[];
%lambdadot2=[];
lambdax=0:0.0001:1;


for lambda=0:0.0001:1
   
    %Psi2=0.328*2180*burckhardt(lambda,2);
    Psi1=0.328*2180*burckhardt(lambda,1);
    lambdadot1=[lambdadot1,-((1-lambda)/(40))*(Psi1-600)];
    %lambdadot2=[lambdadot2,-((1-lambda)/(40))*(Psi2-600)];
    Psi=[Psi,Psi1];
end


L1=[lambdax;lambdadot1];
L2=[lambdax;zeros(1,10001)];
P=InterX(L1,L2);

figure('name','Stability Points')
plot(lambdax,Psi)
hold on
for i=1:size(P,2)
vline(P(1,i),'r')
hold on
end 
hline(600)

xlabel('$\lambda[-]$','Interpreter','latex')
ylabel('$\Psi\left(\lambda\right)[Nm]$','Interpreter','latex')

