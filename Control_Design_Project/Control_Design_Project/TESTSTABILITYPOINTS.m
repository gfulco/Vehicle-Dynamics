clear all
close all
clc

lambdadot=[];
lambdax=0:0.0001:1;


for lambda=0:0.0001:1
    Psi2=0.328*2180*burckhardt(lambda,1);
    lambdadot=[lambdadot,-((1-lambda)/(40))*(Psi2-600)];

end


L1=[lambdax;lambdadot];
L2=[lambdax;zeros(1,10001)];
[m,p]=min(lambdadot);


P=InterX(L1,L2);

figure('name','Stability Points')
plot(lambdax,lambdadot)
hold on
for i=1:(size(P,2))
    vline(P(1,i),'r')
    hold on
    plot(P(1,i),0,'r*')
    hold on
end
plot(lambdax,zeros(1,10001),'--')
hold on
vline(lambdax(p),'g')
xlabel('$\lambda$','Interpreter','latex')
ylabel('$\dot{\lambda}$','Interpreter','latex')
