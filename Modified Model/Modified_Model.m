VarN = [0.0 0.25 0.5  0.75 1];
figure;
for i = 1:length(VarN)
ref.nH = VarN(i);
[t,y] = ode45(@SEIR,[0,200],[100,20,10,0,1000,20,30],odeset('RelTol',1e-2,'AbsTol',1e-4),ref);
%% Plotting Equations
% Human Variables
plot(t,y(:,3),'LineWidth',5)
hold on;
xlabel('t')
ylabel('Number of Humans')
title('A graph of Number of Infected Humans against Time for varying ITN use rates(nH)')
grid ;
grid minor;

end
legend("0", "0.25","0.5","0.75","1",'Location','Best');
%% Extended Ross Model System of Equations
function sol = SEIR(t,Initial,ref)
Sh = Initial(1) ;
Eh = Initial(2);
Ih = Initial(3);
Rh = Initial(4);
Sm = Initial(5);
Em = Initial(6);
Im = Initial(7);
triH = 0.0004;
betaH = 0.2;
uH = 0.00006;
sigmaH = 0.001;
alphaH = 0.06;
r = 0.04;
b = 0.15;
vm = 0.3;
triM = 0.07;
betaM = 0.09;
uM = 0.067;
sigmaM = 0.01;
alphaM = 0.055;
omega = 0.0014;
pM = 0.5;

vh = 0.7;
dSh = triH - ((b*betaH*Sh*Im)/(1+(vh*Im))) - ((uH+ref.nH)* Sh)+ (omega*Rh);
dEh = ((b*betaH*Sh*Im)/(1+(vh*Im))) - ((alphaH+uH)*(Eh));
dIh = (alphaH*Eh)-((r+uH+sigmaH)*Ih);
dRh = (r*Ih)-((uH+omega)*(Rh))+(ref.nH*Sh);
dSm = triM - ((b*betaM*Sm*Ih)/(1+(vm*Ih))) - ((uM+pM)* Sm);
dEm = ((b*betaM*Sm*Ih)/(1+(vm*Ih))) - ((alphaM+uM+pM)*(Em));
dIm = (alphaM*Em)-((uM+sigmaM+pM)*Im);
sol = [dSh; dEh; dIh; dRh; dSm; dEm; dIm];
end