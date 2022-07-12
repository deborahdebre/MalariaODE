%% Spring 2021 Differential Equations Final Project
%
%  Team Members :
%
% * Deborah Enya Esi Debre
% * Joel Adu-Kwarteng
%% Ross Model
[tR,yR] = ode45(@Ross,[0,500],[0.01,0.01],odeset());
InitR = [0.01,0.01];

%% Plotting Equations
figure;
hold on;
plot(tR,yR(:,1),'LineWidth',5)
plot(tR,yR(:,2),'LineWidth',5)
legend('Infected Humans ( I_{h} )','Infected Mosquitoes ( I_{m} )','Location','best');
xlabel('t')
ylabel('Number')
title('A graph of Number of specie against Time')
grid ;
grid minor;
hold off;

%% Ross Model System of Equations
function Rsol = Ross(tR,InitR)

IhR = InitR(1);
ImR = InitR(2);
a = 0.4;
m = 0.6;
b = 0.3;
c = 0.5;
r = 0.002;
uR = 0.2;

dImR = ((a*c*IhR)*(1-ImR))-(uR*ImR);
dIhR = ((a*b*m*ImR)*(1-IhR)) - (r*IhR);

Rsol = [dIhR; dImR];
end