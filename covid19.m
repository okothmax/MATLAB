function [] = corona()
clear all;
clc;

omega = 1/10;
delta = 30/100;
epsilonR = 8/10;

gammaA = 1/21;
gammaM = 1/14;
gammaH = 1/14;

criticalFrac = 0.67/100;
severeFrac = 2.6/100;

zeta = criticalFrac;
phi = 1/10;
kappa = severeFrac*omega
ICUDeathTime = 7;
lambdaC = (1/ICUDeathTime)*criticalFrac

RoA = 1.5;
RoM = 1.5;
beta1 = 1 * RoA * gammaA/delta;
beta2 = 2 * RoM * (gamma + kappa)/(1 - delta);
Ro = RoA + RoM;

function du = covid19(~,u)
    du = zeros(8,1);
    du(1) = beta1 * u(1) * u(3)/N - beta2 * u(1) * u(4)/N + epsilonR;
    du(2) = beta1 * u(1) * u(3)/N + beta2 * u(1) * u(4)/N - omega * u(2);
    du(3) = delta * omega * u(2) - gammaA * u(3);
    du(4) = (1 - delta) * omega * u(2) - (kappa + gammaM) * u(4);
    du(5) = kappa * u(4) + phi * u(6) - (zeta + gammaH) * u(5);
    du(6) = zeta * u(5) - (phi + lambdaC) * u(6);
    du(7) = gammaA * u(3) + gammaM * u(4) + gammaH * u(5) epsilonR;
    du(8) = lambdaC;
end
N = 76 * 10 ^ 7; 
Eo = 0.1; 
Ao = 0; 
Mo = 0; 
Ho = 0; 
Co = 0; 
Ro = 0; 
Do = 0;
So = N - Eo - Ao - Mo - Ho - Co - Ro - Do; 
options = odeset('RelTol', 1e-4, 'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4]);
[T1,U1] = ode45(@covid19, 0:1:365, [So Eo Ao Mo Ho Co Ro Do], options);
Date1 = datetime(2020,3,13) + caldays(1:length(T1));

figure(1)
plot(Date1,U1(:,1), 'b', Date1, U1(:,2), 'm', Date1, U1(:,3), 'r', Date1, U1(:,4), 'k', Date1, U1(:,5), '--r', Date1, U1(:,6),'--k', Date1, U1(:,7), 'g', Date1, U1(:,8), 'c', 'linewidth',1.2)
xlabel('Date');
ylabel('Humans');
legend('Susceptible', 'Exposed', 'Asymptomatic', 'Mild symptomatic', 'Hospitalized', 'Critical', 'Recovered', 'Fatalities', 'location', 'best')
grid on
grid minor
end