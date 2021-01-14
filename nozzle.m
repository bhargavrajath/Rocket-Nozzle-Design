
F = 20; % required thrust in kg
Pc = 10; % chamber pressure in bar/atm
Tc = 1500; % combustion temperature in kelvin
Dc = 0.05; % chamber dia
rho = 1800; % solid propellant density kg/m^3
k = 1.05; % ratio of specific heats for the exhaust gas
M = 42; % molecular mass of combustion products
Pe = 1; % ambient pressure in atm/bar
a = 0.0665; % solid propellant burn rate constant
n = 0.323; % solid propellant burn rate pressure exponent
R = 8314.3; % Universal gas constant
g = 9.806; % gravitational acceleration in m/s^2
con = 60; % nozzle convergence half angle in degrees
div = 15; % nozzle divergence half angle in degrees

ve = sqrt((2*k/(k-1))*(R*Tc/M)*(1-((Pe/Pc)^((k-1)/k)))); % exhaust exit velocity
mdot = F*g/ve; % mass flow rate in kg/s (= fuel consumption rate at constant chamber pressure)
Me = sqrt((2/(k-1))*(((Pc/Pe)^((k-1)/k))-1)); % exit mach
At = (mdot/(Pc*10^5))*sqrt((R*Tc/(M*k))*((0.5*(k+1))^((k+1)/(k-1)))); % throat cross sectional area
Ae = (At/Me)*sqrt(((1 + 0.5*(k-1)*(Me^2))/(1 + 0.5*(k-1)))^((k+1)/(k-1))); % exit cross sectional area
Te = Tc/(1+(0.5*(k-1)*Me^2)); % Exhaust temperature

Dt = 2*sqrt(At/pi); % throat dia
De = 2*sqrt(Ae/pi); % exit dia
r = a*(Pc*10^5)^n; % burn rate perpendicular to burning surface
Ab = mdot/(rho*r); % burning surface area required to maintain mdot
Lc = 0.5*(Dc - Dt)/tand(con); % length of convergent section
Ld = 0.5*(De - Dt)/tand(div); % length of divergent section

disp('Burn area required (m^2) ='),disp(Ab);
disp('Mass flow rate (kg/s) ='),disp(mdot);
disp('Exhaust velocity (m/s) ='),disp(ve);
disp('Exhaust temperature (K) ='),disp(Te);
disp('Throat dia (cm) ='),disp(Dt*100);
disp('Exit dia (cm) ='),disp(De*100);
disp('Convergence length (cm) ='),disp(Lc*100);
disp('Divergence length (cm) ='),disp(Ld*100);

if(Dc >= De)
    shift = Dc;
else
    shift = De;
end

yu = 0.5*[shift+Dc,shift+Dc,shift+Dt,shift+De];
yl = 0.5*[shift-Dc,shift-Dc,shift-Dt,shift-De];
yax = 0.5*shift*ones(1,4);
L = 0.1*(Lc+Ld);
x = [0,L,L+Lc,L+Lc+Ld];

if(Dc > L+Lc+Ld)
    ax = Dc;
else
    ax = L+Lc+Ld;
end
ax = ax + 0.1*ax;
figure(1);
plot(x,yu,'r'),hold all;
plot(x,yl,'r'),plot(x,yax,'k'),axis([0,ax,0,ax]),grid on;
xlabel('length (m)'),ylabel('width (m)'),title('Rocket nozzle profile');

figure (2);
gam = 1:0.01:2;
vel = sqrt((2*gam./(gam-1)).*(R*Tc/M).*(1-((Pe/Pc).^((gam-1)./gam))));
subplot(311),plot(gam,vel),grid on,ylabel('m/s');
title('Exit velocity as a function of gamma');

To = 273:0.1:3000;
vel = sqrt((2*k/(k-1))*(R*To/M)*(1-((Pe/Pc)^((k-1)/k))));
subplot(312),plot(To,vel),grid on,ylabel('m/s');
title('Exit velocity as a function of combustion temperature(K)');

Mdot = 0:0.1:5;
vel = F*g./Mdot;
Po = 1-(((vel.^2)./(R*Tc))*0.5*(k-1)/k);
Po = Pe*Po.^(k/(1-k));
subplot(313),plot(Mdot,Po),grid on,ylabel('bar'),xlabel('kg/s');
title('Chamber pressure required for various mass flow rates for achieving same thrust');
