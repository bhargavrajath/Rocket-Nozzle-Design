import sympy as sp
import numpy as np
import math
import matplotlib.pyplot as plt

F = 20 # required thrust in kg
Pc = 10 # chamber pressure in bar/atm
Tc = 1500 # combustion temperature in kelvin
Dc = 0.05 # chamber dia
rho = 1800 # solid propellant density kg/m^3
k = 1.05 # ratio of specific heats for the exhaust gas
M = 42 # molecular mass of combustion products
Pe = 1 # ambient pressure in atm/bar
a = 0.0665 # solid propellant burn rate constant
n = 0.323 # solid propellant burn rate pressure exponent
R = 8314.3 # Universal gas constant
g = 9.806 # gravitational acceleration in m/s^2
con = 60 # nozzle convergence half angle in degrees
div = 15 # nozzle divergence half angle in degrees

ve = np.sqrt((2*k/(k-1))*(R*Tc/M)*(1-((Pe/Pc)**((k-1)/k)))) # exhaust exit velocity
mdot = F*g/ve # mass flow rate in kg/s (= fuel consumption rate at constant chamber pressure)
Me = np.sqrt((2/(k-1))*(((Pc/Pe)**((k-1)/k))-1)) # exit mach
At = (mdot/(Pc*10**5))*np.sqrt((R*Tc/(M*k))*((0.5*(k+1))**((k+1)/(k-1)))) # throat cross sectional area
Ae = (At/Me)*np.sqrt(((1 + 0.5*(k-1)*(Me**2))/(1 + 0.5*(k-1)))**((k+1)/(k-1))) # exit cross sectional area
Te = Tc/(1+(0.5*(k-1)*Me**2)) # Exhaust temperature

def tand(x):
    return sp.cos(x * sp.pi / 180)

Dt = 2*math.sqrt(At/sp.pi) # throat dia
De = 2*math.sqrt(Ae/sp.pi) # exit dia
r = a*(Pc*10**5)**n # burn rate perpendicular to burning surface
Ab = mdot/(rho*r) # burning surface area required to maintain mdot
Lc = 0.5*(Dc - Dt)/tand(con) # length of convergent section
Ld = 0.5*(De - Dt)/tand(div) # length of divergent section

print('Burn area required (m^2) =',Ab)
print('Mass flow rate (kg/s) =',mdot)
print('Exhaust velocity (m/s) =',ve)
print('Exhaust temperature (K) =',Te)
print('Throat dia (cm) =',Dt*100)
print('Exit dia (cm) =',De*100)
print('Convergence length (cm) =',Lc*100)
print('Divergence length (cm) =',Ld*100)

if(Dc >= De):
    shift = Dc
else:
    shift = De

yu = 0.5*np.array([shift+Dc,shift+Dc,shift+Dt,shift+De])
yl = 0.5*np.array([shift-Dc,shift-Dc,shift-Dt,shift-De])
yax = 0.5*shift*np.ones(4)
L = 0.1*(Lc+Ld)
x = np.array([0,L,L+Lc,L+Lc+Ld])

if(Dc > L+Lc+Ld):
    ax = Dc
else:
    ax = L+Lc+Ld

ax = ax + 0.1*ax

plt.figure(1)
plt.plot(x,yu,'r')
plt.plot(x,yl,'r')
plt.plot(x,yax,'k')
plt.axis([0,ax,0,ax])
plt.grid(True)
plt.xlabel('length (m)')
plt.ylabel('width (m)')
plt.title('Nozzle profile')

plt.figure(2)
gam = np.linspace(1,2,100)
vel = np.sqrt((2*gam/(gam-1))*(R*Tc/M)*(1-((Pe/Pc)**((gam-1)/gam))))
plt.subplot(311)
plt.plot(gam,vel)
plt.grid(True)
plt.ylabel('m/s')
plt.title('Exit velocity as a function of gamma')

To = np.linspace(273,3000,1000)
vel = np.sqrt((2*k/(k-1))*(R*To/M)*(1-((Pe/Pc)**((k-1)/k))))
plt.subplot(312)
plt.plot(To,vel)
plt.grid(True)
plt.ylabel('m/s')
plt.title('Exit velocity as a function of combustion temperature(K)')

Mdot = np.linspace(0,5,50)
vel = F*g/Mdot
Po = 1-(((vel**2)/(R*Tc))*0.5*(k-1)/k)
Po = Pe*Po**(k/(1-k))
plt.subplot(313)
plt.plot(Mdot,Po)
plt.grid(True)
plt.ylabel('bar')
plt.xlabel('kg/s')
plt.title('Chamber pressure required for various mass flow rates for achieving same thrust')

plt.show()

