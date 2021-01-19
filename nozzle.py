"""

   	Burn area
		|
		|
		V
    Fuel Mdot <-------.
	|				  |
	|			  Burn rate
	|                 ^
	V				  |  (We assume equilibrium i.e., neutral burn)
	Chamber pressure--'
			|
			|nozzle throat
			V
		Exhaust_mdot
			|
			|Nozzle exit
			V
		Exhaust velocity
			|
			|
			V
		  Thrust
		  
Author: Rajath Bhargav

"""

import numpy as np
import matplotlib.pyplot as plt

ar = 0.01 # Burn area in m^2
Tc = 1500.0 # combustion temperature in kelvin
Dc = 0.05 # chamber dia in m
rho = 1500.0 # solid propellant density kg/m^3
k = 1.05 # ratio of specific heats for the exhaust gas
M = 42.0 # molecular mass of combustion products
Pe = 1.0 # exit pressure in atm/bar
R = 8314.3 # Universal gas constant
g = 9.806 # gravitational acceleration in m/s^2
con = 60 # convergence half angle in degrees
div = 10 # divergence half angle in degrees

Pc = np.linspace(2,10,100) # range for combustion chamber pressure at equilibrium in bar

# Robert's law: burn rate r = a*(Pc)^n
a = 0.01 # solid propellant burn rate constant
n = 0.1 # solid propellant burn rate pressure exponent
r = a*(Pc*10**5)**n # burn rate perpendicular to burning surface

mdot = r*rho*ar # mass flow rate in kg/s at equilibrium

ve = np.sqrt((2*k/(k-1))*(R*Tc/M)*(1-((Pe/Pc)**((k-1)/k)))) # exhaust exit velocity
Me = np.sqrt((2/(k-1))*(((Pc/Pe)**((k-1)/k))-1)) # exit mach
At = (mdot/(Pc*10**5))*np.sqrt((R*Tc/(M*k))*((0.5*(k+1))**((k+1)/(k-1)))) # throat cross sectional area
Ae = (At/Me)*np.sqrt(((1 + 0.5*(k-1)*(Me**2))/(1 + 0.5*(k-1)))**((k+1)/(k-1))) # exit cross sectional area
Te = Tc/(1+(0.5*(k-1)*Me**2)) # Exhaust temperature

F = mdot*ve # Thrust generated

Dt = 2*np.sqrt(At/np.pi) # throat dia
De = 2*np.sqrt(Ae/np.pi) # exit dia
Lc = 0.5*(Dc - Dt)/np.tan(con*np.pi/180) # length of convergent section
Ld = 0.5*(De - Dt)/np.tan(div*np.pi/180) # length of divergent section

print('Mean Chamber pressure (atm/bar) =',np.average(Pc))
print('Average Mass flow rate (kg/s) =',np.average(mdot))
print('Mean Exhaust velocity (m/s) =',np.average(ve))
print('Mean Exhaust temperature (K) =',np.average(Te))
print('Mean Thrust generated (kg-force) =',np.average(F/g))
print('Mean Throat dia (cm) =',np.average(Dt*100))
print('Mean Exit dia (cm) =',np.average(De*100))

plt.subplot(221)
plt.plot(Pc,F/g,'r')
plt.grid(True)
plt.ylabel('kg-force (kg)')
plt.title('Thrust v/s chamber pressure')

plt.subplot(222)
plt.plot(Pc,np.divide(Ae,At),'g')
plt.grid(True)
plt.ylabel('Ae/At')
plt.title('Nozzle expansion ratio required v/s chamber pressure')

plt.subplot(223)
plt.plot(Pc,mdot,'b')
plt.grid(True)
plt.ylabel('mass flow (kg/s)')
plt.title('Fuel consumption v/s chamber pressure')

plt.subplot(224)
plt.plot(Pc,Te,'y')
plt.grid(True)
plt.ylabel('kelvin (K)')
plt.title('Exhaust temperature v/s chamber pressure')

plt.show()
