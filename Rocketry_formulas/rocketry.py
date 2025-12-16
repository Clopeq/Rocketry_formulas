# 3rd December 2025

from scipy.optimize import bisect
import scipy.constants as const
import numpy as np

def calculate_expansion_ratio_mach(k: float, Mach_number: float):
   """
   Calculate nozzle area ratio for given exhaust machn number
   
   :param k: Isentropic exponent
   :type k: float
   :param Mach_number: Exhaust mach number
   :type Mach_number: float
   """

   gamma1 = (k-1)/2
   gamma2 = (k+1)/(2*(k-1))
   expansion_ratio = 1/Mach_number * ((1 + gamma1*Mach_number**2)/(1 + gamma1))**gamma2
   return expansion_ratio

def calculate_expansion_ratio_pressure(k: float, pressure_chamber: float , pressure_exit: float = 101325):
   """
   Calculate nozzle area ratio for given pressure ratio
   
   :param k: Isentropic exponent
   :type k: float
   :param pressure_total: [Pa] Total (chamber) pressure
   :type pressure_total: float
   :param pressure_exit: [Pa] Exit pressure
   :type pressure_exit: float
   """

   gamma1 = ((k+1)/2)**(1/(k-1))
   gamma2 = (k+1)/(k-1)
   expansion_ratio = gamma1 * (pressure_exit/pressure_chamber)**(1/k) * np.sqrt(gamma2*(1-(pressure_exit/pressure_chamber)**((k-1)/k)))
   return expansion_ratio

def calculate_cstar(area_throat: float, pressure_chamber: float, flow_rate: float):
   """
   Calculate characteristic velocity using it's definition
   
   :param area_throat: [m2] Nozzle critical area
   :type area_throat: float
   :param pressure_chamber: [Pa] Chamber pressure
   :type pressure_chamber: float
   :param flow_rate: [kg/s] Total propellant flow rate
   :type flow_rate: float
   """

   cstar = area_throat*pressure_chamber/flow_rate
   return cstar


def calculate_cstar_ideal(k: float, Rs: float , temperature_combustion: float):
   """
   Calculate ideal (theoretical) characteristic velocity
   
   :param k: [-] isentropic exponent
   :type k: float
   :param Rs: [J/kg-K] Specific gas constant
   :type Rs: float
   :param temperature_combustion: [K] Combustion temperature (total temperature)
   :type temperature_combustion: float
   """
   gamma = k * (2/(k+1))**((k+1)/(k-1))
   cstar = np.sqrt(Rs*temperature_combustion/gamma)
   return cstar

def calculate_Cf(thrust: float, area_throat: float, pressure_chamber: float):
   """
   Calculate thrust coefficient using it's definition
   
   :param thrust: [N] Rocket thrust
   :type thrust: float
   :param area_throat: [m2] Nozzle critical area
   :type area_throat: float
   :param pressure_chamber: [Pa] Chamber pressure
   """

   Cf = thrust/area_throat/pressure_chamber
   return Cf

def calculate_Cf_ideal(k: float, pressure_chamber: float , pressure_exit: float = 101325):
   """
   Calculate ideal (theoretical) thrust coefficient for optimally expanded nozzle
   
   :param k: [-] isentropic exponent
   :type k: float
   :param pressure_chamber: [Pa] Chamber pressure (total pressure)
   :type pressure_chamber: float
   :param pressure_exit: [Pa] Nozzle exit pressure (ambient)
   :type pressure_exit: float
   """

   gamma1 = 2*k**2/(k-1) * (2/(k+1))**((k+1)/(k-1))
   gamma2 = (k-1)/k
   Cf = np.sqrt(gamma1 * (1-(pressure_exit/pressure_chamber)**gamma2))
   return Cf

def calculate_Cf_fixed_nozzle(k: float, area_exit: float, area_throat: float, pressure_chamber: float, pressure_exit: float = 101325, pressure_ambient: float = 101325):
   """
   Calculate ideal (theoretical) thrust coefficient for fixed expantion ratio nozzle.

   note: the nozzle does not have to be optimally expanded
   
   :param k: [-] isentropic exponent
   :type k: float
   :param area_exit: [m2] Nozzle exit cross section area
   :type area_exit: float
   :param area_throat: [m2] Nozzle throat cross section area
   :type area_throat: float
   :param pressure_chamber: [Pa] Chamber pressure
   :type pressure_chamber: float
   :param pressure_exit: [Pa] Exit pressure
   :type pressure_exit: float
   :param pressure_ambient: [Pa] Ambient pressure
   :type pressure_ambient: float
   """

   Cf = calculate_Cf_ideal(k, pressure_exit, pressure_chamber)
   pressure_coefficient = area_exit/area_throat * (pressure_exit-pressure_ambient)/pressure_chamber

   return Cf + pressure_coefficient

def calculate_impulse(thrust: list, time: list):
   """
   Calculate the total impulse from thrust curve.

   thrust and time have to be of the same length
   
   :param thrust: [N] Thrust values
   :type thrust: list
   :param time: [s] Time stamp values
   :type time: list
   """

   impulse = 0.0

   for i in len(1, thrust):
      impulse += 0.5*(thrust[i]+thrust[i-1]) * (time[i]-time[i-1])

   return impulse

def calculate_isp(thrust: float, flow_rate: float):
   """
   Calculate specific impulse from it's definition

   :param thrust: [N] Motor's thrust
   :type thrust: float
   :param flow_rate: [kg/s] Total proprellant flow rate
   :type flow_rate: float
   """

   isp = thrust/const.g/flow_rate
   return isp

def calculate_isp_vaq(velocity_effective: float):
   """
   Calculate specific impulse from effective exhaust velocity
   
   :param velocity_effective: [m/s] Effective exhaust velocity
   :type velocity_effective: float
   """

   isp = velocity_effective/const.g
   return isp

def calculate_isp_impulse(total_impulse: float, mass_consumed: float):
   """
   Calculate specific impulse from total delivered impulse and propellnt burned
   
   :param total_impulse: [Ns] Total impulse
   :type total_impulse: float
   :param mass_consumed: [kg] Propellant consumed
   :type mass_consumed: float
   """

   isp = total_impulse/const.g/mass_consumed
   return isp

def calculate_isp_cstar_Cf(cstar: float, Cf: float):
   """
   Calculate specific impulse from characteristic velocity and thrust coefficient
   
   :param cstar: [m/s] Characteristic velocity
   :type cstar: float
   :param Cf: [-] Thrust coefficient
   :type Cf: float
   """

   isp = cstar*Cf/const.g
   return isp

def calculate_isp_ideal(k: float, Rs: float, temperature_combustion: float, pressure_chamber: float, pressure_exit: float = 101325):
   """
   Calculate ideal (theoretical) specific impuls for optimally expanded nozzle.

   for vacuum performance pass 0 as exit pressure
   
   :param k: [-] isentropic exponent
   :type k: float
   :param Rs: [J/kg-K] Specific gas constant
   :type Rs: float
   :param temperature_combustion: [K] Combustion temperature
   :type temperature_combustion: float
   :param pressure_chamber: [Pa] Combusiton chamber pressure
   :type pressure_chamber: float
   :param pressure_exit: [Pa] Nozzle exit pressure
   :type pressure_exit: float
   """

   gamma1 = 2*k/(k-1)
   gamma2 = (k-1)/k

   isp = np.sqrt(Rs*temperature_combustion*gamma1 * (1-(pressure_exit/pressure_chamber)**gamma2))/const.g
   return isp

def calculate_pressure_exit(k: float, pressure_chamber: float, expansion_ratio: float):
   """
   Calculate static exit pressure
   
   :param k: [-] isentropic exponent
   :type k: float
   :param pressure_chamber: [Pa] Combusiton chamber pressure
   :type pressure_chamber: float
   :param expansion_ratio: [-] Nozzle area ratio
   :type expansion_ratio: float
   """
   
   f = lambda M: (calculate_expansion_ratio_mach(k, M) - expansion_ratio)
   mach = bisect(f, 1, 10)

   gamma1 = (k-1)/2
   gamma2 = k/(k-1)

   pressure_exit = pressure_chamber / (1+gamma1*mach**2)**gamma2
   return pressure_exit


def calculate_temperature_exit(k: float, temperature_combustion: float, expansion_ratio: float):
   """
   Calculate exhaust gas temperature at the nozzle exit
   
   :param k: [-] isentropic exponent
   :type k: float
   :param temperature_combustion: [K] Comsustion temperature
   :type temperature_combustion: float
   :param expansion_ratio: [-] Nozzle area ratio
   :type expansion_ratio: float
   """

   f = lambda M: (calculate_expansion_ratio_mach(k, M) - expansion_ratio)
   mach = bisect(f, 1, 10)

   gamma1 = (k-1)/2

   temperature_exit = temperature_combustion / (1+gamma1*mach**2)
   return temperature_exit

def calculate_density_exit(k: float, density_chamber: float, expansion_ratio: float): 
   """
   Calculate exhaust gas density at the nozzle exit
   
   :param k: [-] isentropic exponent
   :type k: float
   :param density_chamber: [kg/m3] Exhaust gas density inside combustion chamber (total density)
   :type density_chamber: float
   :param expansion_ratio: [-] Nozzle area ratio
   :type expansion_ratio: float
   """

   f = lambda M: (calculate_expansion_ratio_mach(k, M) - expansion_ratio)
   mach = bisect(f, 1, 10)

   gamma1 = (k-1)/2
   gamma2 = 1/(k-1)

   density_exit = density_chamber / (1+gamma1*mach**2)**gamma2
   return density_exit

def calculate_pressure_throat(k: float, pressure_chamber: float):
   """
   Calculate pressure at nozzle throat
   
   :param k: [-] isentropic exponent
   :type k: float
   :param pressure_chamber: [Pa] Combustion chamber pressure
   :type pressure_chamber: float
   """

   gamma = (2/(k+1))**(k/(k-1))
   pressure_throat = pressure_chamber*gamma
   return pressure_throat

def calculate_tempreature_throat(k: float, temperature_combustion: float):
   """
   Calculate exhaust gas temperature at the nozzle throat
   
   :param k: [-] isentropic exponent
   :type k: float
   :param temperature_combustion: [K] Comsustion temperature
   :type temperature_combustion: float
   """

   gamma = (2/(k+1))
   tempreature_throat = temperature_combustion*gamma
   return tempreature_throat

def calculate_density_throat(k: float, density_chamber: float):
   """
   Calculate exhaust gas density at the nozzle throat
   
   :param k: [-] isentropic exponent
   :type k: float
   :param density_chamber: kg/m3] Exhaust gas density inside combustion chamber (total density)
   :type density_chamber: float
   """
   gamma = (2/(k+1))**(1/(k-1))
   density_throat = density_chamber*gamma
   return density_throat

def calculate_thrust_ideal(k: float, area_throat: float, area_exit: float, pressure_chamber: float, pressure_exit: float=101325, pressure_ambient: float=101325):
   """
   Calculate rocket motor thrust

   
   :param k: [-] isentropic exponent
   :type k: float
   :param area_throat: [m2] Nozzle throat area
   :type area_throat: float
   :param area_exit: [m2] Nozzle exit area
   :type area_exit: float
   :param pressure_chamber: [Pa] Combustion chamber pressure
   :type pressure_chamber: float
   :param pressure_exit: [Pa] Nozzle exit pressure
   :type pressure_exit: float
   :param pressure_ambient: [Pa] Ambient pressure
   :type pressure_ambient: float
   """

   gamma1 = 2*k**2/(k-1) * (2/(k+1))**((k+1)/(k-1))
   gamma2 = (k-1)/k
   pressure_ratio = pressure_exit/pressure_chamber
   thrust = area_throat * pressure_chamber * np.sqrt(gamma1 * (1-pressure_ratio**gamma2)) + area_exit*(pressure_exit-pressure_ambient)
   return thrust

def calculate_velocity_exhaust(k: float, Rs: float, temperature_combustion: float, pressure_chamber: float, pressure_exit: float = 101325):
   """
   Calculate exhaust gas velocity
   
   :param k: [-] isentropic exponent
   :type k: float
   :param Rs: [J/kg-K] Specific gas constant
   :type Rs: float
   :param temperature_combustion: [K] Combustion temperature
   :type temperature_combustion: float
   :param pressure_chamber: [Pa] Combustion chamber pressure
   :type pressure_chamber: float
   :param pressure_exit: Pa] Nozzle exit pressure
   :type pressure_exit: float
   """

   gamma1 = 2*k/(k-1)
   gamma2 = (k+1)/k
   velocity_exhaust = np.sqrt(Rs*temperature_combustion*gamma1*(1-(pressure_exit/pressure_chamber)**gamma2))
   return velocity_exhaust

def calculate_flow_rate(k: float, Rs: float, area_throat: float, temperature_combustion:float, pressure_chamber: float):
   """
   Calculate propellant mass flow throught the nozzle
   
   :param k: [-] isentropic exponent
   :type k: float
   :param Rs: [J/kg-K] Specific gas constant
   :type Rs: float
   :param area_throat: [m2] Nozzle throat area
   :type area_throat: float
   :param temperature_combustion: [K] Combustion temperature
   :type temperature_combustion: float
   :param pressure_chamber: [K] Combustion chamber pressure
   :type pressure_chamber: float
   """

   gamma1 = (2/(k+1))**(1/2*(k+1)/(k-1))
   gamma2 = np.sqrt(k/Rs/temperature_combustion)
   flow_rate = pressure_chamber*area_throat*gamma1*gamma2
   return flow_rate
