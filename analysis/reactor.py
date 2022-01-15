import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
import pandas as pd
from scipy.optimize import fsolve


model_file = '/home/moon/autokmc/rmg2kmc/examples/co_oxidation/chem_annotated.cti'

# Define reactor conditions
# initial molar concentrations
# CO_0 = 1.0  # kmol / m^3
# O2_0 = 1.0  # kmol / m^3
# Ar_0 = 1.0  # kmol / m^3
REACTOR_VOLUME = 1.0  # m^3
REACTOR_TEMPERATURE = 900  # K
REACTOR_PRESSURE = 100000.0  # 1 bar = 100000 Pa
MAX_SIMULATION_TIME = 1.0  # seconds
CONCENTRATIONS = {
    'CO(4)': 0.1,
    'O2(2)': 0.1,
    'Ar': 0.8,
}

x_CO = CONCENTRATIONS['CO(4)']
x_O2 = CONCENTRATIONS['O2(2)']
x_Ar = CONCENTRATIONS['Ar']


print(f'Reading in mechanism file {model_file}')
gas = ct.Solution(model_file, "gas")
surf = ct.Interface(model_file, "surface1", [gas])

# initialize T and P
gas.TPX = REACTOR_TEMPERATURE, REACTOR_PRESSURE, CONCENTRATIONS
surf.TP = REACTOR_TEMPERATURE, REACTOR_PRESSURE

# get the molecular weights from the gas object
MW_O2 = gas.molecular_weights[gas.species_names.index('O2(2)')]
MW_CO = gas.molecular_weights[gas.species_names.index('CO(4)')]
MW_Ar = gas.molecular_weights[gas.species_names.index('Ar')]

# Catalyst settings
# TODO figure out where any of this comes from
catalyst_weight = 4.24e-3
cat_site_per_wt = 5 * 61.67 * 1e-6 * 1e3  # [mol/kg] 1e-6mol/micromole, 1000g/kg
site_density = (
    surf.site_density * 1000
)  # [mol/m^2]cantera uses kmol/m^2, convert to mol/m^2
cat_area = (catalyst_weight * cat_site_per_wt) / site_density  # [m^3]

# TODO make sure starting coverage matches corresponding KMC simulation. Here the surface starts clean
surf.coverages = "X(1):1.0"


gas_reactor = ct.IdealGasReactor(gas)
gas_reactor.volume = REACTOR_VOLUME
surface_reactor = ct.ReactorSurface(surf, gas_reactor, A=cat_area)


# set up mass flow controllers
inlet = ct.Reservoir(gas)
exhaust = ct.Reservoir(gas)

FC_temp = 900
volume_flow = 1.0  # TODO volume flow?
molar_flow = volume_flow * ct.one_atm / (8.3145 * FC_temp)  # [mol/s]
mass_flow = molar_flow * (
    x_CO * MW_CO + x_O2 * MW_O2 + x_Ar * MW_Ar
)  # [kg/s]
mfc = ct.MassFlowController(inlet, gas_reactor, mdot=mass_flow)

outlet_mfc = ct.PressureController(gas_reactor, exhaust, master=mfc, K=0.01)


# Create the reactor network
tank = ct.IdealGasReactor(gas, energy='off', volume=REACTOR_VOLUME)
reactor_network = ct.ReactorNet([tank])

# Now compile a list of all variables for which we will store data
column_names = [tank.component_name(item) for item in range(tank.n_vars)]
column_names = ['pressure'] + column_names
time_history = pd.DataFrame(columns=column_names)

# Run the simulation
t = 0
counter = 1
y_A = []
t_steps = []
while t < MAX_SIMULATION_TIME:
    t = reactor_network.step()

    # store every 10th value
    if True:  # counter % 10 == 0:
        state = np.hstack([tank.thermo.P, tank.mass, tank.volume, tank.T, tank.thermo.X])
        y_A.append(tank.thermo.X[4])
        t_steps.append(t)

        # Update the dataframe
        time_history.loc[t] = state
    counter += 1

N = -1

# plot concentration over time
plt.plot(time_history.index[0:N], time_history['CO2(3)'].values[0:N])
plt.xlabel('time')
plt.ylabel('concentration')
plt.show()
plt.clf()
