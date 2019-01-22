from biosysid import * 
# Import a model to create a biosysid Model object
model = import_sbml('sbml_example.xml')
# Supports SBML rate rules, initial assignments, 
# and other SBML level 3 core features (packages not supported)
# (Optional) model = import_bioscrape('bs.xml')
# (Optional) model = import_ode(ode_function)

# Print out parameters (dict)
# model.get_parameters()

# Print out species (dict)
# model.get_species()

# Create biosysid data objects
# properties_dict = {'Property1': 'Value1','Property2': 'Value2'}
# data = import_timeseries('data.csv',time_column = 1, value_column = 2, properties = properties_dict)
data = import_timeseries('test_data.csv', time_column = 2, value_column = 4, properties = {3 : 51})

# Plot experimental data, give a True value to the optional argument plot_show
# import_timeseries('test_data.csv', time_column = 2, value_column = 4, plot_show = True)

# Initial simulation data
import numpy as np
import matplotlib.pyplot as plt  
t_exp = data.get_keys()
timepoints = np.linspace(t_exp[0],t_exp[-1],len(t_exp))
bs_data, m = model.simulate_bioscrape(timepoints, 'deterministic')
simtime = timepoints
simdata = bs_data[:,m.get_species_index('c1')]
# simtime, simdata = model.simulate_roadrunner(timepoints)

# Create artificial data for model
lines = []
exp_value = []
for t,d in zip(simtime, simdata):
    v = d + np.random.normal(0,0.5)
    exp_value.append(v)
    lines.append([t,v])
import csv 
with open('test_data_noisy.csv','w') as f:
    writer = csv.writer(f)
    writer.writerows(lines)
    f.close()

plt.plot(simtime, exp_value)
plt.plot(simtime, simdata)
plt.title('Simulated model and artificial data')
plt.show()

data = import_timeseries('test_data_noisy.csv', time_column = 1, value_column = 2)
#Analysis
list_of_measurements = ['c1'] #outputs being measured
list_of_timepoints = [t_exp] # time range for which the output is measured
list_of_intervals = [0.1] # time interval between each measurements above
# initial_params = {'kc' : 0.6, 'k1' : 1} #optionally set initial guesses for parameters here
# initial_species = {'A' : 500} # optionally set initial conditions for species
'''
TODO: Implement analysis_sens
TODO: Implement analysis_ident
'''
#### model.analysis_sens(list_of_measurements, list_of_timepoints, list_of_intervals, initial_params, initial_species)
# Returns parameters that are identifiable
#### params = model.analysis_ident(list_of_measurements, list_of_timepoints, list_of_intervals, initial_params, initial_species)

# Parameter identification
priors = {'kc' : [1e-3, 1e3],'k1' : [1e-2, 1e5]}
model.parameters = {'kc':6, 'k1':1}
# model.species = ['c1']
params = model.parameters
fit_model, id_params = model.run_mcmc(params, priors, timepoints, nwalkers = 4* len(params), nsteps = 200, nsamples = 100, measurements = ['c1'], exp_data = data, plot_show = False)
# fit_model is biosysid model object with identified parameters substituted
res_orig, m_orig = model.simulate_bioscrape(timepoints, 'stochastic')
res, m = fit_model.simulate_bioscrape(timepoints, 'stochastic')

plt.figure()
simtime = timepoints
simdata_orig = res_orig[:,m_orig.get_species_index('c1')]
simdata = res[:,m.get_species_index('c1')]
plt.plot(simtime, exp_value, label ='exp')
plt.plot(simtime, simdata_orig, label = 'original model')
plt.plot(simtime, simdata, label = 'identified model')
plt.legend()
plt.title('Identified model with data')
plt.show()

#
# Write the fitted model to SBML/bioscrape
fit_model.export_sbml('sbml_fit.xml')
# Optionally, write a bioscrape file.
fit_model.export_bioscrape('bs_fit.xml')