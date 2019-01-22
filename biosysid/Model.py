import warnings
import numpy as np
import matplotlib.pyplot as plt

class Model(object):
    def __init__(self, name, parameters = {}, species = {}, sbml = None, bioscrape = None, ode = None):
        '''
        name : string representing the name of the model
        parameters : dictionary consisting of all model parameter names and values
        species : dictionary consisting of all output species names and amounts? TODO (decide)
        '''
        self.name = name
        self.parameters = parameters
        self.species = species
        if sbml:
            self.sbml = sbml
        elif bioscrape:
            self.bioscrape = bioscrape
        elif ode:
            self.ode = ode
        else:
            raise ValueError('Only SBML, bioscrape, or Python ODE models are accepted.')
        
    def run_mcmc(self, params, priors, timepoints, nwalkers, nsteps, nsamples, measurements, exp_data,
                simulator = 'bioscrape', sim_type = 'stochastic', cost = 'L2norm', penalty = 1, plot_show = False):
        '''

        cost: L2norm for L2 norm or inf for infinity norm cost
        '''
        def log_prior(param_dict, priors):
            for key,value in param_dict.items():
                range = priors[key]
                if value > max(range) or value < min(range):
                    return False
            return True
        
        def log_likelihood(log_params):
            param_dict = {}
            params_exp = np.exp(log_params)
            for key, p in zip(self.parameters.keys(),params_exp):
                param_dict[key] = p
            # Check prior
            if log_prior(param_dict, priors) == False:
                return -np.inf

            if simulator == 'bioscrape':
                try:
                    import bioscrape
                except:
                    print('bioscrape package must be installed to use bioscrape simulator')
                if self.sbml:
                    try:
                        import libsbml
                    except:
                        print('python-libsbml must be installed to use SBML models')
                    filename = 'temp_simulate.xml'
                    libsbml.writeSBML(self.sbml,filename)
                    self.bioscrape = bioscrape.types.read_model_from_sbml(filename)
                outputs = []
                for species in measurements:
                    outputs.append(self.bioscrape.get_species_index(species))
                results = np.zeros((len(timepoints), len(outputs), nsamples))
                for sample in range(nsamples):
                    sim, m = self.simulate_bioscrape(timepoints, sim_type)
                    for i in range(len(outputs)):
                        out = outputs[i]
                        results[:,i,sample] = sim[:,out]
            else:
                raise ValueError('Other simulators not implemented.')

            total_error = 0
            for i in range(len(outputs)):
                for j in range(nsamples):
                    d1 = results[:,i,j]
                    diff = np.abs(d1 - exp_data.get_values()) 
                    if cost == 'inf':
                        infinity_error = np.max(diff)
                        total_error += infinity_error**2
                    elif cost == 'L2norm':
                        L2_norm_error = diff**2
                        L2_norm_error = np.linalg.norm(diff)
                        total_error += L2_norm_error
            return -total_error*penalty
        
        # Run emcee
        try:
            import emcee
        except:
            print('emcee package not installed')
        ndim = len(self.parameters)
        p0 = []
        for walker in range(nwalkers):
            plist = []
            ploglist = []
            for key, value in self.parameters.items():
                pinit = np.random.normal(value, 0.5*value)
                plist.append(pinit)
                ploglist.append(np.log(pinit))
            p0.append(np.array(plist))
            # print(p0)
        # print('going to run emcee now')    
            print('Sample log-like: {0}'.format(log_likelihood(np.array(ploglist))))

        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_likelihood)
        # for i, junk in enumerate(sampler.sample(p0, iterations=nsteps)):
        #     print('Step %d' % i)
        sampler.run_mcmc(p0, nsteps)    
        # Write results
        import csv
        with open('mcmc_results.csv','w', newline = "") as f:
            writer = csv.writer(f)
            writer.writerows(sampler.flatchain)
            # f.write(str(sampler.flatchain))
            f.close()
            
        print('Successfully completed MCMC parameter identification procedure.')

        best_p = []
        for i in range(len(self.parameters)):
            my_list = [tup[i] for tup in sampler.flatchain]
            new_list = []
            for x in my_list:
                if x > 0:
                    new_list.append(x)
            if plot_show:
                n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
            else:
                fig = plt.figure()
                n, bins, patches = plt.hist(new_list, density = True, histtype = "bar")
                plt.close(fig)
            # n, bins, patches = plt.hist(new_list, density = True, bins = 80, histtype = "bar")
            # Find best p
            best_p_ind = np.where(n == np.max(n))
            best_p.append(bins[best_p_ind])
            # Plot
            if plot_show:
                plt.savefig('parameter - ' + str(list(self.parameters.keys())[i]) +' .svg')
                plt.show()

        # Write fitted model
        best_p = list(best_p)
        fitted_model = self
        params_names = list(fitted_model.parameters.keys())
        params = {}
        for i in range(len(params_names)):
            p_name = params_names[i]
            p_sampled_value = best_p[i]
            params[p_name] = p_sampled_value

        fitted_model.parameters = params

        # Simulate again
        fitted_model.simulate_bioscrape(timepoints, 'stochastic', species_to_plot = ['c1'], plot_show = plot_show)
        return fitted_model, params


    def simulate_bioscrape(self, timepoints, type = 'deterministic', species_to_plot = [], plot_show = False):
        ''' 
        To simulate using bioscrape.
        Returns the data for all species and bioscrape model object which can be used to find out species indexes.
        NOTE : Needs bioscrape package installed to simulate. 
        TODO : species_to_plot not implemented. 
        TODO : Returns result and model
        '''
        try:
            import bioscrape
        except:
            print('bioscrape is not installed.')

        if self.sbml:
            try:
                import libsbml
            except:
                print('libsbml-python must be installed to use SBML models')
            # If SBML model,
            filename = 'temp_simulate.xml'
            libsbml.writeSBML(self.sbml, filename) 
            m = bioscrape.types.read_model_from_sbml(filename)
            s = bioscrape.simulator.ModelCSimInterface(m)
            if type == 'deterministic':
                s.py_prep_deterministic_simulation()
                s.py_set_initial_time(timepoints[0])
                sim = bioscrape.simulator.DeterministicSimulator()
                result = sim.py_simulate(s, timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            elif type == 'stochastic':
                warnings.warn('For stochastic simulation of SBML models using bioscrape, it is highly recommended to NOT use reversible reactions as the SSA algorithm might not work for such cases.')
                sim = bioscrape.simulator.SSASimulator()
                s.py_set_initial_time(timepoints[0])
                result = sim.py_simulate(s,timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            else:
                raise ValueError('Optional argument "type" must be either deterministic or stochastic')

        elif self.bioscrape:
            # If bioscrape model
            m = self.bioscrape
            s = bioscrape.simulator.ModelCSimInterface(m)
            if type == 'deterministic':
                s.py_prep_deterministic_simulation()
                s.py_set_initial_time(timepoints[0])
                sim = bioscrape.simulator.DeterministicSimulator()
                result = sim.py_simulate(s, timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            elif type == 'stochastic':
                warnings.warn('For stochastic simulation of SBML models using bioscrape, it is highly recommended to NOT use reversible reactions as the SSA algorithm might not work for such cases.')
                sim = bioscrape.simulator.SSASimulator()
                s.py_set_initial_time(timepoints[0])
                result = sim.py_simulate(s,timepoints)
                result = result.py_get_result()
                if plot_show:
                    for species in species_to_plot:
                        ind = m.get_species_index(species)
                        plt.plot(timepoints,result[:,ind])
                    plt.title(str(species_to_plot) + ' vs time')
                    plt.show()
                return result, m
            else:
                raise ValueError('Optional argument "type" must be either deterministic or stochastic')
            
    def simulate_roadrunner(self, timepoints, species_to_plot = []):
        ''' 
        To simulate using roadrunner.
        Returns the data for all species and bioscrape model object which can be used to find out species indexes.
        NOTE : Needs roadrunner package installed to simulate. 
        TODO : species_to_plot not implemented. 
        TODO : bioscrape.convert_to_sbml not implemented (possibly available in later versions of bioscrape)
        '''
        try:
            import roadrunner
        except:
            print('roadrunner is not installed.')

        if self.sbml:
            try:
                import libsbml
            except:
                print('libsbml-python must be installed to use SBML models')
            filename = 'temp_simulate.xml'
            libsbml.writeSBML(self.sbml, filename) 
            rr = roadrunner.RoadRunner(filename)
            if species_to_plot:
                rr.timeCourseSelections = ['time', species_to_plot]
            result = rr.simulate(timepoints[0],timepoints[-1],len(timepoints))
            res_ar = np.array(result)
            return res_ar[:,0],res_ar[:,1]
        elif self.bioscrape:
            # If bioscrape model
            '''
            TODO
            '''

            # sbml = bioscrape.types.convert_to_sbml()
            # self.sbml = sbml
            self.simulate_roadrunner(timepoints, species_to_plot)


    def export_sbml(self, filename):
        try:
            import libsbml
        except:
            print('libsbml-python must be installed to use SBML models')
        if self.sbml:
            model = self.sbml.getModel()
            params = self.parameters
            for pid,pval in params.items():
                if isinstance(pval, (list, np.ndarray)):
                    pval = pval[0]
                model.getElementBySId(pid).setValue(float(pval))
            libsbml.writeSBML(self.sbml,filename)
        elif self.bioscrape:
            try:
                import bioscrape
            except:
                print('bioscrape package must be installed to use bioscrape models')
            model_doc = self.bioscrape.convert_to_sbml()

            libsbml.writeSBML(model_doc, filename)
        else:
            raise ValueError('Model must be SBML or bioscrape XML. Other models not supported.')
        
    def export_bioscrape(self, filename):
        try:
            import bioscrape
        except:
            print('bioscrape package must be installed to use bioscrape')
        if self.sbml:
            try:
                import libsbml
            except:
                print('libsbml-python must be installed to use SBML models')
            libsbml.writeSBML(self.sbml,'sbml_model.xml')
            try:
                self.bioscrape = bioscrape.types.read_model_from_sbml('sbml_model.xml')
            except:
                import os
                os.remove('sbml_to_txt_temp.xml')
                self.bioscrape = bioscrape.types.read_model_from_sbml('sbml_model.xml')
            print('bioscrape text XML written to sbml_to_txt_temp.xml')
        elif self.bioscrape:
            # self.bioscrape.write_xml(filename)
            # raise NotImplemented('In future, this will be implemented using using bioscrape API')
            print('Not implemented yet. Export to bioscrape')
        else:
            raise ValueError('Model must be SBML or biosrape XML. Other models not supported')