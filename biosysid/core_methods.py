from biosysid.Model import Model
from biosysid.Data import Data

def import_sbml(filename):
    try:
        import libsbml
    except:
        print('Unable to import libsbml. Make sure that python-libsbml package is installed and working')
    reader = libsbml.SBMLReader()
    doc = reader.readSBML(filename)
    model = Model(filename, sbml = doc)
    params_dict = {}
    species_dict = {}
    # Get global parameters
    for param in doc.getModel().getListOfParameters():
        params_dict[param.getId()] = param.getValue()
    count = 0
    # Get local parameters
    for reaction in doc.getModel().getListOfReactions():
            for param in reaction.getKineticLaw().getListOfLocalParameters():
                count = count + 1
                params_dict[param.getId() + '_local_' + reaction.getId() + str(count)] = param.getValue()
    model.parameters = params_dict

    for species in doc.getModel().getListOfSpecies():
        species_dict[species.getId()] = species.getInitialAmount()

    model.species = species_dict
    return model

def import_timeseries(filename, time_column, value_column, properties = {}, plot_show = False):
    '''
    filename : csv file with columns for data values 
    (The column numbers start at 1)
    time_column : the column number in the file that has all the time series indexes that you want to import
    value_column : the column number in the file that has all the corresponding values that you want to import 
    properties : Optional dictionary to specify other properties that the imported data must satisfy. For example, 
    properties = {3 : 'abc'}, would only impor those rows that have column 3 value equal to 'abc'
    '''
    try:
        import csv
        from operator import itemgetter
        from itertools import groupby
        import math
    except:
        print('Packages not found. Make sure csv, operator, itertool, and math are installed.')

    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        data_dict = {}
        for row in csv_reader:
            if row and row[time_column - 1] and row[value_column - 1]:
                if properties:
                    for col, value in properties.items():
                        if row[col - 1] == str(value):
                            data_dict[float(row[time_column - 1])] = float(row[value_column - 1])
                        else:
                            break 
                else:
                    cell_t = row[time_column - 1]
                    cell_v = row[value_column - 1]
                    temp_str_t = cell_t.replace('.','',1).replace('e','',1).replace('-','',1)
                    temp_str_v = cell_v.replace('.','',1).replace('e','',1).replace('-','',1)
                    if temp_str_t.isdigit() and temp_str_v.isdigit():
                        data_dict[float(cell_t)] = float(cell_v)
        data_obj = Data(filename, 'timeseries', data_dict)
        time = list(data_obj.get_keys())
        values = list(data_obj.get_values())
        if plot_show:
            try:
                import matplotlib.pyplot as plt
            except:
                raise Exception('matplotlib not installed.')
            max_time = math.floor(max(time))
            index = []
            for i in range(len(time)):
                if int(math.floor(float(time[i]))) == max_time:
                    index.append(i)
            final_index = []
            for k, g in groupby(enumerate(index), lambda x:x[1]-x[0]):
                map_t = map(itemgetter(1), g)
                final_index.append(max(map_t)+1)
            init_time_index = 0
            for i in final_index:
                plt.plot(time[init_time_index:i],values[init_time_index:i])
                plt.show()
                init_time_index = i
    return data_obj