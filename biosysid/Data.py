class Data(object):
    def __init__(self, name, type, data = {}):
        '''
        name : string representing the name of the data set 
        type : type of data set - whether time series (Use string 'timeseries') or distributional (use 'distrib')
        data : dictionary with key and value for the data
        '''
        self.name = name
        self.type = type
        self.data = data 

    def get_keys(self):
        '''
        Returns the key list of the data dictionary
        '''
        return list(self.data.keys())

    def get_values(self):
        '''
        Returns the values list of the data dictionary
        '''
        return list(self.data.values())
 