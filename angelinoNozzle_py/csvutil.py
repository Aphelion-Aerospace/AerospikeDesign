import csv
import pickle


def read_csv(filename):
    """ Function that reads csv data and returns data in an array

        INPUTS:
            filename (str): name of csv file to be read

        OUTPUTS:
            data (list): csv data

    """
    with open(filename,'r') as dest_f:
        data_iter = csv.reader(dest_f, delimiter = '    ',  quotechar = '"')
        data = [data for data in data_iter]

    return data

def read_pickle(filename):
    """ Function that reads pickled data and returns it in memory

        INPUTS:
            filename (str): name of pickle file to be read

        OUTPUTS:
            unserialized_data (?): data read from pickle file, can be nearly any form

    """
    with open(filename, 'rb') as handle:
        unserialized_data = pickle.load(handle)

    return unserialized_data

def save_pickle(filename,data):
    """ Function that saves data to a pickle file under the name filename

        INPUTS:
            filename (str): name data will be saved under

        OUTPUTS:
            None

    """
    with open(filename,'wb') as handle:
        pickle.dump(data,handle,protocol = pickle.HIGHEST_PROTOCOL)    


def table_to_npwheaders(data):
    """ Function reads data to a dictionary wear the keys are the headers of the table and the data values are floats

        INPUTS:
            data (list): list where headers are in the first row and all other data is type() = float

        OUTPUTS:
            data_dict (dict): dictionary where keys are header names and all data is type() = float
            
    """
    data = np.asarray(data)

    cols = data[0,:]

    data_dict = {}

    for i in range(len(cols)):
        data_column = [float(el) for el in data[1:,i]]
        data_dict[cols[i]] = data_column 

    return data_dict

## DEPRECATED
def strpickle_to_numpypickle(filename):
    data = read_pickle(filename)
    data = np.asarray(data)

    cols = data[0,:]


    cfd_dict = {}

    for i in range(len(cols)):
        data_column = [float(el) for el in data[1:,i]]
        cfd_dict[cols[i]] = data_column
    
    with open(filename[:-7]+'_np.pickle','wb') as handle:
        pickle.dump(cfd_dict,handle,protocol = pickle.HIGHEST_PROTOCOL)

def serialize_csv(filename):
    # read csv data
    with open(filename,'r') as dest_f:
        data_iter = csv.reader(dest_f, delimiter = '    ',  quotechar = '"')
        data = [data for data in data_iter]
    # serialize csv data
    with open(filename[:-4] + '.pickle','wb') as handle:
        pickle.dump(data, handle, protocol = pickle.HIGHEST_PROTOCOL)
