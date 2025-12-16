# 3rd December 2025

import csv
import numpy as np

def open_csv(file_path: str):
    """
    Opens the .csv file and returns dictionary, where the first row of .csv are the dictionary headers and subsequent data is converted into list
    
    :param file_path: Path to the .csv file
    :type file_path: str
    """


    data = {}

    with open(file_path, newline='') as file:
        csv_data = csv.reader(file)


        for header in next(csv_data):   # initialize headers in the dictionary
            data[header] = []

        for row in csv_data:            # Iterate every row
            i=0
            for column in data:         # in each row iterate through each header
                data[column].append(float(row[i]))
                i += 1

    return data #

def detect_start(data: list):
    """
    Returns the index of the beggining of the burn
    
    :param data: data of interest (np. pressure)
    :type data: list
    """


    maxsofar = data[0]
    minsofar = data[0]
    avgsofar = data[0]
    stdevsofar = 9999999
    sum = 0
    i = 1
    sqsum = 0

    for item in data:
        if item > maxsofar:
            maxsofar = item
        if item < minsofar:
            minsofar = item
        sum += item
        avgsofar = sum/i

        sqsum = 0
        for i2 in range(i):
            sqsum += (data[i2] - avgsofar)**2
        stdevsofar = np.sqrt(sqsum/i)
        
        if (abs(maxsofar-avgsofar) - abs(minsofar-avgsofar) > 5*stdevsofar) and i>1:
            return i-1
        i += 1
    return -1

def detect_stop(data: list):
    """
    Returns the index of the end of the burn
    
    :param data: data of interest (np. pressure)
    :type data: list
    """
    reverse_data = []

    for i in range(len(data)):
        reverse_data.append(data[len(data)-1-i])

    return len(data)-1 - detect_start(reverse_data)


def downsample(data: list, number_of_samples: int):
    """
    the function reduces the number of elements in data (list) to match the specified number of elements (number_of_samples)
    
    :param data: Data to be reducet
    :type data: list
    :param number_of_samples: length of ouput list
    :type number_of_samples: int
    """

    len_data = len(data)
    skip = int(len_data/number_of_samples)
    new_data = []

    for i in range(number_of_samples):
        new_data.append(data[i*skip])

    return new_data


def average_smoothing(data: list, block_size: int = 10):
    """
    Apply average smoothing to the data
    
    :param data: Data to be smoothed
    :type data: list
    :param block_size: Number of samples to calculate average from
    :type block_size: int
    """
    output = []
    data_len = len(data)
    for i in range(data_len-1):
        if i+block_size > data_len-1:
            output.append(np.average(data[i:-1]))
        else:
            output.append(np.average(data[i:i+block_size]))
    
    output.append(data[-1])

    return output