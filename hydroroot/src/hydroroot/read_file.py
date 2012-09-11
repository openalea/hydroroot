import numpy as np

def readCSVFile(filename):
    """
    Read and extract data from a csv file, supposed that the data is stored in 2 columns.

    :Parameters:
    -    `filename`    - file name to read

    :Returns:
    `data`    -    record array of (x, y) values, column headers recorded in dtype
    """
    data = np.recfromcsv(filename,delimiter=';')
    #fobj = file(filename, 'r')
    #readerObj = csv.reader(fobj, delimiter = '\t')
    #mainList = [list(row) for row in readerObj]
    #angles = [[int(angle) for angle in row if angle != "" ] for row in mainList]
    #fobj.close()
    return data