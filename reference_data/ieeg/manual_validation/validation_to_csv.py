import os
import pandas as PD
from sys import argv

if __name__ == '__main__':
    
    # Read in the data
    DF_dict = PD.read_excel(argv[1],sheet_name=None)
    
    # Create directory structure as needed
    try:
        os.mkdir("DATA/")
    except FileExistsError:
        pass
    for isheet in DF_dict.keys():
        try:
            os.mkdir('DATA/%s/' %(isheet))
        except FileExistsError:
            pass
        
    # Loop over sheets and save their output to csv
    for isheet in DF_dict.keys():
        ofile = "DATA/%s/%s.csv" %(isheet,isheet)
        if not os.path.exists(ofile):
            DF_dict[isheet].to_csv(ofile)
        else:
            userin = input("%s already exists. Overwrite? (y/n)")
            if userin.lower() == 'y':
                DF_dict[isheet].to_csv(ofile)