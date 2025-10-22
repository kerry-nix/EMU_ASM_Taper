import pandas as PD
from glob import glob

def check_ieeg_data(infile,data_pointer_loc='../../../../data_pointers/ieeg/cached_ieeg_data.csv'):
    """
    Determines if the user requested ieeg file already exists in Lief.

    Parameters
    ----------
    infile : STR
        User provided identifier of the file to be downloaded.
    data_pointer_loc : STR, optional
        Absolute or relative path to data pointers. The default is '../../../data_pointers/cached_ieeg_data.csv'.

    Returns
    -------
    Boolean flag.

    """
    
    pointer_DF = PD.read_csv(data_pointer_loc)
    filelist   = pointer_DF['filepath'].values
    
    if infile in filelist:
        return True
    else:
        return False
    
def check_image_data(infile,data_pointer_loc='../../../data_pointers'):
    """
    Determines if the user requested imaging file already exists in Lief.

    Parameters
    ----------
    infile : STR
        User provided identifier of the file to be downloaded.
    data_pointer_loc : STR, optional
        Absolute or relative path to data pointers. The default is '../../../data_pointers'.

    Returns
    -------
    Boolean flag.

    """
    
    filelist = glob(data_pointer_loc+'/**/*fits',recursive=True)
    filelist = [ifile.split('/')[-1] for ifile in filelist]
    
    if infile in filelist:
        return True
    else:
        return False
    
def check_applewatch_data(infile,data_pointer_loc='../../../data_pointers'):
    """
    Determines if the user requested applewatch file already exists in Lief.

    Parameters
    ----------
    infile : STR
        User provided identifier of the file to be downloaded.
    data_pointer_loc : STR, optional
        Absolute or relative path to data pointers. The default is '../../../data_pointers'.

    Returns
    -------
    Boolean flag.

    """
    
    filelist = glob(data_pointer_loc+'/**/*json',recursive=True)
    filelist = [ifile.split('/')[-1] for ifile in filelist]
    
    if infile in filelist:
        return True
    else:
        return False