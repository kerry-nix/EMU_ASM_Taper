# Standard library import
import sys

# User requested tools for feature selection
from CNT_research_tools import clean_labels as CL

def main(DF,preprocesses=None):
    """
    Create a dictionary with feature selection for ieeg data. 

    Parameters
    ----------
    Dataframe : dataframe structure
        Dataframe with channel names as column headers.
    preprocesses : list, optional
        Preprocessing steps. Currently available:
            CL = Clean labels
        The default is None. If None, all preprocessing steps. Case Sensitive.
        Steps performed in order of provided list.
    channels : list, optional
        Channels to analyze.
        The default is None. If None, all channels.

    Returns
    -------
    Dictionary with requested features.

    """
    
    # Setup variables to go through functions in order
    if preprocesses == None:
        preprocessing_steps = ['CL']
    else:
        preprocessing_steps = preprocesses
    
    # Loop over preprocessing steps
    for ival in preprocessing_steps:
        if ival == 'CL':
            #ieeg_obj = ieeg_obj.copy()
            DF.columns = CL.clean_labels(DF)
    return DF

if __name__ == '__main__':
    
    main()