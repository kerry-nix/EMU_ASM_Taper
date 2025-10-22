# Standard library imports
import sys

# User requested tools for feature selection
from CNT_research_tools import bandpower as BP
from CNT_research_tools import line_length as LL

def main(DF,samp_freq,features=None,channels=None,bp_freq_min=60,bp_freq_max=120):
    """
    Create a dictionary with feature selection for ieeg data. 

    Parameters
    ----------
    DF : dataframe structure
        Dataframe to derive features from.
    features : list, optional
        Feature selection. Currently available:
            LL = Line Length,
            BP = Band Power
            T0 = Start time of seizure
            Spike
        The default is None. If None, all features. Case Sensitive.
    channels : list, optional
        Channels to analyze.
        The default is None. If None, all channels.

    Returns
    -------
    Dictionary with requested features.

    """
    
    # Initialize variables
    feature_dict = {}
    if channels == None:
        channel_list = DF.columns
    else:
        channel_list = channels
    if features == None:
        feature_list = ['LL','BP']
    else:
        feature_list = features
    
    # Work through by feature than by channel request
    if 'LL' in feature_list:
        feature_dict['LL'] = {}
        for ichannel in channel_list:
            feature_dict['LL'][ichannel] = LL.line_length(DF[ichannel].values)
    if 'BP' in feature_list:
        feature_dict['BP'] = {}
        for ichannel in channel_list:
            feature_dict['BP'][ichannel] = BP.bandpower(DF[ichannel].values,samp_freq,[bp_freq_min,bp_freq_max])    
    return feature_dict

if __name__ == '__main__':
    
    main()