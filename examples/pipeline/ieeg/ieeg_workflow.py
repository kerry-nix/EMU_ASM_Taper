# General imports
import sys
import argparse

# User library import
from pipeline_ieeg.data_pull import pipeline_datapull_ieeg as PDI
from pipeline_ieeg.data_quality import dataframe_properties_check as DPC
from pipeline_ieeg.preprocessing import pipeline_preprocessing_ieeg as PPI
from pipeline_ieeg.feature_selection import pipeline_feature_selection_ieeg as PFSI

def main():
    """
    Calls series of commands for a simple ieeg pipeline.

    Returns
    -------
    None.

    """
    
    # Command line options needed to obtain data.
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--user', required=True, help='username')
    parser.add_argument('-p', '--password', help='password (will be prompted if omitted)')
    parser.add_argument('--dataset', help='dataset name')
    parser.add_argument('--start', type=int, help='start offset in usec')
    parser.add_argument('--duration', type=int, help='number of usec to request')
    parser.add_argument('--local_path', default=None, type=str, help='Path to local data to ingest manually. Default=None.')
    parser.add_argument('--silent', dest='verbose', default=True, action='store_false', help='Silent Verbose Output. Default=False.')
    parser.add_argument('--nchan', type=int, help='Number of channels')
    args = parser.parse_args()
    
    # Data ingestion
    DF,fs = PDI.main(args)
    
    # Data quality check
    qflag = DPC.main(DF,16,verbose=args.verbose)
    if qflag and args.verbose:
        print("All data quality checks came back True. Proceeding to next step.")
    
    # Data preprocessing
    DF = PPI.main(DF)
    
    # Feature selection
    feature_dict = PFSI.main(DF,fs)
    
    return DF,fs,feature_dict

if __name__ == '__main__':
    DF,fs,feature_dict = main()
