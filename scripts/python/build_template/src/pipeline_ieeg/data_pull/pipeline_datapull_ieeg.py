# Standard library import
import os
import sys
import pickle
import getpass
import argparse
import subprocess

# User library import
import pipeline_ieeg.data_pull.check_data_repository as CDR
from CNT_research_tools.get_iEEG_data import get_iEEG_data

def main(args,lief_parent_dir="/cache/dev/bjprager/data/",user_dir='../../../../user_data/'):
    """
    Checks if data requested already exists in user data.
    If not, check if data is in lief cache. If so, return commands to download from lief.
    If not, download requested data from ieeg.org .

    Returns
    -------
    None.

    """
    
    # Figure out expected filename based on parameters
    if args.local_path == None:
        localfile = '%s%s-%d_%d/data.pickle' %(user_dir,args.dataset,args.start,args.duration)
    else:
        localfile = args.local_path

    # Check if data exists
    if os.path.exists(localfile):
        print("Reading local copy.")
        DF,fs = pickle.load(open(localfile,'rb'))
        return DF,fs
    else:
        
        # Make path to data through borel
        ifile = lief_parent_dir+'%s-%d_%d/data.pickle' %(args.dataset,args.start,args.duration)

        # Check lief for data or download data as needed
        if CDR.check_ieeg_data(ifile):
            print('Reading in Lief copy.')

            # Copy lief file to local environment
            try:
                os.mkdir('%s%s-%d_%d/' %(user_dir,args.dataset,args.start,args.duration))            
            except FileExistsError:
                pass
            subprocess.run(["rsync", "-avze", "ssh", "bjprager@borel.seas.upenn.edu:%s" %(ifile),"%s" %(localfile)])

            # Read in the data
            DF,fs = pickle.load(open(localfile,'rb'))
        else:
            
            print("Connecting to ieeg.org for data.")
            # Grab password if needed
            if not args.password:
                args.password = getpass.getpass()
            
            # Connet to ieeg api for data
            endtime = args.start+args.duration 
            DF,fs   = get_iEEG_data(args.user, args.password, args.dataset, args.start, endtime)
            
        # Save the data to the user_data directory
        try:
            os.mkdir('%s%s-%d_%d/' %(user_dir,args.dataset,args.start,args.duration))
        except FileExistsError:
                pass
        pickle.dump((DF,fs),open('%s%s-%d_%d/data.pickle' %(user_dir,args.dataset,args.start,args.duration),'wb'))
            
        return DF,fs

if __name__ == '__main__':

    # Command line options needed to obtain data.
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--user', required=True, help='username')
    parser.add_argument('-p', '--password', help='password (will be prompted if omitted)')
    parser.add_argument('--dataset', help='dataset name')
    parser.add_argument('--start', type=int, help='start offset in usec')
    parser.add_argument('--duration', type=int, help='number of usec to request')
    parser.add_argument('--local_path', default=None, type=str, help='Path to local data to ingest manually. Default=None.')
    args = parser.parse_args()
    
    DF,fs = main()
