import getpass
import argparse
import pandas as PD
from sys import argv
from ieeg.auth import Session

def ieeg():
    """
    Returns sample dataframe from ieeg.org for quick testing purposes.
    """
    
    # Command line options needed to obtain data.
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--user', required=True, help='username')
    parser.add_argument('-p', '--password', help='password (will be prompted if omitted)')
    parser.add_argument('--dataset', default='I004_A0003_D001', help='dataset name')
    parser.add_argument('--start', default=13090000, type=int, help='start offset in usec')
    parser.add_argument('--duration', default=100000, type=int, help='number of usec to request')
    args = parser.parse_args()

    if not args.password:
        args.password = getpass.getpass()

    with Session(args.user, args.password) as session:
        dataset_name = args.dataset
        dataset      = session.open_dataset(dataset_name)
        channels     = list(range(len(dataset.ch_labels)))
        raw_data     = dataset.get_data(args.start, args.duration, channels)
        labels       = dataset.get_channel_labels()
        session.close_dataset(dataset_name)
    return PD.DataFrame(raw_data,columns=labels)

def main():

    if argv[2] == 'imaging':
        pass
    elif argv[2] == 'applewatch':
        pass
    else:
        DF = ieeg()

if __name__ == "__main__":
    
    main()
