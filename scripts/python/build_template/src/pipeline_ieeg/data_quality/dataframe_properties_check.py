# Standard library import
import os
import sys

# User library import
import CNT_unit_tests as unit_tests

def main(DF,N,M=-1,dtype='float64',verbose=True):
    """
    Ensure dataframe properties are correct. Main function runs all checks.

    Parameters
    ----------
    DF : Pandas dataframe
        Dataframe contain ieeg data. SamplesxChannel format.
    N : Integer
        Number of channels expected.
    M : Integer, optional
        Number of expected samples. Set to -1 to skip this check. The default is -1.
    dtype : String, optional
        Expected datatype for the dataframe entries.
    verbose : TYPE, optional
        Print test results to stdout. The default is True.

    Returns
    -------
    Boolean flag.

    """
    
    # Run through all tests
    qflag  = True
    CLS    = unit_tests.array_unit_tests.TestArrayProperties(DF,verbose)
    qflag *= CLS.test_dim(M,N)
    qflag *= CLS.test_type(dtype)
    qflag *= CLS.test_inf()
    qflag *= CLS.test_nan()
    return qflag

if __name__ == '__main__':
    qflag = main()