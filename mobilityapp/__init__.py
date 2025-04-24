from mobilityapp.app import run

import sys

from mobilityapp.supportfunctions import (perform_deriv_fit, 
                                perform_drude_fit, 
                                perform_Rs_fit, 
                                perform_entire_prodecure,
                                manual_inflection,
                                compute_asym_uncertainties,
                                compute_mu_uncertainties)

def main():
    if len(sys.argv)>1:
        scaling=float(sys.argv[1])
    else:
        scaling=0
    run(scaling=scaling)