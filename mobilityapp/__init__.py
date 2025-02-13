from .app import run

import sys

from .supportfunctions import (perform_deriv_fit, 
                                perform_drude_fit, 
                                perform_Rs_fit, 
                                perform_entire_prodecure,
                                manual_inflection)

def main():
    if len(sys.argv)>1:
        scaling=float(sys.argv[1])
    else:
        scaling=0
    run(scaling=scaling)