from .app import run

from .supportfunctions import (perform_deriv_fit, 
                                perform_drude_fit, 
                                perform_Rs_fit, 
                                perform_entire_prodecure,
                                manual_inflection)

def main():
    run()