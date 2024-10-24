# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 15:13:41 2023

@author: au156185
"""
import argparse

from flowsim import flowsim

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description =
       'Runs script to compute transient response of hydrailic head or ' +
       'groundwater flux')
    parser.add_argument('--yaml', default='flowsim.yaml', 
                        help='YAML is name of yaml file')
    parser.add_argument('--log', default='flowsim.log', 
                        help='LOG is name of log file')
    
    args = vars(parser.parse_args())
    
    flowsim.run_model(args['yaml'], args['log'])
pass