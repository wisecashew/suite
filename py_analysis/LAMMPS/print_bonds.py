#!/home/satyend/.conda/envs/phase/bin/python

import lmp_data
import lmp_settings
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Add the impropers to a pnipam datafile.")
parser.add_argument("--datafile",      dest='df',    type=str, action='store', help="enter address of datafile.")
args = parser.parse_args()


if __name__=="__main__":

	my_data = lmp_data.Data(args.df)
	
	
