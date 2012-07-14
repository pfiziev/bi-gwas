import sys
import array
from ctypes import *

def find_potential_snps(matrix):
    total_snps = len(matrix)
    total_inds = len(matrix[0])
    potential_snps = pointer(c_int(0))
    total_potential_snps = c_int(0)


    return cdll.LoadLibrary("./libpreprocessing.so").find_potential_snps(matrix, total_snps, total_inds, )


