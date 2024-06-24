from smps import *
import sys

if __name__ == "__main__":
    lpfile = "vac250a"
    beta = .05
    
    [model, A, b, costs, sense, x] = corfile(lpfile + '.cor')
    
    [m, A, b, c, sense, x] = stofile(lpfile + '.sto', model, beta, sense, x)
    
    model.write(lpfile+'i.lp')
    print("Done!")
