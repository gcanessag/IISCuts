#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 14:12:37 2017

@author: confitao
"""

from gurobipy import *

def parser(file):
    data = open(file)
    
    A = {}
    c = {}
    
    f = 1

    i = 1
    j = 1
    f = 1
    C = 0
    R = 0
    
    line = data.readline().split()
    while line != []:
        if C:
            for r in line:
                if j <= ncols:
                    c[j] = int(r)
                j += 1
        elif R:
            if first_R:
                ncol_r = int(line[0])
                j = 1
                first_R = 0
            else:
                for r in line:
                    A[i,int(r)] = 1
                    j += 1
                if j > ncol_r:
                    first_R = 1
                    i += 1
            
        if f:
            nrows = int(line[0])
            ncols = int(line[1])
            f = 0
        elif line[0] == 'Costs':
            C = 1
        elif line[0] == 'Rows':
            C = 0
            R = 1
            first_R = 1
        
        line = data.readline().split()
    
    return(A, c, ncols, nrows)