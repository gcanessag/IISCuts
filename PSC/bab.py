#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 17:56:48 2017

@author: confitao
"""

from gurobipy import *
from numpy import *
from par import *
from time import *
import sys

def epsfun(eps, model, c, u):
    x = {}
    for v in model.getVars():
        if v.VarName[0] == 'x':
            x[int(v.VarName[1:])] = v
        else:
            e = v
    model.addConstr(u - quicksum(c[i] * x[i] for i in tuplelist(c)) <= e, 'LB_Obj')
    model.addConstr(u - quicksum(c[i] * x[i] for i in tuplelist(c)) >= eps, 'LB_eps')
    
    model.update()
    model.optimize()
    
    if model.status == GRB.status.OPTIMAL:
#        model.write("eps.lp")
        return(e.X)
    else:
        return(eps)
    
def iiscut(model, where):
    #try:
    if where == GRB.callback.MIPNODE:
        tnow = time() 
        u = []
        zs = {}
        xs = {}
        for i in tuplelist(model._x):
            xs[i] = model.cbGetNodeRel(model._x[i])
        for i in tuplelist(model._z):
            zs[i] = model.cbGetNodeRel(model._z[i])
        U = sorted(zs, key = zs.__getitem__, reverse = True)
        u = []
        for i in range(int(ceil(model._beta * len(zs)))):
            if U[i] >= model._tol:
                u.append(U[i])
        ub = model.cbGet(GRB.callback.MIPNODE_OBJBST)
        if ub > 0 and ub < model._bounds:
#            if int(model.cbGet(GRB.callback.MIPNODE_NODCNT)) > 0:
#                eps = 0
#            else:
            eps = epsfun(model._eps, model._epsm.copy(), model._c, model._bounds)
            model._epsarr.append(eps)
            newmodel = check(u, ub, model._checkm, model._c, eps)
#            newmodel.params.timelimit = 5
            newmodel.optimize()
            if newmodel.status == GRB.status.OPTIMAL:
                if newmodel.ObjVal < model._bounds:
                    for i in tuplelist(model._z):
                        if i not in u:
                            model.cbSetSolution(model._z[i], 0)
                        else:
                            model.cbSetSolution(model._z[i], 1)
                    for (i,r) in enumerate(newmodel.getVars()):
                        model.cbSetSolution(model._x[i+1], r.X)
                    model._bounds = newmodel.ObjVal
                    print('Adding new node, using epsilon = %.3f. New upper bound: %.2f.' % (eps, newmodel.ObjVal))
            elif newmodel.status == GRB.status.INFEASIBLE or newmodel.status == GRB.status.INF_OR_UNBD:
                if len(u) > 0 or len(u) <= int(ceil(model._beta * len(zs))):
                    if len(u) < int(ceil(model._beta * len(zs))):
                        for i in range(int(ceil(model._beta * len(zs)) - len(u2))):
                            u.append(0)
                    if tuple(u) not in tuplelist(model._iiss):
                        timeiis = time()
                        newmodel.computeIIS()
                        newcut = []
                        for c in newmodel.getConstrs():
                            if c.ConstrName[0] == 's' and c.IISConstr:
#                                k = int("".join(x for x in c.ConstrName if x.isdigit()))
                                newcut.append(int(c.ConstrName[5:]))
#                                newcut.append(k)
                        newcut = sorted(newcut,key=newcut.count,reverse=True)
                        #newcut = [newcut[i] for i in range(min(model._clength,len(newcut)))]
                        #newcut = set(newcut)
                        cut = []
                        i = 1
                        while len(cut) < min(model._cutlen,len(newcut)) and len(newcut) > 0:
                            if newcut[i] != newcut[i-1]:
                                cut.append(newcut[i-1])
                            if i < len(newcut):
                                i += 1
                            else:
                                break
                        newcut = cut
                        if len(newcut) > 0:
                            print('Adding cut: %s >= 1, using epsilon = %.3f.' % (quicksum(model._z[i] for i in newcut), eps))
                            model.cbLazy(quicksum(model._z[i] for i in newcut if i != 0) >= 1)
                        model._iiss[tuple(u)] = newcut
                    else:
                        if len(model._iiss[tuple(u)]) > 0:
                            print('Adding cut: %s >= 1, using epsilon = %.3f.' % (quicksum(model._z[i] for i in model._iiss[tuple(u)]), eps))
                            model.cbLazy(quicksum(model._z[i] for i in model._iiss[tuple(u)] if i != 0) >= 1)
        model._Stime += time() - tnow

def check(u, ub, checkm, c, ep):
    m = checkm.copy()
    m.addConstr(quicksum(c[i] * m.getVars()[i-1] for i in tuplelist(c)) <= ub - ep, 'UB')
    for i in m.getConstrs():
        if int(i.ConstrName[5:]) in u:
            m.remove(i)
    return(m)

if __name__ == "__main__":
    #Call: S N cut (cutlen)
    objs = []
    Alens = []
    nodes = []
    seed = int(sys.argv[1])
#    seeds = [40635.0, 65.0, 0.0, 18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, 0.0, 4.0, 5468.0, 0.0, 2266.0, 0.0, 98.0, 0.0, 0.0, 3249.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 604.0, 0.0, 0.0, 0.0, 7.0, 0.0, 0.0, 0.0, 0.0, 37769.0, 80.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6866.0, 0.0, 0.0, 0.0, 0.0, 0.0, 18839.0, 0.0, 0.0, 8547.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 244.0, 0.0, 0.0, 2298.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 32920.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 3447.0, 93.0, 73.0, 3864.0, 0.0, 0.0, 67251.0, 0.0, 0.0, 0.0, 38657.0, 0.0, 1050.0, 5.0, 0.0, 0.0, 0.0, 3.0, 17281.0, 0.0, 6611.0, 0.0, 0.0, 5784.0, 0.0, 619.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 13.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 49.0, 5180.0, 0.0, 3582.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 87.0, 0.0, 45076.0, 0.0, 0.0, 0.0, 3261.0, 255.0, 0.0, 0.0, 18320.0, 0.0, 0.0, 0.0, 8190.0, 0.0, 0.0, 0.0, 0.0, 27.0, 111.0, 0.0, 0.0, 0.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4264.0, 0.0, 6.0, 0.0, 0.0, 0.0, 23.0, 260.0]
#    for i,s in enumerate(seeds):
    random.seed(seed)
    file = "scp41.txt" #sys.argv[1]
    cut = int(sys.argv[3])
    if cut:
        cutlen = int(sys.argv[4])
    else:
        cutlen = 0
    N = int(sys.argv[2])
    beta = .1            
    q0 = .0001
    ind = 1
    
    print('Loading file %s...' % file)
    T, c, ncols, nrows = parser(file)
    if ind:
        print("Loading independent model...")
        p = [1 - q0 ** (1/i) for i in range(1,201)]
    else:
        print("Loading dependent model...")
        p = random.random(200)/10
    h = {}
    s = {}
    if ind:
        for n in range(1,N+1):
            i = 1
            Ys = random.binomial(1,p)
            for y in Ys:
                if i not in tuplelist(s):
                    s[i] = 0
                if y:
                    h[i,n] = y
                    s[i] += 1
                i += 1
#        else:
#            for n in range(1,N+1):
#                i = 1
#                for r in range(1,21):
#                    E = []
#                    Ys = random.binomial(1,p[20*(r-1):20*r])
#                    for j,y in enumerate(Ys):
#                        if j != len(Ys)-1:
#                            E.append(max(Ys[j],Ys[j+1]))
#                        else:
#                            E.append(max(Ys[j], Ys[0]))
#                    for e in E:
#                        if e:
#                            h[i,n] = e
#                        i += 1
    
    A = []
    A = [i for i in tuplelist(s) if s[i] > floor(beta*N)]
                    
    if cut:
        print('Loading PSP...')
        psp = Model()
        psp.params.outputflag = 0
        psp.params.logfile = ''
        psp.params.timelimit = 7200
        psp.params.threads = 8
        x2 = {}
        for i in tuplelist(c):
            x2[i] = psp.addVar(obj = c[i], name = "x%s" % i, vtype = GRB.BINARY)
        psp.update()
#        for n in range(1,N+1):
#            for i in range(1,nrows+1):
        for i in A:
            psp.addConstr(quicksum(x2[j] for (g,j) in tuplelist(T).select(i,'*')) >= 1, 'c%3.0f,%s' % (i,n))
        for (i,n) in tuplelist(h):
            if i not in A:
                psp.addConstr(quicksum(x2[j] for (g,j) in tuplelist(T).select(i,'*')) >= 1, 's%3.0f,%s' % (i,n))
        psp.update()
    #    psp.write('psp.lp')
    
    print('Loading model...')
    model = Model()
    model.params.logfile = ''
    model.params.seed = 1
    model.params.threads = 8
    model.params.timelimit = 7200
    model.params.cuts = 0
#    model.params.presolve = 0
    x = {}
    y = {}
    z = {}
    Nset = range(1,N+1)
    Rset = range(1,nrows+1)
#    for i in tuplelist(c):
#        x[i] = model.addVar(obj = c[i], name = "x%s" % i, vtype = GRB.BINARY)
#    for n in range(1,N+1):
#        z[n] = model.addVar(name = "z%s" % n, vtype = GRB.BINARY)
#        for i in range(1,nrows+1):
#            y[i,n] = model.addVar(name = "y%s,%s" % (i,n), vtype = GRB.BINARY)
    
    x = model.addVars(c, obj = c, name = 'x', vtype = GRB.BINARY)
    z = model.addVars(Nset, name = 'z', vtype = GRB.BINARY)
    #y = model.addVars(Rset,Nset, name = 'y', vtype = GRB.BINARY)
    model.update()
    
#    for n in  range(1,N+1):
#        for i in range(1,nrows+1):
    for i in A:
        model.addConstr(quicksum(x[j] for (g,j) in tuplelist(T).select(i,'*')) >= 1, 'c%3.0f,%s' % (i,n))
    for (i,n) in tuplelist(h):
        if i not in A:
            model.addConstr(quicksum(x[j] for (g,j) in tuplelist(T).select(i,'*')) + z[n] >= 1, 's%3.0f,%s' % (i,n))
#        model.addConstr(y[i,n] + z[n] >= b[i,n], 'M%3.0f,%s' % (i,n))
#    model.addConstrs((quicksum(x[j] for (g,j) in tuplelist(T).select(i,'*')) + z[n] >= b[i,n] for i in Rset for n in Nset), name = 's')
        
    #model.addConstrs((y[i,n] + z[n] >= b[i,n] for (i,n) in tuplelist(b)), 'M')
    model.addConstr(quicksum(z[n] for n in tuplelist(z)) <= floor(beta * N), 'knapsack')
    model.update()
    
    epsm = Model()
    epsm.params.threads = 1
    epsm.params.outputflag = 0
    epsm.params.logfile = ''
#    epsm.params.timelimit = 10
    x3 = {}
    for i in tuplelist(c):
        x3[i] = epsm.addVar(name = "x%s" % i, vtype = GRB.BINARY)
    e = epsm.addVar(name = 'eps')
    epsm.update()
    epsm.setObjective(e)
    for i in A:
        epsm.addConstr(quicksum(x3[j] for (g,j) in tuplelist(T).select(i,'*')) >= 1, 'c%3.0f,%s' % (i,n))
    epsm.update()
        
    eps = 1e-1
    model._Stime = 0
    rtime = time()
    if cut:
        model._tol = 1 - 1e-3
        model._beta = beta
        model._eps = 1e-1
        model._epsm = epsm
        model._c = c
        model._x = x
        model._y = y
        model._z = z
        model._epsarr = []
        model._rootnode = 0
        model._checkm = psp
        model.params.LazyConstraints = 1
        model._cutlen = cutlen
        model._bounds = float("inf")
        model._iiss = {}
        model.optimize(iiscut)
    else:
        model.optimize()
    Alens.append(len(A))
    rtime = time() - rtime
    if model.status == GRB.status.OPTIMAL:
        model.write('optimal.lp')
        model.write('optimal.sol')
    print('LOG: %s %s %s %s %s %s %s %s %s %s\n' % (seed, beta, eps, cut, cutlen, model.ObjVal, model.MipGap, model.NodeCount, model._Stime, rtime))
        #with open('log' + m + '_gurobi.txt', 'a') as tfile:
        #    tfile.write('%s %s %s %s %s %s %s %s %s\n' % (m, beta, eps, uc, xi, model.ObjVal, model.MipGap, model.NodeCount, model.RunTime))
    if cut:
        print(mean(model._epsarr), var(model._epsarr))
