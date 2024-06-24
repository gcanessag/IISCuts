from gurobipy import *

def check(u, ub, A, b, T, r, xs, c, checkm, ep):
    m = checkm.copy()
    m.addConstr(quicksum(c[i] * m.getVars()[i-1] for i in tuplelist(c)) <= ub - ep, 'UB')
    for i in m.getConstrs():
        if i.ConstrName[0] == 's' and i.ConstrName[2:] in u:
            m.remove(i)
    #x = {}
    #for i in tuplelist(xs):
    #    x[i] = m.addVar(obj = xs[i].Obj, name = 'x%s' % i, lb = -float('inf'))
    #m.update()
    #for i in tuplelist(b):
    #    m.addConstr(quicksum(A[i,j] * x[j] for j in tuplelist(x) if (i,j) in tuplelist(A)) <= b[i], 'c%s' % (i))
    #for i in tuplelist(r):
    #    if i not in n.u:
            #m.addConstr(quicksum(T[i,j] * x[j] for j in tuplelist(x) if (i,j) in tuplelist(T)) <= r[i], 'sc%s' % (i))
    #        m.addConstr(c_T[i] <= r[i], 'sc%s' % i)
    
    #for c in m.getConstrs():
    #    if c.ConstrName == 'knapsack' or c.ConstrName[0:2] == 'c_' or c.ConstrName == 'cut' or c.ConstrName[:
    #        m.remove(c)
    #for x in m.getVars():
    #    if x.VarName[0] == 'z':
    #        m.remove(x)

    #m.addConstr(quicksum(x.Obj * x for x in m.getVars() if x.VarName[0] == 'x') <= ub - n.r, '(3c)')
    #m.addConstr(quicksum(x for x in m.getVars() if x.VarName[0] =='z') <= 0, '(3d)')
    #c = node(m, n.lvl, n.r, n.sp, n.x)
    #e = quicksum(p[i] * n.zs[i] for i in tuplelist(n.zs))
    #print e.getValue(), beta
    #if e.getValue() <= beta:
    #    c = 0
    #else:
    #    c = 1
    return(m)
