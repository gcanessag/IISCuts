from gurobipy import *


def check(u, ub, A, b, T, r, xs, c, checkm, ep):
    m = checkm.copy()
    m.addConstr(quicksum(c[i] * m.getVars()[i-1] for i in tuplelist(c)) <= ub - ep, 'UB')
    for i in m.getConstrs():
        k = int("".join(x for x in i.ConstrName if x.isdigit()))
        if i.ConstrName[0] == 's' and k in u:
            m.remove(i)
    return m
