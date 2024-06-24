from math import *
from gurobipy import *


def cuts(u, zs, ub, b, r, cutm, ep, cb=0, zval=None):
    if zval is None:
        zval = {}
    if u:
        cut = cutm.copy()
        y1 = {}
        y2 = {}

        for v in cut.getVars():
            if v.VarName[1] == '1':
                y1[int(v.VarName[3:])] = v
            elif v.VarName[1] == '2':
                y2[int(v.VarName[3:])] = v
            else:
                y3 = v

        for i in tuplelist(y2):
            if i in u:
                y2[i].ub = 0
            else:
                if cb:
                    y2[i].Obj = zval[i]
                else:
                    y2[i].Obj = zs[i].X

        cut.addConstr(
            quicksum(y1[i] * b[i] for i in tuplelist(b)) + quicksum(y2[i] * r[i] for i in tuplelist(y2)) + y3 * (
                        ub - ep) <= -1, 'c2')
        cut.update()
        cut.optimize()
        cut.write('cut.lp')
        expr = {}
        if cut.status == GRB.status.OPTIMAL:
            for v in tuplelist(y2):
                if y2[v].X != 0:
                    expr[v] = 1
        return expr
    else:
        return {}
