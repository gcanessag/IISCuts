from math import *
from gurobipy import *


def leaf(node):
    m1 = node.model.copy()
    m2 = node.model.copy()

    n = len(node.xs) - 1

    if node.flag != 'Inf':
        v1 = m1.getVars()
        m1.addConstr(v1[node.fi + n] <= floor(node.zs[node.fi].X), 'c_ub_%i' % node.lvl)
        v2 = m2.getVars()
        m2.addConstr(v2[node.fi + n] >= ceil(node.zs[node.fi].X), 'c_lb_%i' % node.lvl)
    
    m1.update()
    m2.update()
    return m1, m2
