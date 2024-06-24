from time import *
from gurobipy import *
from smps import *
from check import *
from math import *
from numpy import *


def epsfun(eps, model, c, u):
    x = {}
    e = None
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
        # model.write("eps.lp")
        return e.X
    else:
        return eps


def bab(m, beta, uc=0, eps=1e-3):
    if uc:
        print('Loading SMPS files for %s instance, using beta = %s and base epsilon = %s.' % (m, beta, eps))
        [model, A, b, costs, sense, x] = corfile('datafiles/' + m + '.cor')

        epsm = model.copy()
        epsm.params.outputflag = 0
        epsm.params.logfile = ''
        #        epsm.params.threads = 4
        e = epsm.addVar(name='eps')
        epsm.update()
        epsm.setObjective(e)
        epsm.update()

        [sp, rw, sw, T, p, z] = stofile('datafiles/' + m + '.sto', model, beta, sense, x)
        print('Load complete!')

        for i in tuplelist(x):
            if i not in tuplelist(costs):
                costs[i] = 0

        print('Loading checking model precomputations.')
        checkm = Model('Checking Model')
        checkm.params.outputflag = 0
        checkm.params.logfile = ''
        checkm.params.seed = 1
        # checkm.params.threads = 4
        # checkm.params.MIPGap = 1
        # checkm.params.timelimit = 0.1
        x2 = checkm.addVars(tuplelist(x), name='x', vtype=GRB.BINARY)
        checkm.update()
        checkm.setObjective(quicksum(costs[i] * x2[i] for i in tuplelist(costs)))
        checkm.addConstrs(
            (quicksum(A[n, j] * x2[j] for (n, j) in tuplelist(A).select(i, '*')) == b[i] for i in tuplelist(b)), "c")
        checkm.addConstrs((quicksum(T[i, j] * x2[j] for j in tuplelist(x2)) <= rw[i] for i in tuplelist(rw)), "sc")
        #        for i in tuplelist(b):
        #            tA = [int(k[1]) for k in tuplelist(A).select(i,'*')]
        #            checkm.addConstr(quicksum(A[i,j] * x2[j] for j in tA) == b[i], 'c%s' %i)
        #        c_T = {}
        #        for i in tuplelist(rw):
        #            checkm.addConstr(quicksum(T[i,j] * x2[j] for j in tuplelist(x2) if (i,j) in tuplelist(T)) <= rw[i],
        #            'sc%s' % i)
        #            tT = [int(k[1]) for k in tuplelist(T).select(i,'*')]
        #            checkm.addConstr(quicksum(T[i,j] * x2[j] for j in tT) <= rw[i], 'sc%s' % i)
        checkm.optimize()
        print('Done!')

    else:
        print('Loading SMPS files for %s instance, using beta = %s.' % (m, beta))
        [model, A, b, costs, sense, x] = corfile(m + '.cor')
        [sp, rw, sw, T, p, z] = stofile(m + '.sto', model, beta, sense, x)
        print('Load complete!')
    model._x = x
    model._z = z
    model._A = A
    model._b = b
    model._T = T
    model._r = rw
    model._c = costs
    model._eps = eps
    model._tol = 1.0 - 1e-3
    model._beta = beta
    model.params.seed = 1
    # model.params.threads = 1
    model._epsarr = []
    model._Stime = 0
    if uc:
        # model.params.presolve = 0
        model.params.LazyConstraints = 1
        model._checkm = checkm
        model._epsm = epsm
        model._bounds = float("inf")
        # model._cutm = cutm
        model._iiss = {}
        rtime = time()
        model.optimize(iiscut)
    else:
        rtime = time()
        model.optimize()
    rtime = time() - rtime

    if model.status == GRB.status.OPTIMAL:
        model.write('optimal.lp')
        model.write('optimal.sol')
    print('LOG: %s %s %s %s %s %s %s %s %s\n' % (m, beta, eps, uc, model.ObjVal, model.MipGap, model.NodeCount,
                                                 model._Stime, rtime))
    # with open('log' + m + '_gurobi.txt', 'a') as tfile:
    #    tfile.write('%s %s %s %s %s %s %s %s %s\n' % (m, beta, eps, uc, xi, model.ObjVal, model.MipGap,
    #    model.NodeCount, model.RunTime))
    if uc:
        print(mean(model._epsarr), var(model._epsarr))


def iiscut(model, where):
    # try:
    if where == GRB.Callback.MIPSOL:
        u = []
        zs = {}
        xs = {}
        for i in tuplelist(model._x):
            xs[i] = model.cbGetSolution(model._x[i])
        for i in tuplelist(model._z):
            zs[i] = model.cbGetSolution(model._z[i])
        U = sorted(zs, key=zs.__getitem__, reverse=True)
        for i in range(int(ceil(model._beta * len(zs)))):
            if U[i] >= model._tol:
                u.append(U[i])
        ub = model.cbGet(GRB.Callback.MIPSOL_OBJBND)
        if 0 < ub < model._bounds:
            #            if int(model.cbGet(GRB.callback.MIPNODE_NODCNT)) > 0:
            #                eps = 0
            #            else:
            eps = epsfun(model._eps, model._epsm.copy(), model._c, model._bounds)
            model._epsarr.append(eps)
            tnow = time()
            newmodel = check(u, ub, model._A, model._b, model._T, model._r, model._x, model._c, model._checkm, eps)
            #            newmodel.params.timelimit = 5
            newmodel.optimize()
            #            print(newmodel.status)
            #            if newmodel.status != GRB.status.INFEASIBLE or newmodel.status != GRB.status.INF_OR_UNBD:
            #                print("Starting checking model optimization...")
            #                timemodel = time()
            #                newmodelUBF = check(u, ub, model._A, model._b, model._T, model._r, model._x, model._c,
            #                model._checkm, eps)
            #                newmodelUBF.optimize()
            #                print(newmodelUBF.status)
            #                newmodelUBF.write("IIS.lp")
            #                newmodelUBF.write("IIS.sol")
            #                print("Done! It took %s seconds." % (time()-timemodel))
            if newmodel.status == GRB.status.OPTIMAL:
                if newmodel.ObjVal < model._bounds:
                    for i in tuplelist(model._z):
                        if i not in u:
                            model.cbSetSolution(model._z[i], 0)
                        else:
                            model.cbSetSolution(model._z[i], 1)
                    for (i, r) in enumerate(newmodel.getVars()):
                        model.cbSetSolution(model._x[i + 1], r.X)
                    model._bounds = newmodel.ObjVal
                    print('Adding new node, using epsilon = %.3f. New upper bound: %.2f.' % (eps, newmodel.ObjVal))
            elif newmodel.status == GRB.status.INFEASIBLE or newmodel.status == GRB.status.INF_OR_UNBD:
                if len(u) > 0 or len(u) <= int(ceil(model._beta * len(zs))):
                    if len(u) < int(ceil(model._beta * len(zs))):
                        for i in range(int(ceil(model._beta * len(zs)) - len(u2))):
                            u.append(0)
                    #                    print(model._iiss)
                    #                    print(tuplelist(model._iiss))
                    #                    print(tuple(u))
                    #                    print(u)
                    if tuple(u) not in tuplelist(model._iiss):
                        #                        print("Searching for an IIS...")
                        timeiis = time()
                        newmodel.computeIIS()
                        #                        print("Done! It took %s seconds." % (time() - timeiis))
                        newcut = []
                        for c in newmodel.getConstrs():
                            if c.ConstrName[0] == 's' and c.IISConstr:
                                k = int("".join(x for x in c.ConstrName if x.isdigit()))
                                newcut.append(k)
                        if len(newcut) > 0:
                            #                            newmodel.write('is_model.lp')
                            #                            newmodel.write('is_model.ilp')
                            print('Adding cut: %s >= 1, using epsilon = %.3f.'
                                  % (quicksum(model._z[i] for i in newcut), eps))
                            # model.terminate()
                            model.cbLazy(quicksum(model._z[i] for i in newcut if i != 0) >= 1)
                            # with open('ncuts.txt', 'a') as tfile:
                            #    tfile.write('%s %s\n' % (len(newcut), newcut))
                        model._iiss[tuple(u)] = newcut
                    else:
                        if len(model._iiss[tuple(u)]) > 0:
                            print('Adding cut: %s >= 1, using epsilon = %.3f.'
                                  % (quicksum(model._z[i] for i in model._iiss[tuple(u)]), eps))
                            model.cbLazy(quicksum(model._z[i] for i in model._iiss[tuple(u)] if i != 0) >= 1)
            model._Stime += time() - tnow
    # except GurobiError as e:
    #     print('GurobiError: %s' % e.message)


if __name__ == "__main__":
    m = sys.argv[1]  # 'vac100a'
    beta = float(sys.argv[2])  # 0.05
    uc = int(sys.argv[3])
    bab(m, beta, uc)
