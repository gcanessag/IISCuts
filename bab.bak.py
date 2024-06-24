from leaf import *
from cut import *
from time import *
from gurobipy import *
from smps import *
from check import *
from math import *

class node:
    def __init__(self, model, lvl, r, sp, x):
        self.lvl = lvl
        self.model = model
        self.model.params.outputflag = 0
        self.model.params.logfile = ""
        self.model.params.threads = 4
        self.model.optimize()
        self.zs = {}
        self.x = x
        self.xs = {}
        self.sp = sp
        self.u = []
        self.l = []
        self.r = r
        maxf = 0
        if model.status == GRB.status.OPTIMAL:
            self.model.write('leaf.lp')
            self.obj = model.objVal
            for v in model.getVars():
                i = int(v.VarName[1:])
                if sp and v.VarName[0] == 'z':
                    self.zs[i] = v
                    if (self.zs[i].X - floor(self.zs[i].X)) > maxf and i not in self.u + self.l:
                        maxf = self.zs[i].X - floor(self.zs[i].X)
                        self.fi = i
                else:
                    self.xs[i] = v
            self.flag = 'I'
            for v in model.getVars():
                if round(v.X - floor(v.X),r) > 0:
                    if self.sp and v.VarName[0] == 'z':
                        self.flag = 'C'
                        break
                    if not self.sp:
                        self.flag = 'C'
                        break
        else:
             self.flag = 'Inf'
             self.obj = float("inf")
        if sp:
            for v in tuplelist(self.zs):
                if round(self.zs[v].X,r) == 1 and int(self.zs[v].VarName[1:]) not in self.u:
                    self.u.append(int(self.zs[v].VarName[1:]))
                elif round(self.zs[v].X,r) == 0 and int(self.zs[v].VarName[1:]) not in self.l:
                    self.l.append(int(self.zs[v].VarName[1:]))

def bab(m, breadth, sp, beta, iis, eps = 1e-5, uc = 0, xi = 0):
    if sp:
        print('Loading SMPS files, using beta = %s and epsilon = %s.' % (beta, eps))
    if breadth:
        med = 'breadth first'
    else:
        med = 'depth first'
    
    [model, A, b, costs, sense, x] = corfile('vac100a.cor')

    if sp != '0':
        [sp, rw, sw, T, p, z] = stofile(m+'.sto', model, beta, sense, x)
    else:
        sp = 0
        T = {}
    print('Load complete!')

    if iis and uc:
        #Cut generating LP precumputations:
        #print 'Loading cut generating LP precalculations.'
        #cutm = Model('cut0')
        #cutm.params.outputflag = 0
        #cutm.params.logfile = ""
        #cutm.params.threads = 4
        #cutm.params.seed = 1
        #y1 = {}
        #y2 = {}
        #y3 = cutm.addVar(obj = 0, lb = 0, name = 'y3')
        #for i in tuplelist(b):
        #    y1[i] = cutm.addVar(obj = 0, lb = 0, name = 'y1_%s' % (i))
        #for i in tuplelist(rw):
        #    y2[i] = cutm.addVar(lb = 0, name = 'y2_%s' % i)
        #cutm.update()
        #for i in tuplelist(x):
        #    if i not in tuplelist(costs):
        #        costs[i] = 0
        #    cutm.addConstr(quicksum(y1[j] * A[j,i] for j in tuplelist(b) if (j,i) in tuplelist(A)) + quicksum(y2[j] * T[j,i] for j in tuplelist(rw) if (j,i) in tuplelist(T)) + y3 * costs[i] == 0, 'c1_%s' % (i))
        #    tA = [int(k[0]) for k in tuplelist(A).select('*',i)]
        #    tT = [int(k[0]) for k in tuplelist(T).select('*',i)]
        #    cutm.addConstr(quicksum(y1[j] * A[j,i] for j in tA) + quicksum(y2[j] * T[j,i] for j in tT) + y3 * costs[i] == 0, 'c1_%s' % i)

        #cutm.update()

        #Checking model precomputations:
        print('Loading checking model precomputations.')
        checkm = Model()
        checkm.params.outputflag = 0
        checkm.params.logfile = ""
        checkm.params.seed = 1
        checkm.params.threads = 4
        x2 = {}
        for i in tuplelist(x):
            x2[i] = checkm.addVar(obj = x[i].Obj, name = 'x%s' % i, lb = -float('inf'), vtype = GRB.BINARY)
        z2 = {}
        for i in tuplelist(z):
            z2[i] = checkm.addVar(ub = 0)
        checkm.update()
        for i in tuplelist(b):
            #checkm.addConstr(quicksum(A[i,j] * x2[j] for j in tuplelist(x2) if (i,j) in tuplelist(A)) <= b[i], 'c%s' % i)
            tA = [int(k[1]) for k in tuplelist(A).select(i,'*')]
            checkm.addConstr(quicksum(A[i,j] * x2[j] for j in tA) <= b[i], 'c%s' %i)
        c_T = {}
        for i in tuplelist(rw):
            #checkm.addConstr(quicksum(T[i,j] * x2[j] for j in tuplelist(x2) if (i,j) in tuplelist(T)) <= rw[i], 'sc%s' % i)
            tT = [int(k[1]) for k in tuplelist(T).select(i,'*')]
            checkm.addConstr(quicksum(T[i,j] * x2[j] for j in tT) <= rw[i], 'sc%s' % i)
        checkm.optimize()
        print('Done!')
    
    if uc == 2:
        it = 1
        id = 1
        r = 5
        nodes = {}
        nnodes = 1
        nodes[1,1] = node(model, it, r, sp, x)
        nodes[1,1].model.write('relaxed.lp')

        newcut = {}
        cut = 0
        
        relax = nodes[1,1].obj
        lb = nodes[1,1].obj
        print('LP Relaxation value = %s.' % lb)
        ub = float('inf')
        flag = 1
        start_time = time()
        now = time()
        now2 = time()
        print('')
        print('Starting B&B algorithm')
        print('Using %s method, using a single thread.' % med)
        if sp:
            print 'Only branching on z variables.'
        if iis:
            print 'Using IIS cuts with epsilon = %s.' % eps
        print '--------------------------------------------------------------------'
        print '        Nodes         |       Best      | Gap  |  Time'
        print ' explored |  left     |  int  | relaxed |      | (secs)'
        print '--------------------------------------------------------------------'
        while flag:
            newcut = {}
            c = 0
            if breadth:
                (I,J) = min(tuplelist(nodes))
                tu = tuplelist(nodes).select(I,'*')
            else:
                (I,J) = max(tuplelist(nodes))
                tu = tuplelist(nodes).select(I,'*')
            for (i,j) in tu:
                nnodes += 1
                if time() - now >= 5:
                    print ' %s\t\t%s  \t%s\t%s\t %s%%\t%s\t(%s n/s)' % (nnodes, len(nodes), round(ub,3), round(lb,3), gap, int(time()-start_time), round(nnodes/(time()-start_time),1))
                    now = time()
                if time() - now2 >= 60:
                    print '--------------------------------------------------------------------'
                    print '        Nodes         |       Best      | Gap  |  Time'
                    print ' explored |  left     |  int  | relaxed |      | (secs)'
                    print '--------------------------------------------------------------------'
                    now2 = time()
                if nodes[i,j].flag == 'I':
                    if ub > nodes[i,j].obj:
                        nodes[i,j].model.write('bab_inc.lp')
                        nodes[i,j].model.write('bab_inc.sol')
                        ub = nodes[i,j].obj
                        xs = nodes[i,j].xs
                        ls = nodes[i,j].lvl
                        if lb == -float("inf"):
                            gap = 100
                        elif lb == ub:
                            gap = 0
                        elif ub == 0:
                            gap = round(-lb * 100,1)
                        else:
                            gap = round((ub - lb) / ub * 100, 1)
                        print '*%s\t\t%s  \t%s\t%s\t %s%%\t%s\t(%s n/s)' % (nnodes, len(nodes), round(ub,3), round(lb,3), gap, int(time()-start_time), round(nnodes/(time()-start_time),1))
                    nodes[i,j].flag = 'A'
                elif nodes[i,j].flag == 'C':
                    c += 1
                    nodes[i,j].flag = 'A'
                    if nodes[i,j].obj < ub:
                        if iis and ub != float('inf'):
                            newmodel = check(nodes[i,j].u, ub, A, b, T, rw, x, costs, checkm, eps)
                            newnode = node(newmodel, nodes[i,j].lvl, r, sp, x)
                            if newnode.flag == 'Inf':
                                newcut = cuts(nodes[i,j].u, nodes[i,j].zs, ub, b, rw, cutm, eps)
                                if newcut != {} and newcut is not None:
                                    cut += 1
                                    for (k,l) in tuplelist(nodes).select(i,'*'):
                                        nodes[k,l].model.addConstr(quicksum(nodes[k,l].zs[v] for v in tuplelist(newcut)) >= 1, 'cut')
                                        nodes[k,l].model.update()
                                        nodes[k,l].model.optimize()
                                
                            else:
                                id += 1
                                nodes[j,id] = newnode
                        if nodes[i,j].flag == 'A':
                            m1, m2 = leaf(nodes[i,j])
                            id += 2
                            nodes[j,id] = node(m1, nodes[i,j].lvl + 1, r, sp, x)
                            nodes[j,id-1] =  node(m2, nodes[i,j].lvl + 1, r, sp, x)
                if lb == -float("inf"):
                    gap = 100
                elif lb == ub:
                    gap = 0
                elif ub == 0:
                    gap = round(-lb * 100,1)
                else:
                    gap = round((ub - lb) / ub * 100, 1)

            nodes = {(k,l): nodes[k,l] for (k,l) in nodes if nodes[k,l].flag != 'Inf' and  nodes[k,l].flag != 'A' and nodes[k,l].obj <= ub}

            it += 1
            
            if len(nodes) == 0 or gap <= 1e-6:
                flag = False
            
            if len(nodes) > 0:
                lb1 = nodes[min(nodes, key = lambda (k,l): nodes[k,l].obj)].obj
                if lb1 < float('Inf'):
                    lb = lb1

        if len(nodes) == 0:
            gap = 0


        if ub < float("inf"):
            print 'An Integer Optimum has been found at iteration %s!' % it
            print 'Level of Depth                        = %s.' % ls
            print 'LP Relaxation optimum objective value = %s.' % relax
            print 'MIP Integer optimum objective value   = %s.' % ub
            print 'Gap                                   = %s%%' % gap
            print 'Nodes explored                        = %s.' % nnodes
            print 'Number of IIS cuts produced           = %s.' % cut
            print 'Total process time                    = %s s.' % (round(time() - start_time,2))
            with open('log' + m + '.txt', 'a') as tfile:
                tfile.write('%s %s %s %s %s %s %s %s %s %s\n' % (m, breadth, beta, eps, iis, ub, gap, nnodes, cut, (round(time() - start_time,2))))
        else:
            print 'No Integer solution found after %s iterations.' % it

    else:
        for i in tuplelist(z):
            z[i].vtype = GRB.BINARY
        if xi:
            for i in tuplelist(x):
                x[i].vtype = GRB.BINARY
        model._x = x
        model._z = z
        model._A = A
        model._b = b
        model._T = T
        model._r = rw
        model._c = costs
        model._eps = eps
        model._tol = 1.0 - 1e-6
        model.params.seed = 1
        model.params.threads = 16
        if uc:
            model.params.presolve = 0
            model.params.LazyConstraints = 1
            model._checkm = checkm
            #model._cutm = cutm
            model.optimize(iiscut)
        else:
            model.optimize()

        if model.status == GRB.status.OPTIMAL:
            model.write('optimal.lp')
            model.write('optimal.sol')
            print 'LOG: %s %s %s %s %s %s %s %s %s\n' % (m, beta, eps, uc, xi, model.ObjVal, model.MipGap, model.NodeCount, model.RunTime)
            with open('log' + m + '_gurobi.txt', 'a') as tfile:
                tfile.write('%s %s %s %s %s %s %s %s %s\n' % (m, beta, eps, uc, xi, model.ObjVal, model.MipGap, model.NodeCount, model.RunTime))


def iiscut(model, where):
    try:
        if where == GRB.callback.MIPNODE:
            u = []
            zs = {}
            for i in tuplelist(model._z):
                zs[i] = model.cbGetNodeRel(model._z[i])
                if zs[i] >= model._tol:
                    u.append(i)
            ub = model.cbGet(GRB.callback.MIPNODE_OBJBST)
            newmodel = check(u, ub, model._A, model._b, model._T, model._r, model._x, model._c, model._checkm, model._eps)
            newmodel.optimize()
            if newmodel.status == GRB.status.OPTIMAL:
                print 'Adding new node.'
                for i in tuplelist(model._z):
                    if i not in u:
                        model.cbSetSolution(model._z[i], 0)
                    else:
                        model.cbSetSolution(model._z[i], 1)
            else:
                if len(u) > 0:
                    #Old method
                    #newcut = {}
                    #newcut = cuts(u, model._z, ub, model._b, model._r, model._cutm, model._eps, 1, zs)

                    #Using conflict refiner
                    newmodel.computeIIS()
                    newcut = []
                    for c in newmodel.getConstrs():
                        if c.ConstrName[0] == 's' and c.IISConstr:
                            newcut.append(int(c.ConstrName[2:]))
                    if len(newcut) > 0:
                        model.cbLazy(quicksum(model._z[i] for i in newcut) >= 1)
                        #with open('ncuts.txt', 'a') as tfile:
                        #    tfile.write('%s %s\n' % (len(newcut), newcut))
    except GurobiError as e:
        print 'GurobiError: %s' % e.message
