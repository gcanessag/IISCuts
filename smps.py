from gurobipy import *


def stofile(sto, model, beta, sense, x):
    sto = open(sto)
    l = sto.readline().split()
    p = {}
    z = {}
    T = {}
    r = {}
    s = {}
    js = {}
    while len(l) != 0:
        if l[0] == 'SC':
            i = int(l[1][4:])
            p[i] = float(l[3])
            z[i] = model.addVar(name='z%s' % (i), vtype=GRB.BINARY)
            l = sto.readline().split()
            while l[0] != 'SC' and l[0] != 'ENDATA':
                if l[0][0] == 'x':
                    j = int(l[0][1:])
                    T[i, j] = float(l[2])
                elif l[0] == 'RHS':
                    s[i] = int(l[1][1:])
                    r[i] = float(l[2])
                l = sto.readline().split()
        else:
            l = sto.readline().split()

    model.update()
    if r == {}:
        for i in tuplelist(z):
            r[i] = 1
            s[i] = max(tuplelist(sense))

    for i in tuplelist(r):
        if sense[s[i]] == 'E':
            model.addConstr(quicksum(T[i, j] * x[j] for j in tuplelist(x)) == r[i], 'sc%s' % (i))
        elif sense[s[i]] == 'L':
            model.addConstr(quicksum(T[i, j] * x[j] for j in tuplelist(x)) - 10 * z[i] <= r[i], 'sc%s' % (i))
        else:
            model.addConstr(quicksum(-T[i, j] * x[j] for j in tuplelist(x)) + 10 * z[i] <= -r[i], 'sc%s' % (i))

    model.addConstr(quicksum(p[i] * z[i] for i in tuplelist(z)) <= beta, 'knapsack')
    sto.close()
    return 1, r, s, T, p, z


def corfile(cor):
    name = cor
    cor = open(cor)
    l = cor.readline().split()
    sense = []
    x = {}
    m = Model(name)
    A = {}
    c = {}
    b = {}
    while len(l) != 0:
        if l[0] == 'ROWS':
            l = cor.readline().split()
            sense = {}
            while l[0] != 'COLUMNS':
                if l[1] != 'obj':
                    sense[int(l[1][1:])] = l[0]
                l = cor.readline().split()
        if l[0] == 'COLUMNS':
            l = cor.readline().split()
            while l[0] != 'RHS':
                i = int(l[0][1:])
                o = 0
                for j in range(len(l)):
                    if l[j][0] == 'c' and int(l[j][1:]) < 31:
                        A[int(l[j][1:]), i] = float(l[j + 1])
                    elif l[j] == 'obj':
                        # if i != io:
                        o = float(l[j + 1])
                        c[i] = o
                if i not in tuplelist(x):
                    x[i] = m.addVar(obj=o, name='x%s' % i, vtype=GRB.BINARY)
                l = cor.readline().split()
        if l[0] == 'RHS':
            while l[0] != 'BOUNDS' and len(l) != 0:
                l = cor.readline().split()
                for j in range(len(l)):
                    if l[j][0] == 'c' and int(l[j][1:]) < 31:
                        b[int(l[j][1:])] = float(l[j + 1])
        # if l[0] == 'BOUNDS':
        #    while l[0] != 'ENDATA':
        #        l.cor.readline().split()
        #        if l[0] == 'UP':
        l = cor.readline().split()

    m.update()

    #    for i in range(1,31):
    #        for j in tuplelist(x):
    #            if (i,j) in tuplelist(A):
    #                A[634+i,j] = -A[i,j]

    bs = tuplelist(sense)
    mb = max(bs)

    for i in bs:
        if i != mb:
            if sense[i] == 'E':
                m.addConstr(quicksum(A[i, j] * x[j] for j in tuplelist(x) if (i, j) in tuplelist(A)) == b[i],
                            'c%s' % i)
            elif sense[i] == 'L':
                m.addConstr(quicksum(A[i, j] * x[j] for j in tuplelist(x) if (i, j) in tuplelist(A)) <= b[i],
                            'c%s' % i)
            else:
                m.addConstr(quicksum(-A[i, j] * x[j] for j in tuplelist(x) if (i, j) in tuplelist(A)) <= -b[i],
                            'c%s' % i)
    m.update()
    cor.close()
    return m, A, b, c, sense, x
