# Modified from:
# https://stackoverflow.com/questions/27586404/how-to-efficiently-get-all-combinations-where-the-sum-is-10-or-below-in-python

import itertools

def parts(a, n, m, x):
    """ Knuth Algorithm H, Combinatorial Algorithms, Pre-Fascicle 3B
        Finds all partitions of n having exactly m elements.
        An upper bound on running time is (3 x number of
        partitions found) + m.  Not recursive!      
    """
    while (1):
        x.append(a[1:m+1])
        while a[2] < a[1]-1:
            a[1] -= 1
            a[2] += 1
            x.append(a[1:m+1])
        j=3
        s = a[1]+a[2]-1
        while a[j] >= a[1]-1:
            s += a[j]
            j += 1
        if j > m:
            break
        z = a[j] + 1
        a[j] = z
        j -= 1
        while j>1:
            a[j] = z
            s -= z
            j -= 1
            a[1] = s

def distinct_perms(partition):
    """ Aaron Williams Algorithm 1, "Loopless Generation of Multiset
        Permutations by Prefix Shifts".  Finds all distinct permutations
        of a list with repeated items.  I don't follow the paper all that
        well, but it _possibly_ has a running time which is proportional
        to the number of permutations (with 3 shift operations for each  
        permutation on average).  Not recursive!
    """

    perms = []
    val = 0
    nxt = 1
    l1 = [[partition[i],i+1] for i in range(len(partition))]
    l1[-1][nxt] = None
    #print(l1)
    head = 0
    i = len(l1)-2
    afteri = i+1
    tmp = []
    tmp += [l1[head][val]]
    c = head
    while l1[c][nxt] != None:
        tmp += [l1[l1[c][nxt]][val]]
        c = l1[c][nxt]
    perms.extend([tmp])
    while (l1[afteri][nxt] != None) or (l1[afteri][val] < l1[head][val]):
        if (l1[afteri][nxt] != None) and (l1[i][val]>=l1[l1[afteri][nxt]][val]):
            beforek = afteri
        else:
            beforek = i
        k = l1[beforek][nxt]
        l1[beforek][nxt] = l1[k][nxt]
        l1[k][nxt] = head
        if l1[k][val] < l1[head][val]:
            i = k
        afteri = l1[i][nxt]
        head = k
        tmp = []
        tmp += [l1[head][val]]
        c = head
        while l1[c][nxt] != None:
            tmp += [l1[l1[c][nxt]][val]]
            c = l1[c][nxt]
        perms.extend([tmp])

    return perms

# Find all partitions of length p or less adding up
# to maxn or less
def getMonom(maxn, p):

    # Special cases (Knuth's algorithm requires n and m >= 2)
    x = [[i] for i in range(maxn+1)]

    # Main cases: runs parts fn (maxn^2+maxn)/2 times
    for i in range(2, maxn+1):
        for j in range(2, min(p+1, i+1)):
            m = j
            n = i
            a = [0, n-m+1] + [1] * (m-1) + [-1] + [0] * (n-m-1)
            parts(a, n, m, x)
    y = []

    # For each partition, add zeros if necessary and then find
    # distinct permutations.  Runs distinct_perms function once
    # for each partition.
    for part in x:
        if len(part) < p:
            y += distinct_perms(part + [0] * (p - len(part)))
        else:
            y += distinct_perms(part)

    return y

