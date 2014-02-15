"""
import, convert, and merge fwd and rev alignments with a growing heuristic

WARNING this file is in a state of flux and it's very important to consider the e-f ordering at all stages


"""

num = 1000

def load_alignments(fn):
    f = open(fn).read().splitlines()
    l = [[int(n) for n in e.split()] for e in f]
    m = [[s,e,f] for [s,e,f] in l if s <= num] # leq because of first line fencepost issue
    z = [[(u-1, v-1) for [s, u, v] in m if s == i] for i in range(num+1)][1:]
    return z

def flip(l):
    return [[(v, u) for (u,v) in m] for m in l]

def aligned():
    pass

def final(a, m, n, u):
    for e in range(m):
        for f in range(n):
            if (not e in [j for (j, k) in a] or not f in [k for (j, k) in a]) and (e, f) in u:
                a.add((e, f))
    return a

def symmetrize(fwd, rev):
    a = fwd.intersection(rev)
    #return a
    u = fwd.union(r)
    neighbors= ((-1,-1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1))
    new = set([])
    m = max(max([e for (e, f) in fwd]), max([e for (f, e) in rev]))
    n = max(max([f for (e, f) in fwd]), max([f for (f, e) in rev]))
    while len(new) > 0:
        new = set()
        for e in range(m):
            for f in range(n):
                if (e, f) in a:
                    for (de, df) in neighbors:
                        ee, ff = e + de, f + df
                        if (not ee in [j for (j, k) in a] or not ff in [k for (j, k) in a]) and (ee, ff) in u:
                            new.add((ee, ff))
        a = a.union(new)
    a = final(a, m, n, u)
    return a


if __name__ == "__main__":
    fwd = [set(l) for l in flip(load_alignments("f"))]
    rev = [set(l) for l in load_alignments("r")]
    for i in range(len(fwd)):
        f, r = fwd[i], rev[i]
        a = symmetrize(f, r)
        a = list(a)
        a.sort()
        print " ".join(["%d-%d" % (u, v) for (u, v) in a])
