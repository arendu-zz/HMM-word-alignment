"""
import, convert, and merge fwd and rev alignments with a growing heuristic

WARNING this file is in a state of flux and it's very important to consider the e-f ordering at all stages


"""
import sys

num = 100000 # careful ...
def load_alignments(fn):
    f = open(fn)
    z = [[]]
    for line in f:
        s, e, f = [int(n) for n in line.split()]
        if len(z) < s:
            z.append([])
        z[s-1].append((e-1,f-1))
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
    return a ### !!!
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
    #a = final(a, m, n, u)
    return a


if __name__ == "__main__":
    ffile, rfile = sys.argv[1], sys.argv[2]
    fwd = [set(l) for l in flip(load_alignments(ffile))]
    sys.stderr.write("fwd loaded")
    rev = [set(l) for l in load_alignments(rfile)]
    sys.stderr.write("rev loaded")
    for i in range(len(fwd)):
        if (i % 1000 == 0): sys.stderr.write(".")
        f, r = fwd[i], rev[i]
        a = symmetrize(f, r)
        a = list(a)
        a.sort()
        print " ".join(["%d-%d" % (u, v) for (u, v) in a])
