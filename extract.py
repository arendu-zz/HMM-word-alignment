"""
convert Adi's model1 etc format to Adam + Matt's format
right now it just reads the alignments from a file called 'u'; TODO proper fileinput
"""

fn = "u"
f = open(fn).read().splitlines()
l = [[int(n) for n in e.split()] for e in f]
m = [[s,e,f] for [s,e,f] in l if s <= 1000] # leq because of first line fencepost issue
z = [[(u-1, v-1) for [s, u, v] in m if s == i] for i in range(1000)]
for l in z[1:]:
    print " ".join(["%d-%d" % (u, v) for (u, v) in l])



