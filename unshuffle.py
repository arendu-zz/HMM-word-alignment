e = open("hansards.e").read().splitlines()
s = open("hansards.e.shuffle").read().splitlines()
a = open("alignments").read().splitlines()

o = open("alignments.shuffled", 'w')

indices = [s.index(w) for w in e]

for i in indices:
    o.write(str(i))
    #o.write(a[i])
    o.write('\n')
