a = open("hmm_77").read().splitlines()
i = open("index").read().splitlines()
index = [int(e) for e in i]


l = []
t = [(val, id) for (id, val) in enumerate(index)].sort()
for val, id in t[:1000]:
    l.append(a[val])

f = open("out", 'w')
f.write("\n".join([l]))
f.close()
