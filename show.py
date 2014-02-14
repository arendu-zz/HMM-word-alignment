import re
def dot(e, f, a):
    #e = e.split()
    #f = f.split()

    

    s = ""
    s += 'graph G{\ngraph [splines=line, rankdir=TB]\nnode[color=white, shape=box, height=.3]\n\n'
    s += "subgraph  {rank=same; "
    for i in range(len(e)):
        s += 'e%02d[label="%s"];' % (i, e[i])
    s += "}\n"
    s += "subgraph  {rank=same; "
    for i in range(len(f)):
        s += 'f%02d[label="%s"];' % (i, f[i])
    s += "}\n"

    #for i in range(m):
    #    s += 'subgraph {rank=same; f%02d[label="%s"]; e%02d[label="%s"]}\n' % \
    #            (i, f[i], i, e[i])

    s += 'subgraph E {\nedge [style=invisible]\n'
    for i in range(len(e)-1):
        j = i + 1
        s += "e%02d--e%02d;\n" % (i, j)
    s += "}\n\n"
    s += 'subgraph F {\nedge [style=invisible]\n'
    for i in range(len(f)-1):
        j = i + 1
        s += "f%02d--f%02d;\n" % (i, j)
    s += "}\n\n"

    s += "edge [tailport=s, headport=n]\n"
    #for x in a.split():
    for x in a:
        t = re.search("-|\?", x).group()[0]
        [i, j] = [r.strip() for r in re.split("-|\?", x)]
        style =   " [style=dotted]" if t == "?" else  ""
        s +=  "e%02d--f%02d%s\n" % (int(j), int(i), style)

    s += "}"
    return s
    


if __name__ == "__main__":
    e = "the white house".split()
    f = "la casa blanca".split()
    a = "0-0 1?2 2-1".split()
    print dot(e, f, a)

