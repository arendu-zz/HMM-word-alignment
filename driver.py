import optparse
import sys
import show
from subprocess import call

optparser = optparse.OptionParser()
optparser.add_option("-d", "--data", dest="train", default="data/hansards", help="Data filename prefix (default=data)")
optparser.add_option("-e", "--english", dest="english", default="e", help="Suffix of English filename (default=e)")
optparser.add_option("-f", "--french", dest="french", default="f", help="Suffix of French filename (default=f)")
optparser.add_option("-a", "--alignments", dest="alignment", default="a", help="Suffix of gold alignments filename (default=a)")
optparser.add_option("-n", "--num_sents", dest="n", default=sys.maxint, type="int", help="Number of alignments to display")
(opts, args) = optparser.parse_args()
f_data = "%s.%s" % (opts.train, opts.french)
e_data = "%s.%s" % (opts.train, opts.english)
a_data = "%s.%s" % (opts.train, opts.alignment)
 
(opts, _) = optparser.parse_args()
bitext = [[sentence.strip().split() for sentence in pair] for pair in zip(open(f_data), open(e_data), open(a_data))[:opts.n]]

#for [e, f, a] in bitext[:36]:
for i in range(37):
    [f, e, a] = bitext[i]
    #print e, f, a
    o = open("temp", 'w')
    o.write(show.dot(e, f, a))
    o.close()
    of = open("vis/s%d.svg" % (i + 1), 'w')
    call(["dot", "-Tsvg", "temp"],stdout=of)
    of.close()
