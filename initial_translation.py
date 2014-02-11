__author__ = 'arenduchintala'

"""
instead of setting translations to 1/n for initial translations, we set it by custom functions
in this case based on editdistance

Then the initial translations are saved to disk.
Model1 can load this saved file are begin doing iterations

"""
import sys
import editdistance as ed


def parseargs(args):
    try:
        source_idx = args.index('-s')
        target_idx = args.index('-t')
        source = args[source_idx + 1]
        target = args[target_idx + 1]
        save_initial_translations = args[args.index('-o') + 1]
        method = args[args.index('-m') + 1]
        return source, target, save_initial_translations, method
    except (ValueError, IndexError) as er:
        print 'Usage: python initialtranslations.py -t [train target] -s [train source] -o [save initial translations] -m [uniform/editdistance]'
        exit()


if __name__ == "__main__":
    initial_translation = {}
    translations = {}
    source, target, save, method = parseargs(sys.argv)
    corpus_source = open(source, 'r').readlines()
    corpus_target = open(target, 'r').readlines()

for k, sentence_source in enumerate(corpus_source):
    sentence_target = corpus_target[k]
    tokens_source = sentence_source.split()
    tokens_source.insert(0, 'NULL')
    tokens_target = sentence_target.split()
    corpus_source[k] = tokens_source
    corpus_target[k] = tokens_target
    for e in tokens_source:
        n_e = initial_translation.get(e, set()) # making a set of all possible target tokens that appear with source
        n_e.update(tokens_target) # adding all tokens of the target sentence a potential translations for source token 'e'
        initial_translation[e] = n_e  #saving in a map
if method == 'uniform':
    for k, v in initial_translation.iteritems():  # walking through the map and setting initial translation probability uniformaly.
        for v_es in v:
            translations[v_es, k] = 1.0 / len(v)
            #print 'initial t:'
            #pp(translations)
else:
    """
    What if we dont set the initial translation  probabilities uniformly?
    look at: http://research.microsoft.com/pubs/150581/acl11.pdf
    """
    add_delta = 1.0
    for k, v in initial_translation.iteritems():
        print 'setting intial for ', k
        edr_k = map(lambda t: ed.edratio(t, k) + add_delta, v)
        sum_edr = sum(edr_k)
        for v_es, edr_es in zip(v, edr_k):
            translations[v_es, k] = edr_es / sum_edr

writer = open(save, 'w')
for k, v in translations.iteritems():
    writer.write(str(' '.join(k)) + '\t' + str(v) + '\n')
writer.flush()
writer.close()

