__author__ = 'arenduchintala'


import _numpypy.multiarray as np
import pdb, sys, codecs
from pprint import pprint as pp

#np.set_printoptions(precision=4, linewidth=180)


def display_best_alignment(ak, en, es):
    lk = len(en)
    mk = len(es)
    k_mat = np.zeros((mk, lk))
    for jk in range(lk):
        for ik in range(mk):
            k_mat[ik][jk] = delta[ak, ik, jk]
    print ' '.join(en)
    print ' '.join(es)
    for ik, max_jk in enumerate(np.argmax(k_mat, 1)):
        print ik, max_jk, corpus_target[ak][ik], corpus_source[ak][max_jk]


def parseargs(args):
    try:
        source_idx = args.index('-s')
        target_idx = args.index('-t')
        source = args[source_idx + 1]
        target = args[target_idx + 1]
        init_translations = args[args.index('-i') + 1]
        save_translations_learned = args[args.index('-p') + 1]
        save_alignment_out = args[args.index('-a') + 1]
        source_alignment_test = args[args.index('-as') + 1]
        target_alignment_test = args[args.index('-at') + 1]
        return source, target,init_translations, save_translations_learned, save_alignment_out, source_alignment_test, target_alignment_test
    except (ValueError, IndexError) as er:
        print 'Usage: python model1.py -t [train target] -s [train source] -i [initial translations] -p [save translations] -a [save alignment test] -as [alignment test source] -at [alignment test target]'
        exit()

if __name__ == "__main__":
    delta = {}
    q = {}
    translations = {}
    counts = {}
    source, target,init_translations, save_translations_learned, save_alignment_out, source_alignment_test, target_alignment_test = parseargs(sys.argv)

    corpus_source = open(source, 'r').readlines()
    corpus_target = open(target, 'r').readlines()

    """
    initialization
    """
    for line in open(init_translations, 'r').readlines():
        [fi, ej, p] = line.split()
        translations[fi, ej] = float(p)

    #print 'initial t:'
    #pp(translations)

    for k, sentence_en in enumerate(corpus_source):
        sentence_es = corpus_target[k]
        tokens_source = sentence_en.split()
        tokens_source.insert(0, 'NULL')
        tokens_target = sentence_es.split()
        corpus_source[k] = tokens_source
        corpus_target[k] = tokens_target
        mk = len(tokens_target)
        lk = len(tokens_source)
        for e, ej in enumerate(tokens_source):
            for f, fi in enumerate(tokens_target):
                q[e, f, lk, mk] = 1.0 / lk

    #print 'initial q:'
    #pp(q)
    """
    EM iterations
    """

    for iter in range(10):
        print "iteration",iter
        counts = dict.fromkeys(counts.iterkeys(), 0.0)
        for k, tokens_source in enumerate(corpus_source):
            #print iter, k, len(delta), len(translations)
            sys.stderr.write('iteration: %d sentence %d len delta %d len translations %d\r' % (iter, k, len(delta), len(translations)))
            sys.stderr.flush()
            tokens_target = corpus_target[k]
            mk = len(tokens_target)
            lk = len(tokens_source)
            qt_mat = np.zeros((mk, lk))
            #print t_mat, t_mat.shape
            for j in range(0, lk):
                for i in range(0, mk):
                    qt_mat[i][j] = q[j, i, lk, mk] * translations[tokens_target[i], tokens_source[j]]
            qt_sum = qt_mat.sum(1)
            #print qt_mat, qt_sum
            for j in range(0, lk):
                for i in range(0, mk):
                    delta[k, i, j] = qt_mat[i][j] / qt_sum[i]
                    counts[tokens_target[i], tokens_source[j]] = counts.get((tokens_target[i], tokens_source[j]), 0.0) + delta[k, i, j]
                    counts[tokens_source[j]] = counts.get(tokens_source[j], 0.0) + delta[k, i, j]
                    counts[j, i, lk, mk] = counts.get((j, i, lk, mk), 0.0) + delta[k, i, j]
                    counts[i, lk, mk] = counts.get((i, lk, mk), 0.0) + delta[k, i, j]

        """
        update translations
        """
        for t_es_i, t_en_j in translations:
            translations[t_es_i, t_en_j] = counts[t_es_i, t_en_j] / counts[t_en_j]
        for qj, qi, ql, qm in q:
            q[qj, qi, ql, qm] = counts[qj, qi, ql, qm] / counts[qi, ql, qm]

    #clear counts to release memory
    counts = {}
    corpus_source ={}
    corpus_target = {}
    print "completed em iterations"
    writer = open(save_alignment_out, 'w')
    dev_source = open(source_alignment_test, 'r').readlines()
    dev_target = open(target_alignment_test, 'r').readlines()
    for dk in range(len(dev_source)):
        tokens_source = dev_source[dk].split()
        tokens_source.insert(0, 'NULL')
        tokens_target = dev_target[dk].split()
        l = len(tokens_source)
        m = len(tokens_target)
        for i, token_es in enumerate(tokens_target):
            max_p = 0.0
            max_j = 0.0
            for j, token_en in enumerate(tokens_source):
                if q[j, i, l, m] * translations[token_es, token_en] > max_p:
                    max_p = q[j, i, l, m] * translations[token_es, token_en]
                    max_j = j
            if max_j > 0:
                writer.write(str(dk + 1) + ' ' + str(max_j) + ' ' + str(i + 1) + '\n')
    writer.flush()
    writer.close()

