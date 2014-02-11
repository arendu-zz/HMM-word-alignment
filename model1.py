__author__ = 'arenduchintala'

"""
parameters in model 1:
delta[k,i,j] = translation[ foreign[k,i], english[k,j]] / sum_j_0toL (translation( foreign[k,i], english[k,j]))
k = kth sentences in the corpus
i = ith word in the kth sentence in the foreign corpus
j = jth word in the kth sentence in the english corpus
L = total number of words in the kth sentence in the english corpus
M = total number of words in the kth sentence in the foreign corpus
"""
"""
counts in model 1:
count[ejk, fik] = count[ejk, fik] + delta[k,i,j]
count[ejk] = count[ejk] + delta[k,i,j]
"""
"""
translation
translation[f,e] = c(f,e) / c(e)
"""
"""
for debugging purposes
https://class.coursera.org/nlangp-001/forum/thread?thread_id=940#post-4052
"""
import sys
import pdb



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
        initial_translations = args[args.index('-i') + 1]
        save_translations_learned = args[args.index('-p') + 1]
        save_alignment_out = args[args.index('-a') + 1]
        source_alignment_test = args[args.index('-as') + 1]
        target_alignment_test = args[args.index('-at') + 1]
        return source, target, initial_translations, save_translations_learned, save_alignment_out, source_alignment_test, target_alignment_test
    except (ValueError, IndexError) as er:
        print 'Usage: python model1.py -t [train target] -s [train source] -i [initial translations] -p [save translations] -a [save alignment test] -as [alignment test source] -at [alignment test target]'
        exit()


if __name__ == "__main__":
    delta = {}
    translations = {}
    counts = {}
    source, target, init_translation, save_trans, ali_out, ali_source, ali_target = parseargs(sys.argv)
    corpus_source = open(source, 'r').readlines()
    corpus_target = open(target, 'r').readlines()
    pairs = [(['NULL'] + s.split(), t.split()) for s, t in zip(corpus_source, corpus_target)]
    init_translation = open(init_translation, 'r').readlines()
    #corpus_source = corpus_source[:100]
    #corpus_target = corpus_target[:100]

    for line in init_translation:
        [t, s, p] = line.split()
        translations[t, s] = float(p)

    """
    EM iterations
    """
    for iter in range(5):
        counts = dict.fromkeys(counts.iterkeys(), 0.0)
        #for k, source_sentence in enumerate(corpus_source):
        for k, (source_tokens, target_tokens) in enumerate(pairs):
            #print iter, k, len(delta), len(translations)
            sys.stdout.write('iteration: %d sentence %d len delta %d len translations %d\r' % (iter, k, len(delta), len(translations)))
            sys.stdout.flush()
            accu = []
            mat = []
            for i in range(0, len(target_tokens)):
                ac = 0.0
                mat_i = []
                for j in range(0, len(source_tokens)):
                    mat_i.append(translations[target_tokens[i], source_tokens[j]])
                    ac = ac + translations[target_tokens[i], source_tokens[j]]
                accu.append(ac)
                mat.append(mat_i)

            #print t_mat
            #print t_sum
            for j in range(0, len(source_tokens)):
                for i in range(0, len(target_tokens)):
                    delta[k, i, j] = mat[i][j] / accu[i]
                    counts[target_tokens[i], source_tokens[j]] = counts.get((target_tokens[i], source_tokens[j]), 0.0) + delta[k, i, j]
                    counts[source_tokens[j]] = counts.get(source_tokens[j], 0.0) + delta[k, i, j]
                    #print tokens_es[i], tokens_en[j], counts[tokens_es[i], tokens_en[j]]
                    #print tokens_en[j], counts[tokens_en[j]]
                    #print 'iteration:', iter, 'sentence', k
        """
        update translations
        """
        for target_i, source_j in translations:
            translations[target_i, source_j] = counts[target_i, source_j] / counts[source_j]

        """
        print 'iter', iter
        print 'delta:'
        pp(delta)
        print 'counts:'
        pp(counts)
        print 'translations:'
        pp(translations)
        """
        """
        check how the alignment looks for a particular training pair, a particular sentence
        """

        """display_best_alignment(1012, corpus_en[1012], corpus_es[1012])
        display_best_alignment(829, corpus_en[829], corpus_es[829])
        display_best_alignment(2204, corpus_en[2204], corpus_es[2204])
        display_best_alignment(4942, corpus_en[4942], corpus_es[4942])"""
    writer = open(save_trans, 'w')
    for k, v in translations.iteritems():
        writer.write(str(' '.join(k)) + '\t' + str(v) + '\n')
    writer.flush()
    writer.close()
    writer = open(ali_out, 'w')

    test_source = open(ali_source, 'r').readlines()
    test_target = open(ali_target, 'r').readlines()
    for dk in range(len(test_source)):
        source_tokens = test_source[dk].split()
        source_tokens.insert(0, 'NULL')
        target_tokens = test_target[dk].split()
        for i, token_target in enumerate(target_tokens):
            max_p = 0.0
            max_j = 0.0
            for j, token_source in enumerate(source_tokens):
                if translations[token_target, token_source] > max_p:
                    max_p = translations[token_target, token_source]
                    max_j = j
            if max_j > 0:
                writer.write(str(dk + 1) + ' ' + str(max_j) + ' ' + str(i + 1) + '\n')
    writer.flush()
    writer.close()

