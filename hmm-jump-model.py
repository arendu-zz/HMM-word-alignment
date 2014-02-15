__author__ = 'arenduchintala'
#import _numpypy.multiarray as np
import logutils as lu
from math import log, exp, fabs
from pprint import pprint
import pdb, sys, codecs
from collections import Counter

#np.set_printoptions(precision=4, linewidth=180)

BOUNDRY_STATE = "###"
alignment_probs = {}
jump_counts = {}
jump_counts_iip = {}
translations_probs = {}


def jump_key(j1, j0):
    if j1 == BOUNDRY_STATE or j0 == BOUNDRY_STATE:
        return 'count', j1, j0
    else:
        return 'count', j1 - j0


def jump_iip_key(L, i_prime):
    return 'I_i_prime', L, i_prime


def get_trelis(target_seq, source_seq):
    target_trelis = [(i - 1, w) for i, w in enumerate(source_seq) if w != BOUNDRY_STATE]
    trelis = [[(BOUNDRY_STATE, BOUNDRY_STATE)]] + [target_trelis] * (len(target_seq) - 2) + [[(BOUNDRY_STATE, BOUNDRY_STATE)]]
    t = []
    for f in target_seq:
        if f == BOUNDRY_STATE:
            t.append([(BOUNDRY_STATE, BOUNDRY_STATE)])
        else:
            tups = [(translations_probs.get((f, e), float('-inf')), f, (i, e)) for i, e in enumerate(source_seq[1:-1])]
            tups.sort(reverse=True)
            
            if len(source_seq) > 5 and len(source_seq) < 15:
                tups = tups[:int(len(tups) * 0.5)]
            elif len(source_seq) >= 15 and len(source_seq) < 35:
                tups = tups[:int(len(tups) * 0.3)]
            elif len(source_seq) >= 35:
                tups = tups[:int(len(tups) * 0.1)]
            
            (ps, fs, ts) = zip(*tups)
            t.append(list(ts))
    #pprint(trelis)
    #pprint(t)
    #pdb.set_trace()
    return t


def get_jump_transition(current_state, prev_state, sent_length):
    jkey = jump_key(current_state, prev_state)
    jiip_key = jump_iip_key(sent_length, prev_state)
    if jkey not in jump_counts or jiip_key not in jump_counts_iip:
        return float('-inf')
    else:
        return jump_counts[jkey] - jump_counts_iip[jiip_key]


def get_transition(current_state, prev_state):
    return alignment_probs.get((current_state, prev_state), float('-inf'))


def get_emission(obs, state):
    return translations_probs.get((obs, state), float('-inf'))


def do_accumilate_posterior_obs(accumilation_dict, obs, aj, ei, posterior_unigram_val):
    # these are actual counts in log space!!
    if isinstance(obs, str) and (not isinstance(aj, tuple)) and isinstance(ei, str):
        if ('count_obs', obs) in accumilation_dict:
            accumilation_dict[('count_obs', obs)] = lu.logadd(accumilation_dict[('count_obs', obs)], posterior_unigram_val)
        else:
            accumilation_dict[('count_obs', obs)] = posterior_unigram_val
        if ('count_state', aj) in accumilation_dict:
            accumilation_dict[('count_state', aj)] = lu.logadd(accumilation_dict[('count_state', aj)], posterior_unigram_val)
        else:
            accumilation_dict[('count_state', aj)] = posterior_unigram_val

        if ('count_emission', obs, ei) in accumilation_dict:
            accumilation_dict[('count_emission', obs, ei)] = lu.logadd(accumilation_dict[('count_emission', obs, ei)],
                                                                       posterior_unigram_val)
        else:
            accumilation_dict[('count_emission', obs, ei)] = posterior_unigram_val
            # doing total counts ...
        if ('any_emission_from', ei) in accumilation_dict:
            accumilation_dict[('any_emission_from', ei)] = lu.logadd(accumilation_dict[('any_emission_from', ei)], posterior_unigram_val)
        else:
            accumilation_dict[('any_emission_from', ei)] = posterior_unigram_val
        return accumilation_dict
    else:
        print 'obs must be string, aj must be str, ei must be string'
        exit()


def do_accumilate_posterior_bigrams_jump(accumilation_dict, aj, aj_1, posterior_bigram_val, sent_length):
    # these are actual counts in log space!!
    if not isinstance(aj, tuple) or isinstance(aj_1, tuple):
        jkey = jump_key(aj, aj_1)
        jiip_key = jump_iip_key(sent_length, aj_1)
        accumilation_dict[jkey] = lu.logadd(accumilation_dict.get(jkey, float('-inf')), posterior_bigram_val)
        accumilation_dict[jiip_key] = accumilation_dict.get(jiip_key, set([]))
        accumilation_dict[jiip_key].add(jkey)
        return accumilation_dict
    else:
        print 'aj and aj_1 should be str ### or int', aj, aj_1
        exit()


def do_accumilate_posterior_bigrams(accumilation_dict, aj, aj_1, posterior_bigram_val):
    # these are actual counts in log space!!
    if not isinstance(aj, tuple) or isinstance(aj_1, tuple):
        if ('count_transition', aj, aj_1) not in accumilation_dict:
            accumilation_dict[('count_transition', aj, aj_1)] = posterior_bigram_val
        else:
            accumilation_dict[('count_transition', aj, aj_1)] = lu.logadd(accumilation_dict[('count_transition', aj, aj_1)],
                                                                          posterior_bigram_val)

        if ('any_transition_from', aj_1) not in accumilation_dict:
            accumilation_dict[('any_transition_from', aj_1)] = posterior_bigram_val
        else:
            accumilation_dict[('any_transition_from', aj_1)] = lu.logadd(accumilation_dict[('any_transition_from', aj_1)],
                                                                         posterior_bigram_val)
        return accumilation_dict
    else:
        print 'aj and aj_1 should be str ### or int', aj, aj_1
        exit()


def do_append_posterior_unigrams(appending_dict, position, state, posterior_unigram_val):
    if position in appending_dict:
        appending_dict[position].append((state, posterior_unigram_val))
    else:
        appending_dict[position] = [(state, posterior_unigram_val)]
    return appending_dict


def flatten_backpointers(bt):
    reverse_bt = []
    while len(bt) > 0:
        x = bt.pop()
        reverse_bt.append(x)
        if len(bt) > 0:
            bt = bt.pop()
    reverse_bt.reverse()
    return reverse_bt


def get_viterbi_and_forward(obs_sequence, trelis, source_len):
    pi = {(0, (BOUNDRY_STATE, BOUNDRY_STATE)): 0.0}
    alpha_pi = {(0, (BOUNDRY_STATE, BOUNDRY_STATE)): 0.0}
    #pi[(0, START_STATE)] = 1.0  # 0,START_STATE
    arg_pi = {(0, (BOUNDRY_STATE, BOUNDRY_STATE)): []}
    for k in range(1, len(obs_sequence)):  # the words are numbered from 1 to n, 0 is special start character
        for v in trelis[k]:  # [1]:
            max_prob_to_bt = {}
            sum_prob_to_bt = []
            target_token = obs_sequence[k]
            source_token = v[1]
            for u in trelis[k - 1]:  # [1]:
                aj = v[0]
                aj_1 = u[0]
                #q = get_transition(aj, aj_1)
                q = get_jump_transition(aj, aj_1, source_len)
                e = get_emission(target_token, source_token)
                #print q, aj, aj_1, source_len, jump_key(aj, aj_1), jump_counts.get(jump_key(aj, aj_1), float('-inf')), jump_counts_iip.get(
                #    jump_iip_key(source_len, aj_1), float('-inf')), jump_iip_key(source_len, aj_1)
                #print k
                #print v, '|', u
                #print aj, '|', aj_1, '=', q
                #print target_token, '|', source_token, '=', e
                p = pi[(k - 1, u)] + q + e
                alpha_p = alpha_pi[(k - 1, u)] + q + e
                if len(arg_pi[(k - 1, u)]) == 0:
                    bt = [u]
                else:
                    bt = [arg_pi[(k - 1, u)], u]
                max_prob_to_bt[p] = bt
                sum_prob_to_bt.append(alpha_p)

            max_bt = max_prob_to_bt[max(max_prob_to_bt)]
            new_pi_key = (k, v)
            pi[new_pi_key] = max(max_prob_to_bt)
            #print 'mu   ', new_pi_key, '=', pi[new_pi_key], exp(pi[new_pi_key])
            alpha_pi[new_pi_key] = lu.logadd_of_list(sum_prob_to_bt)
            #print 'alpha', new_pi_key, '=', alpha_pi[new_pi_key], exp(alpha_pi[new_pi_key])
            arg_pi[new_pi_key] = max_bt

    max_bt = max_prob_to_bt[max(max_prob_to_bt)]
    max_p = max(max_prob_to_bt)
    max_bt = flatten_backpointers(max_bt)
    return max_bt, max_p, alpha_pi


def get_backwards(obs, trelis, alpha_pi, source_len=None):
    n = len(obs) - 1  # index of last word
    beta_pi = {(n, (BOUNDRY_STATE, BOUNDRY_STATE)): 0.0}
    S = alpha_pi[(n, (BOUNDRY_STATE, BOUNDRY_STATE))]  # from line 13 in pseudo code
    p_unigrams = {}
    p_obs = {}
    p_trans = {}
    for k in range(n, 0, -1):
        for v in trelis[k]:
            pb = beta_pi[(k, v)]
            aj = v[0]
            source_token = v[1]
            posterior_unigram_val = beta_pi[(k, v)] + alpha_pi[(k, v)] - S
            p_obs = do_accumilate_posterior_obs(p_obs, obs[k], aj, source_token, posterior_unigram_val)
            p_unigrams = do_append_posterior_unigrams(p_unigrams, k, v, posterior_unigram_val)

            for u in trelis[k - 1]:
                #print 'reverse transition', 'k', k, 'u', u, '->', 'v', v
                aj_1 = u[0]
                q = get_jump_transition(aj, aj_1, source_len)
                #q = get_transition(aj, aj_1)
                target_token = obs[k]
                e = get_emission(target_token, source_token)
                p = q + e
                beta_p = pb + p
                new_pi_key = (k - 1, u)
                if new_pi_key not in beta_pi:  # implements lines 16
                    beta_pi[new_pi_key] = beta_p
                else:
                    beta_pi[new_pi_key] = lu.logadd(beta_pi[new_pi_key], beta_p)
                    #print 'beta     ', new_pi_key, '=', beta_pi[new_pi_key], exp(beta_pi[new_pi_key])
                posterior_bigram_val = alpha_pi[(k - 1, u)] + p + beta_pi[(k, v)] - S
                #posterior_bigram_val = "%.3f" % (exp(alpha_pi[(k - 1, u)] + p + beta_pi[(k, v)] - S))
                #p_trans = do_accumilate_posterior_bigrams(p_trans, aj, aj_1, posterior_bigram_val)
                p_trans = do_accumilate_posterior_bigrams_jump(p_trans, aj, aj_1, posterior_bigram_val, source_len)
    return p_unigrams, p_trans, p_obs, S, beta_pi


def format_alignments(init_aligns):
    aligns = {}
    for line in init_aligns:
        [snum, inum, jnum] = line.split()
        s_index = int(snum)
        a = aligns.get(snum, [BOUNDRY_STATE])
        j_index = int(jnum)
        i_index = int(inum)
        if len(a) - 1 < j_index:
            pad = [0] * (j_index - (len(a) - 1))
            a = a + pad
            a[j_index] = i_index
        else:
            a[j_index] = i_index
        aligns[snum] = a
    A = []
    for k, a in sorted(aligns.iteritems()):
        A.append(a + [BOUNDRY_STATE])
    return A


def get_jump_mle(alignments_split, source_split):
    jump_counts = {}
    jump_keys_by_sentence_len = {}
    for a, s in zip(alignments_split, source_split):
        alignment_bigrams = [(a[i], a[i - 1]) for i in range(1, len(a))]
        for j1, j0 in alignment_bigrams:
            jkey = jump_key(j1, j0)
            jiip_key = jump_iip_key(len(s), j0)
            jump_counts[jkey] = lu.logadd(jump_counts.get(jkey, float('-inf')), 0.0)
            jump_keys_by_sentence_len[jiip_key] = jump_keys_by_sentence_len.get(jiip_key, set([]))
            jump_keys_by_sentence_len[jiip_key].add(jkey)
    jiip = {}
    for jiip_key, jkeys in jump_keys_by_sentence_len.iteritems():
        jkeys_val = [jump_counts[jk] for jk in jkeys]
        jiip[jiip_key] = lu.logadd_of_list(jkeys_val)
    return jump_counts, jiip


def get_alignment_mle(init_alignments):
    alignment_mles = {}
    alignments_count = {}
    for a in init_alignments:
        num_alignment_states = len(a) - 2
        alignment_bigrams = [(a[i], a[i - 1]) for i in range(1, len(a))]

        for ab in alignment_bigrams:
            alignments_count[ab[1]] = alignments_count.get(ab[1], 0) + 1
            alignments_count[ab] = alignments_count.get(ab, 0) + 1

    for ap in alignments_count:
        if isinstance(ap, tuple):
            (aj, aj_1) = ap
            alignment_mles[ap] = log(float(alignments_count[ap]) / float(alignments_count[aj_1]))
        else:
            alignment_mles[ap] = log(float(alignments_count[ap]) / float(num_alignment_states))

    return alignment_mles


def update_jump_alignment_mle(posterior_alignment_counts):
    for key in posterior_alignment_counts:
        if key[0] == jump_key(0, 0)[0]:
            jump_counts[key] = posterior_alignment_counts[key]

    for key in posterior_alignment_counts:
        if key[0] == jump_iip_key(0, 0)[0]:
            jkey_vals = [jump_counts[jk] for jk in posterior_alignment_counts[key]]
            log_add_jkey_vals = lu.logadd_of_list(jkey_vals)
            jump_counts_iip[key] = log_add_jkey_vals


def update_alignment_mle(posterior_alignment_counts):
    for k in posterior_alignment_counts:
        if k[0] == 'count_transition':
            (comment, aj, aj_1) = k
            count_ajaj_1 = posterior_alignment_counts['count_transition', aj, aj_1]
            count_aj_1 = posterior_alignment_counts['any_transition_from', aj_1]
            if count_aj_1 == float('-inf'):
                alignment_probs[aj, aj_1] = float('-inf')
            else:
                alignment_probs[aj, aj_1] = count_ajaj_1 - count_aj_1
        else:
            (comment, aj_1) = k
            count_aj_1 = posterior_alignment_counts['any_transition_from', aj_1]
            alignment_probs[aj_1] = count_aj_1


def update_translation_mle(posterior_emission_counts):
    for f, e in translations_probs:
        try:
            counts_fe = posterior_emission_counts['count_emission', f, e]
            count_e = posterior_emission_counts['any_emission_from', e]
            if count_e == float('-inf'):
                translations_probs[f, e] = float('-inf')
            else:
                translations_probs[f, e] = counts_fe - count_e
        except KeyError:
            pass


def get_translation_mle(init_trans):
    translations_mle = {}
    for line in init_trans:
        [fi, ej, p] = line.split()
        if float(p) != 0.0:
            translations_mle[fi, ej] = log(float(p))
        else:
            translations_mle[fi, ej] = float('-inf')
    translations_mle[BOUNDRY_STATE, BOUNDRY_STATE] = 0.0
    return translations_mle


def parseargs(args):
    try:
        source_idx = args.index('-s')
        target_idx = args.index('-t')
        source = args[source_idx + 1]
        target = args[target_idx + 1]
        init_translations = args[args.index('-it') + 1]
        init_align = args[args.index('-ia') + 1]
        save_translations_learned = args[args.index('-p') + 1]
        save_alignment_out = args[args.index('-a') + 1]
        source_alignment_test = args[args.index('-as') + 1]
        target_alignment_test = args[args.index('-at') + 1]
        return source, target, init_translations, init_align, save_translations_learned, save_alignment_out, source_alignment_test, target_alignment_test
    except (ValueError, IndexError) as er:
        print 'Usage: python model1.py -t [train target] -s [train source] -ia [init alignments] -it [initial translations]' \
              ' -p [save translations] ' \
              '-a [save alignment test] -as [alignment test source] -at [alignment test target]'
        #return 'dummy.en', 'dummy.es', 'dummy.trans', 'dummy.align', 'hmm.trans', 'hmm.align', 'dummy.en', 'dummy.es'
        #return 'corpus.en', 'corpus.es', 'model1-fwd-out.trans', 'model1-fwd-out.align', 'hmm.trans', 'hmm.align', 'dev.en', 'dev.es'
        exit()


def accumilate(accumilator, addition):
    for k, val in addition.iteritems():
        if isinstance(val, float):
            s = lu.logadd(accumilator.get(k, float('-inf')), val)
            accumilator[k] = s
        elif isinstance(val, set):
            accumilator[k] = accumilator.get(k, set([]))
            accumilator[k].update(val)
    return accumilator


if __name__ == "__main__":

    source, target, init_translations, init_alignments, save_translations_learned, save_alignment_out, source_alignment_test, target_alignment_test = parseargs(
        sys.argv)

    corpus_source = open(source, 'r').readlines()
    corpus_target = open(target, 'r').readlines()
    z = [(s, t) for s, t in zip(corpus_source, corpus_target) if (s.strip() != '' and t.strip() != '')]
    cs, ct = zip(*z)
    corpus_source = list(cs)
    corpus_target = list(ct)
    alignment_split = format_alignments(open(init_alignments, 'r').readlines())
    init_translations = open(init_translations, 'r').readlines()

    source_split = [[BOUNDRY_STATE, 'NULL'] + i.split() + [BOUNDRY_STATE] for i in corpus_source]

    target_split = [[BOUNDRY_STATE] + i.split() + [BOUNDRY_STATE] for i in corpus_target]
    jump_counts, jump_counts_iip = get_jump_mle(alignment_split, source_split)
    alignment_probs = get_alignment_mle(alignment_split)
    translations_probs = get_translation_mle(init_translations)

    trelis_split = []
    for e, f in zip(source_split, target_split):
        t = get_trelis(f, e)
        trelis_split.append(t)

    for i in range(10):
        accu_alpha = 0.0
        accu_mu = 0.0
        posterior_transitions_accumilation = {}
        posterior_emission_accumilation = {}
        final_alignments = []
        for idx, (e, f, t) in enumerate(zip(source_split, target_split, trelis_split)):
            max_bt, max_p, alpha_pi = get_viterbi_and_forward(f, t, len(e))
            #print 'max_p', max_p
            posterior_uni, posterior_trans, posterior_emission, S, beta_pi = get_backwards(f, t, alpha_pi, len(e))
            accu_alpha += S
            accu_mu += max_p
            posterior_transitions_accumilation = accumilate(posterior_transitions_accumilation, posterior_trans)
            posterior_emission_accumilation = accumilate(posterior_emission_accumilation, posterior_emission)
            [out_alignments, out_emissions] = zip(*max_bt)
            final_alignments = final_alignments + list(out_alignments)

        update_translation_mle(posterior_emission_accumilation)
        #update_alignment_mle(posterior_transitions_accumilation)
        update_jump_alignment_mle(posterior_transitions_accumilation)
        print 'iteration', i, 'mu', accu_mu, 'alpha', accu_alpha
        #pdb.set_trace()
        writer = open(save_alignment_out + '-' + str(i), 'w')
        ia = 0
        for aj in final_alignments:
            if aj == '###':
                ia += 1
                w = 1
            else:
                if aj != 0:
                    writer.write(str(ia) + ' ' + str(aj) + ' ' + str(w) + '\n')
                w += 1
        writer.flush()
        writer.close()







