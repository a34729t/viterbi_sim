# Simulated Annealing simulation to find optimal viterbi parameters

import sys, subprocess, math, random
from math import log, exp
from operator import itemgetter

MIN_P = 0.001
MAX_P = 0.999

# TODO: Figure out how to do N-fold cross-validation
# TODO: More data
# TODO: Figure out why max_null_cost_metric turns other chars into '0'. Doesn't actually affect params though!

def main():
    viterbi_c_path = sys.argv[1]
    data_path = sys.argv[2]
    # Two different cost functions
    cost_function = max_null_cost_metric
    # cost_function = levenshtein
    
    # Load the data file
    # format: received data, true transmitted data (at origin)
    data = [
            ('XXXXXXXAXAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
             'XXXXXXXAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),
            ('XXXXXXXXXAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
             'XXXXXXXXXAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),
            ('XXXXXXXXXXXXAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
             'XXXXXXXXXAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),
            ('XXXXXXXAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXBXXXBXXXXXXXXXXX',
             'XXXXXXXAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXBBBBBBXXXXXXXXX'),
            ('XXXXXXXXXXXXXXXXCXXCCXXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
             'XXXXXXXXXXXXXXXXCCCCCCCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),             
            ('XXXXXXXXXXAXXXXAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
            'XXXXXXXXXXAAAAAAXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),
            ('XXXXXXXXXXXXXXXXXXXXXXXXXBXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
            'XXXXXXXXXXXXXXXXXXXXXXXXBBBXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'),
            ('XXXXXXCXCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAXXXXXXXXXXXX',
            'XXXXXXXCCCCCXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAXXXXXXXXXXXXX'),
            # null case
            ('XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX',
             'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
            ]
    
    # Probabilities
    probs = {
        # emission
        'emission_null_to_null': 0.75,
        'emission_null_to_signal': 0.30,
        'emission_signal_to_null': 0.05,
        'emission_signal_to_signal': 0.95,
        'emission_signal_to_other_signal': 0.10,
        # transition
        'transition_null_to_null': 0.8,
        'transition_null_to_signal': 0.2,
        'transition_signal_to_null': 0.05,
        'transition_signal_to_signal': 0.90,
        'transition_signal_to_other_signal': 0.1
    }
    
    # Simulated Annealing
   
    # Do an initial run
    total_cost_old = 100000000000000
    temp = 1.0
    temp_decrement = 0.01
    iter_max = 200
    iters = 0
    while True:
        print "Temp =", temp
        iters += 1
        
        # Generate a random increase or decrease for prob of a key
        key = random.choice(probs.keys())
        old_p = probs[key]
        while True:
            new_p = round(random.uniform(MIN_P, MAX_P), 3)
            if abs(new_p - old_p) < (temp * 0.5):
                probs[key] = new_p
                print "\tchanging", key, "from",old_p,"to", probs[key]
                break
        
        
        # Run against all test data and compute cost
        total_cost = run_viterbi(viterbi_c_path, data, probs, cost_function)
        print "\ttotal_cost =", total_cost, "vs old_cost =", total_cost_old
        
        # Do we keep old_p?
        if total_cost_old < total_cost:
            probs[key] = old_p
        else:
            # Keep new_p
            # Decrease temp
            diff = (total_cost_old - total_cost)
            if diff > 10: diff = 10
            temp -= (temp_decrement * diff)
            total_cost_old = total_cost
        
        # Finish simulation?
        if temp <= 0 or iters > iter_max: break
        
    # Run simulation and print out
    run_viterbi(viterbi_c_path, data, probs, cost_function, True)
    
    # Print probabilities
    print ''
    print '---------------------------------------'
    print '        Optimized Probabilities        '
    print '---------------------------------------'
    sorted_probs = sorted(probs.iteritems(), key=lambda (k,v): itemgetter(1)(k))
    print_tuples(sorted_probs)
    
def run_viterbi(viterbi_c_path, data, probs, cost_function, debug=False):
    total_cost = 0
    for (test, gold) in data:
        # NOTE: We use log-probabilities to avoid underflow
        args = [
                viterbi_c_path,
                test,
                # emission probabilities
                str(log(probs['emission_null_to_null'])),
                str(log(probs['emission_null_to_signal'])),
                str(log(probs['emission_signal_to_null'])),
                str(log(probs['emission_signal_to_signal'])),
                str(log(probs['emission_signal_to_other_signal'])),
                # transition probabilities
                str(log(probs['transition_null_to_null'])),
                str(log(probs['transition_null_to_signal'])),
                str(log(probs['transition_signal_to_null'])),
                str(log(probs['transition_signal_to_signal'])),
                str(log(probs['transition_signal_to_other_signal']))
                ]
        output = subprocess.check_output(args) # the reconstructed data (viterbi)
        total_cost += cost_function(output, gold)
        if debug: print output, gold
    return total_cost

def levenshtein(s1, s2):
    # from http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    if len(s1) < len(s2):
        return levenshtein(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]
    
def max_null_cost_metric(s1, s2):
    # Assign a massive cost to getting a null char instead of true char
    
    # get histogram of chars for each str
    # see http://stackoverflow.com/questions/8952101/best-way-to-count-char-occurences-in-a-string
    s1_chars = dict((c, s1.count(c)) for c in s1)
    s2_chars = dict((c, s2.count(c)) for c in s2)
    
    # compare counts - create a dictionary of the counts of s1 - s2
    common_chars = {};
    for c in s1_chars: common_chars[c] = s1_chars[c]
    for c in s2_chars:
        if c in common_chars:
            common_chars[c] -= s2_chars[c]
        else:
            common_chars[c] = -1 * s2_chars[c]
          
    # Assign more value to null char
    for c in common_chars:
        if c == 'X': common_chars[c] *= 10
    
    print 'sum()=', sum(common_chars.values())
    return abs(sum(common_chars.values()))

def print_tuples(tuples, padding = 3):
    # Assume tuples are all the same length!!!
    # Get max length of each index
    length = len(tuples[0])
    max_length = [0] * length
    for i in range(0, length):
        for t in tuples:
            if len(str(t[i])) > max_length[i]:
                max_length[i] = len(str(t[i]))

    # Print out each element with padding + difference of max_length and length of item
    for t in tuples:
        buf = ''
        for i in range(0, length):
            buf += str(t[i])
            for j in range(0, max_length[i] + padding - len(str(t[i]))):
                buf += ' '
        print buf

if __name__ == "__main__": main()