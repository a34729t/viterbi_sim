# Simulated Annealing simulation to find optimal viterbi parameters

import sys, subprocess, math, random
from math import log, exp
# TODO:
# - Figure out similarity score for cost(test, gold)
# - Set up simulated annealing simulation for probabilities

MIN_P = 0.001
MAX_P = 0.999

def main():
    viterbi_c_path = sys.argv[1]
    data_path = sys.argv[2]
    
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
    iter_max = 100
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
        total_cost = run_viterbi(viterbi_c_path, data, probs)
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
    run_viterbi(viterbi_c_path, data, probs, True)
    
def run_viterbi(viterbi_c_path, data, probs, debug=False):
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
        total_cost += levenshtein(output, gold)
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

if __name__ == "__main__": main()