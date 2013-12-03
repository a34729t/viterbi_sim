// viterbi with clean up, 1D array, ability to change probabilities via cmd line.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

typedef enum { false, true } bool;
typedef double (*prob_function_def)(char, char);

int n_states = 17;
// int n_states = 2;
int n_observations = 17;
char states[17] = 
    { 'X', '0', '1', '2', '3', '4', '5', '6', '7', 
    '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
// char states[2] = 
//     { 'X', '6'};
char observations[17] = 
    { 'X', '0', '1', '2', '3', '4', '5', '6', '7', 
    '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

void viterbi_log (  
                char *samples,
                int n_samples,
                char *best_path, // same length as samples
                prob_function_def prior,
                prob_function_def transition,
                prob_function_def emission,
                bool debug
             )
{
    if (debug) printf("\nviterbi...\n");
    
    // Variables for loops and stuff
    int i, j, t, max_state_index;
    char state_i, state_j, sample_t;
    double trans_p, max_p;
    
    // Data structures
    double *viterbi_table = malloc(sizeof(double) * n_samples * n_states);
    int *best_path_table = malloc(sizeof(int) * n_samples * n_states);
    
    
    // Initialize first column of viterbi table
    if (debug) printf("\nIntialization:\n");
    for (i = 0; i < n_states; i++)
    {
        state_i = states[i];
        sample_t = samples[0];
        viterbi_table[i]
            = prior(state_i, 0) + emission(sample_t, state_i);
        if (debug)
        {
            printf("\t");
            printf("viterbi_table[%d] = log(prior[%c]) + log(emission[%c][%c]) = %e\n", 
                i, state_i, sample_t, state_i, viterbi_table[i]);
        }
    }
    
    // Forward: Find most probable state transitions given data samples
    if (debug) printf("\nForward:\n");
    for (t = 1; t < n_samples; t++)
    {
        sample_t = samples[t];
        if (debug) printf("t=%d => sample=%c\n", t, sample_t);
        for (i = 0; i < n_states; i++)
        {
            state_i = states[i];
            max_state_index = 0;
            max_p = -DBL_MAX;
            
            for (j = 0; j < n_states; j++)
            {
                state_j = states[j];
                trans_p = viterbi_table[((t - 1) * n_states) + j]
                    + transition(state_i, state_j)
                    + emission(sample_t, state_j);
                
                if (trans_p > max_p)
                {
                    max_state_index = j;
                    max_p = trans_p;
                }
                
                if (debug)
                {
                    printf("\t\t");
                    printf("viterbi_table[%d] + ", ((t - 1) * n_states) + j);
                    printf("transition[%c][%c] + ", state_i, state_j);
                    printf("emission[%c][%c] = ", sample_t, state_j);
                    printf("%e + ", viterbi_table[((t - 1) * n_states) + j]);
                    printf("%e + ", transition(state_i, state_j));
                    printf("%e = ", emission(sample_t, state_j));
                    printf("%e\n", trans_p);
                }
            }
            
            viterbi_table[t * n_states + i] = max_p;
            best_path_table[t * n_states + i] = max_state_index;
            if (debug)
            {
                printf("\tviterbi_table[%d] = %d => %e\n", t * n_states + i, best_path_table[t * n_states + i], max_p);
                printf("\tbest_path_table[%d] = %d => %d\n", t * n_states + i, best_path_table[t * n_states + i], max_state_index);
            }
            
        }
        
        // print out rows of viterbi and best path table
        if (debug)
        {
            printf("viterbi, [ ");
            for (i = 0; i < n_states; i++)
            {
                printf("%e ", viterbi_table[t * n_states + i]);
            }
            printf("]\n");
            printf("best_path, [ ");
            for (i = 0; i < n_states; i++)
            {
                printf("[%d], %d ", t * n_states + i, best_path_table[t * n_states + i]);
            }
            printf("]\n");
        }
    }
    
    // Backward: Recover the most probable path
    if (debug) printf("\nBackward:\n");
    
    // Find the most likely ending state
    max_state_index = 0;
    max_p = -DBL_MAX;
    for (i = 0; i < n_states; i++)
    {
        // printf("i=%d, v=%e\n", i, viterbi_table[(n_samples - 1) * n_states + i]);
        if (viterbi_table[(n_samples - 1) * n_states + i] > max_p)
        {
            max_p = viterbi_table[(n_samples-1) * n_states + i];
            max_state_index = i;
        }
    }
    best_path[n_samples - 1] = states[max_state_index];

    // Trace the most probable path backwards
    for (t = n_samples - 1; t > 0; t--)
    {
        // printf("best_path_table[%d * n_states + %d]=%d\n", 
        //     t, max_state_index, best_path_table[t * n_states + max_state_index]);
        best_path[t-1] = states[best_path_table[t * n_states + max_state_index]];
        max_state_index = best_path_table[t * n_states + max_state_index];
    }

    // Dealloc the data structures
    free(viterbi_table);
    free(best_path_table);
}

double prior_prob (char state_i, char state_j) {
    // We ignore state_j, we just want to use the same typedef as the other 
    // prob functions
    return -1.3862943611198906;
}

double transition_null_to_null;
double transition_null_to_signal;
double transition_signal_to_null;
double transition_signal_to_signal;
double transition_signal_to_other_signal;

double transition_prob (char state_i, char state_j) {
    if (state_i == 0 && state_j == 0)
        return transition_null_to_null;
    else if (state_i == 0 && state_j != 0)
        return transition_null_to_signal;
    else if (state_i != 0 && state_j == 0)
        return transition_signal_to_null;
    else if (state_i != 0 && state_i == state_j)
        return transition_signal_to_signal;
    else if (state_i != 0 && state_i != state_j)
        return transition_signal_to_other_signal;
}

double emission_null_to_null;
double emission_null_to_signal;
double emission_signal_to_null;
double emission_signal_to_signal;
double emission_signal_to_other_signal;

double emission_prob (char observation, char state) {
    // The main idea here is that it is very unlikely to observe a non-null when
    // we are in the null state, and vice versa.
    if (observation == 'X' && state == observation)
        return emission_null_to_null;
    if (observation == 'X' && state != 'X')
        return emission_null_to_signal;
    if (observation != 'X' && state == 'X')
        return emission_signal_to_null;
    if (observation != 'X' && state == observation)
        return emission_signal_to_signal;
    if (observation != 'X' && observation != 'X')
        return emission_signal_to_other_signal;
}

int main(int argc, char* argv[])
{
    bool debug = false;
    char* samples = argv[1];
    int n_samples = strlen(samples);
    char *best_path = malloc(sizeof(int) * n_samples); // output
    
    // Set up emission probabilities
    sscanf(argv[2], "%lf", &emission_null_to_null);
    sscanf(argv[3], "%lf", &emission_null_to_signal);
    sscanf(argv[4], "%lf", &emission_signal_to_null);
    sscanf(argv[5], "%lf", &emission_signal_to_signal);
    sscanf(argv[6], "%lf", &emission_signal_to_other_signal);
    
    // Set up transmission probabilities
    sscanf(argv[7], "%lf", &transition_null_to_null);
    sscanf(argv[8], "%lf", &transition_null_to_signal);
    sscanf(argv[9], "%lf", &transition_signal_to_null);
    sscanf(argv[10], "%lf", &transition_signal_to_signal);
    sscanf(argv[11], "%lf", &transition_signal_to_other_signal);
    
    // viterbi with summed log(prob) to avoid underflows
    viterbi_log(samples, n_samples, best_path, 
        &prior_prob, &transition_prob, &emission_prob, debug);
    
    printf("%s", best_path);
    return 0;
}