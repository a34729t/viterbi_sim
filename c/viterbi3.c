// cleaned up viterbi, just log version
// use a 1D array for viterbi table, with mapping:
//      i, j => i * n + j (m rows, n cols)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

typedef enum { false, true } bool;
typedef double (*prob_function_def)(char, char, bool);

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
            = prior(state_i, 0, true) + emission(sample_t, state_i, true);
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
                    + transition(state_i, state_j, true)
                    + emission(sample_t, state_j, true);
                
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
                    printf("%e + ", transition(state_i, state_j, true));
                    printf("%e = ", emission(sample_t, state_j, true));
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

double prior_prob (char state_i, char state_j, bool log_prob) {
    // We ignore state_j, we just want to use the same typedef as the other 
    // prob functions
    if (!log_prob)
        return 0.25;
    else
        return -1.3862943611198906;
}


double transition_prob (char state_i, char state_j, bool log_prob) {
    if (!log_prob)
    {
        if (state_i == 0 && state_j == 0)
            return 0.9; // 0.6
        else if (state_i == 0 || state_j == 0)
            return 0.1; // 0.4
        else if (state_i == state_j)
            return 0.9;
        else
            return 0.1;
    }
    else
    {
        if (state_i == 0 && state_j == 0)
            return -0.10536051565782628; // 0.9
            // return -0.5108256237659907; // 0.6
            // return -0.30102999566; // 0.5
        else if (state_i == 0 || state_j == 0)
            return -2.3025850929940455; // 0.1
            // return -0.916290731874155; // 0.4
            // return -0.5108256237659907; // 0.6
            // return -0.30102999566; // 0.5
        else if (state_i == state_j)
            return -0.10536051565782628;
        else
            return -2.3025850929940455;      
    }
    
}

double emission_prob (char observation, char state, bool log_prob) {
    if (!log_prob)
    {
        if (state == observation)
            return 0.8;
        else
            return 0.2;
    }
    else
    {
        // The main idea here is that it is very unlikely to observe a non-null when
        // we are in the null state, and vice versa.
        if (observation == 'X' && state != 'X')
            return -2.3025850929940455; // 0.1
        if (observation != 'X' && state == 'X')
            return -2.3025850929940455; // 0.1
        if (state == observation)
            return -0.096910013; // 0.8
        else
            return -1.6094379124341003; // 0.2
    }
}


void dtmf_example ()
{
    // Basic vars
    // NOTE: The null state (no tone) is X
    
    // TODO:
    // 1) The below is wrong, should be 0-16
    //      t=1, [ 0 1 2 3 4 5 6 7 8 9 10 11 12 0 1 2 3 ]
    // meh maybe not, seems to work in some cases
    
    bool debug = true;
    // char *samples = "66XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX66"; // works
    // char *samples = "XXXXXXXXXXXXXXXXXXX66XX66XXXXXXXXXXXXXXXXXX"; // works
    // char *samples = "XXXXXXXXXXXXXXXXXXX6666XXXXXXXXXXXXXXXXXX"; // works
    // char *samples = "XX66XX"; // works
    char *samples = "XXXXXXXXXXXXXX6X66XXXXXXXXXXXXXXXXXXXXX"; // works
    int n_samples = strlen(samples);
    char *best_path = malloc(sizeof(int) * n_samples); // this is our output
    
    // viterbi with summed log(prob) to avoid underflows
    viterbi_log(samples, n_samples, best_path, 
        &prior_prob, &transition_prob, &emission_prob, debug);

    printf("observed  = %s\n", samples);
    printf("corrected = ");
    int t;
    for (t = 0; t < n_samples; t++)
        printf("%c", best_path[t]);
    printf("\n");
}

int main()
{
    dtmf_example();
    return 0;
}