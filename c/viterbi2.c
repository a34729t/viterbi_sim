// cleaned up viterbi, just log version

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>

typedef enum { false, true } bool;
typedef double (*prob_function_def)(char, char, bool);

int n_states = 17;
int n_observations = 17;
char states[17] = 
    { 'X', '0', '1', '2', '3', '4', '5', '6', '7', 
    '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};
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
    printf("\nviterbi...\n");
    
    // Variables for loops and stuff
    int i, j, t, max_state_index;
    char state_i, state_j, sample_t, max_state;
    double trans_p, max_p;
    
    // Data structures
    double **viterbi_table = malloc(sizeof * viterbi_table * n_states);
    int **best_path_table = malloc(sizeof * best_path_table * n_states);
    
    // Allocate data structures
    if (viterbi_table)
    {
        for (i = 0; i < n_states; i++)
            viterbi_table[i] = malloc(sizeof * viterbi_table[i] * n_states);
    }
    if (best_path_table)
    {
        for (i = 0; i < n_states; i++)
            best_path_table[i] = malloc(sizeof * best_path_table[i] * n_states);
    }
    
    
    // Initialize first column of viterbi table
    if (debug) printf("\nIntialization:\n");
    for (i = 0; i < n_states; i++)
    {
        state_i = states[i];
        sample_t = samples[0];
        viterbi_table[i][0]
            = prior(state_i, 0, true) + emission(sample_t, state_i, true);
        if (debug)
        {
            printf("\t");
            printf("log(prior[%c]) + log(emission[%c][%c]) = %e\n", 
                state_i, sample_t, state_i, viterbi_table[i][0]);
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
            
            if (debug) printf("\ti=%c\n", state_i);
            
            for (j = 0; j < n_states; j++)
            {
                state_j = states[j];
                trans_p = viterbi_table[j][t - 1]
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
                    printf("viterbi_table[%d][%d] + ", j, t-1);
                    printf("transition[%c][%c] + ", state_i, state_j);
                    printf("emission[%c][%c] = ", sample_t, state_j);
                    printf("%e + ", viterbi_table[j][t - 1]);
                    printf("%e + ", transition(state_i, state_j, true));
                    printf("%e = ", emission(sample_t, state_j, true));
                    printf("%e\n", trans_p);
                }
            }
            
            viterbi_table[i][t] = max_p;
            best_path_table[i][t] = max_state_index;
            printf("\tviterbi_table[%d][%d] = %e\n", i, t, viterbi_table[i][t]);
            printf("\tbest_path_table[%d][%d] = %d\n", i, t, best_path_table[i][t]);
            
        }
    }
    
    // Backward: Recover the most probable path
    if (debug) printf("\nBackward:\n");
    
    // Get the final state, according to max probability in the viterbi table
    max_state_index = 0;
    max_p = -DBL_MAX;
    for (i = 0; i < n_states; i++)
    {
        if (viterbi_table[i][n_samples-1] > max_p)
        {
            max_p = viterbi_table[i][n_samples-1];
            max_state_index = i;
        }
    }
    best_path[n_samples - 1] = states[max_state_index];
    
    printf("states[max_state_index]=%c", states[max_state_index]);
    
    // Walk backwards and find the best state path
    for (t = n_samples - 1; t > 0; t--)
    {
        best_path[t-1] = states[best_path_table[max_state_index][t]];
        max_state_index = best_path_table[max_state_index][t];
    }
    
    if (debug)
    {
        for (t = 0; t < n_samples; t++)
        {
            printf("best_path[%d] = %c\n", t, best_path[t]);
        } 
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
        else if (state_i == 0 || state_j == 0)
            return -2.3025850929940455; // 0.1
            // return -0.916290731874155; // 0.4
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
        if (state == observation)
            return -0.2231435513142097;
        else
            return -1.6094379124341003;      
    }
}


void dtmf_example ()
{
    // Basic vars
    // NOTE: The null state (no tone) is X
    
    bool debug = true;
    
    // int n_samples = 6;
    // char samples[6] = { 'X', 'X', 'X', '1', '0', '0'};
    char *samples = "XXXXXXXXXXXXXX6X66XXXXXXXXXXXXXXXXXXXXX";
    int n_samples = strlen(samples);
    char *best_path = malloc(sizeof(int) * n_samples); // this is our output
    
    // viterbi with summed log(prob) to avoid underflows
    viterbi_log(samples, n_samples, best_path, 
        &prior_prob, &transition_prob, &emission_prob, debug);

    printf("\nbest path = { ");
    int t;
    for (t = 0; t < n_samples; t++)
        printf("%d ", best_path[t]);
    printf("}\n");
}

int main(int argc, char* argv)
{
    dtmf_example();
    return 0;
}