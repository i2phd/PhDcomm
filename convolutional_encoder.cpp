#include "convolutional_encoder.h"
#include <limits> // Per numeric_limits

// Funzione di supporto per contare i bit impostati a 1 (popcount)
int popcount(int n)
{
    int count = 0;
    while (n > 0)
    {
        n &= (n - 1);
        count++;
    }
    return count;
}

// Implementazione della funzione di codifica
void encode(const int* input_bits, int input_size, int* output_bits)
{
    if (input_size > MAX_INPUT_BITS)
    {
        // Se l'input supera il buffer massimo, non fare nulla.
        return;
    }

    int shift_register = 0;
    int output_idx = 0;

    // Codifica i bit di input
    for (int i = 0; i < input_size; ++i)
    {
        int bit = input_bits[i];
        int full_register_state = (bit << M) | shift_register;
        output_bits[output_idx++] = popcount(full_register_state & G1_OCT) % 2;
        output_bits[output_idx++] = popcount(full_register_state & G2_OCT) % 2;
        shift_register >>= 1;
        shift_register |= (bit << (M - 1));
    }

    // Aggiungi M bit di zero per il flushing
    for (int i = 0; i < M; ++i)
    {
        int bit = 0;
        int full_register_state = (bit << M) | shift_register;
        output_bits[output_idx++] = popcount(full_register_state & G1_OCT) % 2;
        output_bits[output_idx++] = popcount(full_register_state & G2_OCT) % 2;
        shift_register >>= 1;
        shift_register |= (bit << (M - 1));
    }
}

// Implementazione della funzione di decodifica SOFT (Algoritmo di Viterbi)
void decode_soft(const float* received_soft_bits, int input_size, int* decoded_bits)
{
    if (input_size > MAX_INPUT_BITS)
    {
        return;
    }

    const int num_states = 1 << M;
    const int num_steps = input_size + M;
    const float neg_infinity = -std::numeric_limits<float>::max();

    // Strutture dati di Viterbi come array statici per evitare lo stack overflow
    static float path_metrics[num_states];
    static float new_path_metrics[num_states];
    static int traceback_path[MAX_DECODED_BUFFER][num_states];
    static float expected_soft_outputs[num_states][2][2];
    static bool is_initialized = false;

    // Inizializza la tabella degli output attesi solo la prima volta
    if (!is_initialized)
    {
        for (int state = 0; state < num_states; ++state)
        {
            for (int input_bit = 0; input_bit < 2; ++input_bit)
            {
                int full_register = (input_bit << M) | state;
                // Mappa bit 0 -> +1.0, bit 1 -> -1.0
                expected_soft_outputs[state][input_bit][0] = 1.0f - 2.0f * (popcount(full_register & G1_OCT) % 2);
                expected_soft_outputs[state][input_bit][1] = 1.0f - 2.0f * (popcount(full_register & G2_OCT) % 2);
            }
        }
        is_initialized = true;
    }

    // Inizializzazione delle metriche di percorso (deve essere fatta a ogni chiamata)
    for (int i = 0; i < num_states; ++i)
    {
        path_metrics[i] = neg_infinity;
    }
    path_metrics[0] = 0.0f;

    // Add-Compare-Select (ACS)
    for (int t = 0; t < num_steps; ++t)
    {
        for (int i = 0; i < num_states; ++i)
        {
            new_path_metrics[i] = neg_infinity;
        }

        const float received_pair[2] = {received_soft_bits[t * 2], received_soft_bits[t * 2 + 1]};

        for (int prev_state = 0; prev_state < num_states; ++prev_state)
        {
            if (path_metrics[prev_state] == neg_infinity)
            {
                continue;
            }

            for (int input_bit = 0; input_bit < 2; ++input_bit)
            {
                int next_state = (prev_state >> 1) | (input_bit << (M - 1));

                float branch_metric = received_pair[0] * expected_soft_outputs[prev_state][input_bit][0] +
                                      received_pair[1] * expected_soft_outputs[prev_state][input_bit][1];

                float total_metric = path_metrics[prev_state] + branch_metric;

                if (total_metric > new_path_metrics[next_state])
                {
                    new_path_metrics[next_state] = total_metric;
                    traceback_path[t][next_state] = prev_state;
                }
            }
        }
        for (int i = 0; i < num_states; ++i)
        {
            path_metrics[i] = new_path_metrics[i];
        }
    }

    // Traceback
    int current_state = 0;

    for (int t = num_steps - 1; t >= 0; --t)
    {
        int prev_state = traceback_path[t][current_state];
        int input_bit = (current_state >> (M - 1)) & 1;

        if (t < input_size)
        {
            decoded_bits[t] = input_bit;
        }
        current_state = prev_state;
    }
}