#include "viterbi_k11.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h> // Per INT_MAX

// --- Costanti e Definizioni Interne ---
#define K 11
#define NUM_STATES (1 << (K - 1)) // 1024
#define POLYA 03345
#define POLYB 03613

// Funzione helper per contare i bit (popcount)
static int popcount(int n) {
    int count = 0;
    while (n > 0) {
        n &= (n - 1);
        count++;
    }
    return count;
}


// --- Implementazione delle Funzioni Pubbliche ---

void encode_k11(const unsigned char* input_data, int nbits, unsigned char* output_symbols) {
    unsigned int sr = 0;
    const int flush_bits = K - 1;

    for (int i = 0; i < nbits + flush_bits; i++) {
        int bit = (i < nbits) ? ((input_data[i >> 3] >> (7 - (i & 7))) & 1) : 0;

        // **CORREZIONE CRITICA:** La logica del registro ora è coerente con il decodificatore.
        // Il bit di input viene combinato con lo stato *prima* dello shift.
        unsigned int temp_reg = (bit << (K-1)) | sr;

        // Calcola i bit di output
        output_symbols[2 * i] = popcount(temp_reg & POLYA) % 2;
        output_symbols[2 * i + 1] = popcount(temp_reg & POLYB) % 2;

        // Aggiorna lo stato del registro (shift a destra, nuovo bit entra a sinistra)
        sr = (sr >> 1) | (bit << (K - 2));
    }
}

int decode_k11(const unsigned char* input_symbols, int nbits, unsigned char* output_data) {
    if (input_symbols == NULL || output_data == NULL || nbits <= 0) return -1;

    const int flush_bits = K - 1;
    const int total_bits = nbits + flush_bits;

    // --- Allocazione Memoria ---
    int* path_metrics = (int*)malloc(NUM_STATES * sizeof(int));
    int* new_path_metrics = (int*)malloc(NUM_STATES * sizeof(int));
    int* traceback_path = (int*)malloc(total_bits * NUM_STATES * sizeof(int));

    if (path_metrics == NULL || new_path_metrics == NULL || traceback_path == NULL) {
        free(path_metrics);
        free(new_path_metrics);
        free(traceback_path);
        return -1;
    }

    // --- Inizializzazione ---
    for (int i = 0; i < NUM_STATES; i++) {
        path_metrics[i] = INT_MAX;
    }
    path_metrics[0] = 0;

    // --- Add-Compare-Select ---
    for (int t = 0; t < total_bits; t++) {
        for (int i = 0; i < NUM_STATES; i++) {
            new_path_metrics[i] = INT_MAX;
        }

        const int received_sym0 = input_symbols[t * 2];
        const int received_sym1 = input_symbols[t * 2 + 1];

        for (int prev_state = 0; prev_state < NUM_STATES; prev_state++) {
            if (path_metrics[prev_state] == INT_MAX) continue;

            for (int input_bit = 0; input_bit < 2; input_bit++) {
                int next_state = (prev_state >> 1) | (input_bit << (K - 2));

                unsigned int temp_reg = (input_bit << (K-1)) | prev_state;
                int expected_bit0 = popcount(temp_reg & POLYA) % 2;
                int expected_bit1 = popcount(temp_reg & POLYB) % 2;

                int metric0 = abs(received_sym0 - (expected_bit0 == 0 ? 255 : 0));
                int metric1 = abs(received_sym1 - (expected_bit1 == 0 ? 255 : 0));
                int branch_metric = metric0 + metric1;

                int total_metric = path_metrics[prev_state] + branch_metric;

                if (total_metric < new_path_metrics[next_state]) {
                    new_path_metrics[next_state] = total_metric;
                    traceback_path[t * NUM_STATES + next_state] = prev_state;
                }
            }
        }
        memcpy(path_metrics, new_path_metrics, NUM_STATES * sizeof(int));
    }

    // --- Traceback ---
    int current_state = 0;
    int min_metric = INT_MAX;
    for(int i=0; i < NUM_STATES; ++i){
        if(path_metrics[i] < min_metric){
            min_metric = path_metrics[i];
            current_state = i;
        }
    }

    memset(output_data, 0, (nbits + 7) / 8);

    for (int i = total_bits - 1; i >= 0; i--) {
        int prev_state = traceback_path[i * NUM_STATES + current_state];
        int decoded_bit = (current_state >> (K - 2)) & 1;

        if (i < nbits) {
            if (decoded_bit) {
                output_data[i >> 3] |= (1 << (7 - (i & 7)));
            }
        }
        current_state = prev_state;
    }

    // --- Pulizia ---
    free(path_metrics);
    free(new_path_metrics);
    free(traceback_path);

    return 0;
}