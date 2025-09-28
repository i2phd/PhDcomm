#ifndef CONVOLUTIONAL_ENCODER_H
#define CONVOLUTIONAL_ENCODER_H

// --- Costanti del Codice Convoluzionale ---
const int K = 11; // Lunghezza di vincolo
const int M = K - 1; // Memoria del codificatore (10)

// Polinomi generatori (ottale)
const int G1_OCT = 03345; // 11011100101 (binario)
const int G2_OCT = 03613; // 11110001011 (binario)

// --- Limiti per Array Statici ---
const int MAX_INPUT_BITS = 1000;
const int MAX_ENCODED_BITS = (MAX_INPUT_BITS + M) * 2;
const int MAX_DECODED_BUFFER = MAX_INPUT_BITS + M;


/**
 * @brief Codifica una sequenza di bit (0/1).
 *
 * @param input_bits Puntatore all'array di bit di input.
 * @param input_size Dimensione dell'array di input (non deve superare MAX_INPUT_BITS).
 * @param output_bits Puntatore all'array di output (dimensione (input_size + M) * 2).
 */
void encode(const int* input_bits, int input_size, int* output_bits);

/**
 * @brief Decodifica una sequenza di simboli soft (es. +1.0, -1.0)
 *        usando l'algoritmo di Viterbi a decisione soft.
 *
 * @param received_soft_bits Puntatore all'array di simboli soft ricevuti.
 * @param input_size La dimensione del messaggio originale (prima della codifica).
 * @param decoded_bits Puntatore all'array di output per i bit decodificati (dimensione input_size).
 */
void decode_soft(const float* received_soft_bits, int input_size, int* decoded_bits);

#endif // CONVOLUTIONAL_ENCODER_H