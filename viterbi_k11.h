#ifndef VITERBI_K11_H
#define VITERBI_K11_H

/**
 * @brief Codifica un blocco di dati usando un codice convoluzionale K=11, R=1/2.
 *
 * @param input_data Puntatore ai dati di input (array di byte).
 * @param nbits Numero di bit da codificare.
 * @param output_symbols Puntatore al buffer di output per i simboli codificati (deve essere preallocato).
 *                       La dimensione richiesta è (nbits + 10) * 2 byte.
 */
void encode_k11(const unsigned char* input_data, int nbits, unsigned char* output_symbols);

/**
 * @brief Decodifica un blocco di simboli usando l'algoritmo di Viterbi per K=11, R=1/2.
 *        Questa funzione gestisce internamente l'allocazione e la deallocazione della memoria.
 *
 * @param input_symbols Puntatore ai simboli ricevuti (0-255, dove 0 è la massima confidenza per il bit 1, 255 per il bit 0).
 * @param nbits Numero di bit di dati *originali* da decodificare.
 * @param output_data Puntatore al buffer di output per i dati decodificati (deve essere preallocato).
 *                    La dimensione richiesta è (nbits + 7) / 8 byte.
 * @return 0 in caso di successo, -1 in caso di errore.
 */
int decode_k11(const unsigned char* input_symbols, int nbits, unsigned char* output_data);

#endif // VITERBI_K11_H