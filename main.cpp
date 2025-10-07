#include <stdio.h>
#include <string.h>
#include "viterbi_k11.h"

int main() {
    const char* original_str = "TEST";
    const int nbytes = strlen(original_str);
    const int nbits = nbytes * 8;

    unsigned char original_data[nbytes];
    memcpy(original_data, original_str, nbytes);

    unsigned char decoded_data[nbytes];
    const int encoded_buffer_size = (nbits + 10) * 2;
    unsigned char encoded_symbols[encoded_buffer_size];
    unsigned char received_symbols[encoded_buffer_size];

    printf("Messaggio Originale: %s\n", original_str);

    // --- 1. Codifica ---
    encode_k11(original_data, nbits, encoded_symbols);
    printf("-> Messaggio codificato con successo.\n");

    // --- 2. Simulazione Canale ---
    // Simula un canale BPSK hard-decision (0->255, 1->0)
    for (int i = 0; i < encoded_buffer_size; ++i) {
        received_symbols[i] = encoded_symbols[i] ? 0 : 255;
    }
    printf("-> Canale simulato (0->255, 1->0).\n");

    // Introduciamo un errore per testare la correzione
    printf("-> Introduzione di un errore nel simbolo 5.\n");
    received_symbols[5] = 255 - received_symbols[5];

    // --- 3. Decodifica ---
    if (decode_k11(received_symbols, nbits, decoded_data) != 0) {
        printf("Errore durante la decodifica.\n");
        return 1;
    }

    // --- 4. Verifica ---
    // Aggiungi un terminatore nullo alla stringa decodificata per stamparla in sicurezza
    char decoded_str[nbytes + 1];
    memcpy(decoded_str, decoded_data, nbytes);
    decoded_str[nbytes] = '\0';

    printf("\n----------------------------------------\n");
    printf("Dati Originali:   %s\n", original_str);
    printf("Dati Decodificati: %s\n", decoded_str);

    if (memcmp(original_data, decoded_data, nbytes) == 0) {
        printf("\nRISULTATO: Successo!\n");
        printf("Il messaggio e' stato decodificato correttamente.\n");
    } else {
        printf("\nRISULTATO: Fallimento.\n");
        printf("Il messaggio decodificato non corrisponde all'originale.\n");
    }
    printf("----------------------------------------\n");

    return 0;
}