#include "convolutional_encoder.h"
#include <cstdio> // Per printf
#include <string> // Mantenuto per std::string nei nomi dei vettori

// Funzione di supporto per stampare un array di bit
void print_array(const std::string& name, const int* arr, int size, int max_elements = 40)
{
    printf("%s (len %d): ", name.c_str(), size);
    for (int i = 0; i < size && i < max_elements; ++i)
    {
        printf("%d", arr[i]);
    }
    if (size > max_elements)
    {
        printf("...");
    }
    printf("\n");
}

// Funzione di supporto per stampare un array di float
void print_soft_array(const std::string& name, const float* arr, int size, int max_elements = 40)
{
    printf("%s (len %d): ", name.c_str(), size);
    for (int i = 0; i < size && i < max_elements; ++i)
    {
        printf("%.2f ", arr[i]);
    }
    if (size > max_elements)
    {
        printf("...");
    }
    printf("\n");
}

// Funzione per confrontare due array
bool compare_arrays(const int* arr1, const int* arr2, int size)
{
    for (int i = 0; i < size; ++i)
    {
        if (arr1[i] != arr2[i])
        {
            return false;
        }
    }
    return true;
}

int main()
{
    // --- Messaggio di Esempio ---
    const int message_size = 16;
    int message_bits[message_size] = {1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1};

    print_array("Messaggio originale:   ", message_bits, message_size);

    // --- 1. Codifica ---
    const int encoded_size = (message_size + M) * 2;
    int encoded_bits[MAX_ENCODED_BITS];

    encode(message_bits, message_size, encoded_bits);
    print_array("Messaggio codificato:  ", encoded_bits, encoded_size);

    // --- 2. Simulazione di Canale Soft (BPSK + Rumore) ---
    float received_soft_bits[MAX_ENCODED_BITS];
    // Mappa 0 -> +1.0, 1 -> -1.0
    for (int i = 0; i < encoded_size; ++i)
    {
        received_soft_bits[i] = 1.0f - 2.0f * encoded_bits[i];
    }

    // Introduciamo un numero maggiore di errori
    int error_positions[] = {5, 12, 25, 38};
    int num_errors = sizeof(error_positions) / sizeof(int);
    printf("\n-> Introduzione di %d errori soft nelle posizioni:", num_errors);
    for (int i = 0; i < num_errors; ++i)
    {
        printf(" %d", error_positions[i]);
    }
    printf("\n");

    for (int i = 0; i < num_errors; ++i)
    {
        int pos = error_positions[i];
        if (pos < encoded_size)
        {
            // Corrompiamo il valore, rendendolo ambiguo
            received_soft_bits[pos] = 0.2f;
        }
    }
    print_soft_array("Simboli soft ricevuti: ", received_soft_bits, encoded_size);

    // --- 3. Decodifica Soft ---
    int decoded_bits[MAX_INPUT_BITS];
    decode_soft(received_soft_bits, message_size, decoded_bits);
    print_array("Messaggio decodificato:", decoded_bits, message_size);

    // --- 4. Verifica ---
    printf("\n----------------------------------------\n");
    if (compare_arrays(message_bits, decoded_bits, message_size))
    {
        printf("RISULTATO: Successo!\n");
        printf("Il messaggio e' stato decodificato correttamente nonostante i %d errori.\n", num_errors);
    }
    else
    {
        printf("RISULTATO: Fallimento.\n");
        printf("Il messaggio decodificato non corrisponde all'originale.\n");
    }
    printf("----------------------------------------\n");

    return 0;
}