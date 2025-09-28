#include "convolutional_encoder.h"
#include <iostream>
#include <string>

// Funzione di supporto per stampare un array di bit
void print_array(const std::string& name, const int* arr, int size, int max_elements = 40)
{
    std::cout << name << " (len " << size << "): ";
    for (int i = 0; i < size && i < max_elements; ++i)
    {
        std::cout << arr[i];
    }
    if (size > max_elements)
    {
        std::cout << "...";
    }
    std::cout << std::endl;
}

// Funzione di supporto per stampare un array di float
void print_soft_array(const std::string& name, const float* arr, int size, int max_elements = 40)
{
    std::cout << name << " (len " << size << "): ";
    for (int i = 0; i < size && i < max_elements; ++i)
    {
        std::cout.precision(2);
        std::cout << std::fixed << arr[i] << " ";
    }
    if (size > max_elements)
    {
        std::cout << "...";
    }
    std::cout << std::endl;
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

    // Introduciamo un errore "soft": un valore viene corrotto ma non invertito completamente
    int error_position = 5;
    if (error_position < encoded_size)
    {
        std::cout << "\n-> Errore soft introdotto al simbolo " << error_position << std::endl;
        // Corrompiamo il valore, rendendolo ambiguo (es. da -1.0 a +0.2)
        received_soft_bits[error_position] = 0.2f;
    }
    print_soft_array("Simboli soft ricevuti: ", received_soft_bits, encoded_size);

    // --- 3. Decodifica Soft ---
    int decoded_bits[MAX_INPUT_BITS];
    decode_soft(received_soft_bits, message_size, decoded_bits);
    print_array("Messaggio decodificato:", decoded_bits, message_size);

    // --- 4. Verifica ---
    std::cout << "\n----------------------------------------" << std::endl;
    if (compare_arrays(message_bits, decoded_bits, message_size))
    {
        std::cout << "RISULTATO: Successo!" << std::endl;
        std::cout << "Il messaggio e' stato decodificato correttamente." << std::endl;
    }
    else
    {
        std::cout << "RISULTATO: Fallimento." << std::endl;
        std::cout << "Il messaggio decodificato non corrisponde all'originale." << std::endl;
    }
    std::cout << "----------------------------------------" << std::endl;

    return 0;
}