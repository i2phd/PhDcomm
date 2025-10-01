// Reed-Solomon Codec - Block size 44, Data 26, Syndrome 18
// Compatibile con Borland C++ Compiler

#include <stdio.h>
#include <string.h>
#include <conio.h>
                                 //======= da correggere ===================

const int N = 44;           // Block size totale
const int K = 26;           // Data symbols (dati utili)
const int T = 18;           // Error correction symbols (syndrome)

// Tabelle per aritmetica GF(256)
unsigned char gf_exp[512];
unsigned char gf_log[256];
int gf_initialized = 0;

// Inizializza le tabelle di logaritmi ed esponenziali per GF(256)
void init_galois_tables()
{
    if (gf_initialized)
        return;

    int x = 1;

    for (int i = 0; i < 255; i++)
    {
        gf_exp[i] = (unsigned char)x;
        gf_log[x] = (unsigned char)i;
        x <<= 1;
        if (x & 0x100)
            x ^= 0x11D;
    }

    for (int i = 255; i < 512; i++)
    {
        gf_exp[i] = gf_exp[i - 255];
    }

    gf_log[0] = 0;
    gf_initialized = 1;
}

// Moltiplicazione in GF(256)
unsigned char gf_mul(unsigned char a, unsigned char b)
{
    if (a == 0 || b == 0)
        return 0;
    return gf_exp[(gf_log[a] + gf_log[b]) % 255];
}

// Divisione in GF(256)
unsigned char gf_div(unsigned char a, unsigned char b)
{
    if (a == 0)
        return 0;
    if (b == 0)
        return 0;
    return gf_exp[(gf_log[a] + 255 - gf_log[b]) % 255];
}

// Valuta un polinomio in un punto
unsigned char poly_eval(unsigned char* poly, int deg, unsigned char x)
{
    unsigned char y = poly[deg];
    for (int i = deg - 1; i >= 0; i--)
    {
        y = gf_mul(y, x) ^ poly[i];
    }
    return y;
}

// CODIFICA Reed-Solomon
void rs_encode(unsigned char* data, unsigned char* codeword)
{
    unsigned char gen[T + 1];
    unsigned char msg[N];

    // Genera il polinomio generatore: prodotto di (x - alpha^i) per i=0..T-1
    memset(gen, 0, T + 1);
    gen[0] = 1;

    for (int i = 0; i < T; i++)
    {
        gen[T - i] = 1;
        for (int j = T - i; j < T; j++)
        {
            gen[j] = gen[j + 1] ^ gf_mul(gen[j], gf_exp[i]);
        }
        gen[T] = gf_mul(gen[T], gf_exp[i]);
    }

    // Prepara il messaggio (dati shiftati a sinistra)
    memset(msg, 0, N);
    memcpy(msg, data, K);

    // Divisione polinomiale per calcolare il resto
    for (int i = 0; i < K; i++)
    {
        unsigned char coef = msg[i];
        if (coef != 0)
        {
            for (int j = 0; j <= T; j++)
            {
                msg[i + j] ^= gf_mul(gen[j], coef);
            }
        }
    }

    // Codeword = dati + resto
    memcpy(codeword, data, K);
    memcpy(codeword + K, msg + K, T);
}

// Calcola le sindromi
int calculate_syndromes(unsigned char* received, unsigned char* syndromes)
{
    int has_error = 0;

    for (int i = 0; i < T; i++)
    {
        syndromes[i] = 0;
        unsigned char alpha_i = gf_exp[i];
        unsigned char alpha_power = 1;

        for (int j = 0; j < N; j++)
        {
            syndromes[i] ^= gf_mul(received[j], alpha_power);
            alpha_power = gf_mul(alpha_power, alpha_i);
        }

        if (syndromes[i] != 0)
            has_error = 1;
    }

    return has_error;
}

// Algoritmo di Berlekamp-Massey
int berlekamp_massey(unsigned char* syndromes, unsigned char* lambda)
{
    unsigned char C[T + 2];
    unsigned char B[T + 2];
    unsigned char temp[T + 2];

    memset(C, 0, T + 2);
    memset(B, 0, T + 2);
    C[0] = 1;
    B[0] = 1;

    int L = 0;
    int m = 1;
    unsigned char b = 1;

    for (int n = 0; n < T; n++)
    {
        unsigned char d = syndromes[n];
        for (int i = 1; i <= L; i++)
        {
            d ^= gf_mul(C[i], syndromes[n - i]);
        }

        if (d == 0)
        {
            m++;
        }
        else
        {
            memcpy(temp, C, T + 2);
            unsigned char scale = gf_div(d, b);

            for (int i = 0; i < T + 2 - m; i++)
            {
                C[m + i] ^= gf_mul(scale, B[i]);
            }

            if (2 * L <= n)
            {
                L = n + 1 - L;
                memcpy(B, temp, T + 2);
                b = d;
                m = 1;
            }
            else
            {
                m++;
            }
        }
    }

    memcpy(lambda, C, T + 2);
    return L;
}

// Ricerca di Chien - trova le radici del polinomio locatore errori
int chien_search(unsigned char* lambda, int L, unsigned char* error_pos)
{
    int num_errors = 0;

    // Prova ogni possibile posizione
    for (int i = 0; i < N; i++)
    {
        unsigned char alpha_inv = gf_exp[(255 - i) % 255];

        // Valuta lambda(alpha^-i)
        unsigned char sum = lambda[0];
        unsigned char x_power = alpha_inv;

        for (int j = 1; j <= L; j++)
        {
            sum ^= gf_mul(lambda[j], x_power);
            x_power = gf_mul(x_power, alpha_inv);
        }

        // Se lambda(alpha^-i) = 0, allora la posizione i ha un errore
        if (sum == 0)
        {
            error_pos[num_errors++] = i;
        }
    }

    return num_errors;
}

// Calcola il polinomio valutatore errori Omega
void compute_omega(unsigned char* syndromes, unsigned char* lambda, int L, unsigned char* omega)
{
    memset(omega, 0, T + 2);

    // Omega(x) = Syndromes(x) * Lambda(x) mod x^T
    for (int i = 0; i < T; i++)
    {
        for (int j = 0; j <= L && j <= i; j++)
        {
            omega[i] ^= gf_mul(syndromes[i - j], lambda[j]);
        }
    }
}

// Formula di Forney per calcolare i valori degli errori
void forney_algorithm(unsigned char* omega, unsigned char* lambda, int L,
                     unsigned char* error_pos, int num_errors, unsigned char* error_values)
{
    for (int k = 0; k < num_errors; k++)
    {
        unsigned char X_inv = gf_exp[(255 - error_pos[k]) % 255];

        // Valuta Omega(X_inv)
        unsigned char omega_val = poly_eval(omega, T - 1, X_inv);

        // Calcola Lambda'(X_inv) - derivata formale
        unsigned char lambda_prime = 0;
        unsigned char x_power = 1;
        for (int i = 1; i <= L; i += 2)
        {
            x_power = gf_mul(x_power, X_inv);
            if (i < L)
                x_power = gf_mul(x_power, X_inv);
            lambda_prime ^= gf_mul(lambda[i], x_power);
        }

        // e_i = Omega(X_inv) / Lambda'(X_inv)
        if (lambda_prime != 0)
        {
            error_values[k] = gf_div(omega_val, lambda_prime);
        }
        else
        {
            error_values[k] = 0;
        }
    }
}

// DECODIFICA Reed-Solomon
int rs_decode(unsigned char* received, unsigned char* decoded)
{
    unsigned char syndromes[T];
    unsigned char lambda[T + 2];
    unsigned char omega[T + 2];
    unsigned char error_pos[T];
    unsigned char error_values[T];

    // Calcola le sindromi
    if (calculate_syndromes(received, syndromes) == 0)
    {
        memcpy(decoded, received, K);
        return 0;
    }

    // Trova il polinomio locatore errori con Berlekamp-Massey
    int L = berlekamp_massey(syndromes, lambda);

    if (L == 0)
    {
        memcpy(decoded, received, K);
        return 0;
    }

    if (L > T / 2)
    {
        return -1;
    }

    // Trova le posizioni degli errori con Chien search
    int num_errors = chien_search(lambda, L, error_pos);

    if (num_errors != L)
    {
        return -1;
    }

    // Calcola i valori degli errori con Forney
    compute_omega(syndromes, lambda, L, omega);
    forney_algorithm(omega, lambda, L, error_pos, num_errors, error_values);

    // Correggi gli errori
    unsigned char corrected[N];
    memcpy(corrected, received, N);

    for (int i = 0; i < num_errors; i++)
    {
        corrected[error_pos[i]] ^= error_values[i];
    }

    // Verifica che la correzione sia corretta ricalcolando le sindromi
    unsigned char check_syndromes[T];
    calculate_syndromes(corrected, check_syndromes);

    int all_zero = 1;
    for (int i = 0; i < T; i++)
    {
        if (check_syndromes[i] != 0)
        {
            all_zero = 0;
            break;
        }
    }

    if (!all_zero)
    {
        return -1;
    }

    memcpy(decoded, corrected, K);

    return num_errors;
}

// Esempio di utilizzo
int main()
{
    init_galois_tables();

    unsigned char data[K];
    unsigned char codeword[N];
    unsigned char decoded[K];

    printf("=== Reed-Solomon (26,44) - Capacita correzione: %d errori ===\n\n", T/2);

    // Inizializza i dati
    for (int i = 0; i < K; i++)
    {
        data[i] = (unsigned char)(0x41 + i);
    }

    printf("Dati originali (%d bytes):\n", K);
    for (int i = 0; i < K; i++)
    {
        printf("%02X ", data[i]);
        if ((i + 1) % 16 == 0) printf("\n");
    }
    printf("\n\n");

    // Codifica
    rs_encode(data, codeword);

    printf("Codeword (%d bytes):\n", N);
    for (int i = 0; i < N; i++)
    {
        printf("%02X ", codeword[i]);
        if ((i + 1) % 16 == 0) printf("\n");
    }
    printf("\n\n");

    // Test senza errori
    printf("--- Test 1: Nessun errore ---\n");
    int result = rs_decode(codeword, decoded);
    printf("Risultato: %d errori\n\n", result);

    // Test con 1 errore
    printf("--- Test 2: 1 errore (posizione 5) ---\n");
    unsigned char test1[N];
    memcpy(test1, codeword, N);
    test1[5] ^= 0xFF;
    result = rs_decode(test1, decoded);
    if (result >= 0)
    {
        printf("Errori corretti: %d\n", result);
        int ok = (memcmp(decoded, data, K) == 0);
        printf("Verifica: %s\n\n", ok ? "OK" : "ERRORE");
    }
    else
    {
        printf("Decodifica FALLITA\n\n");
    }

    // Test con 3 errori
    printf("--- Test 3: 3 errori (posizioni 5, 15, 25) ---\n");
    unsigned char test2[N];
    memcpy(test2, codeword, N);
    test2[5] ^= 0xFF;
    test2[15] ^= 0xAA;
    test2[25] ^= 0x55;
    result = rs_decode(test2, decoded);
    if (result >= 0)
    {
        printf("Errori corretti: %d\n", result);
        int ok = (memcmp(decoded, data, K) == 0);
        printf("Verifica: %s\n\n", ok ? "OK" : "ERRORE");
    }
    else
    {
        printf("Decodifica FALLITA\n\n");
    }

    // Test con 9 errori (limite teorico)
    printf("--- Test 4: 9 errori (limite teorico) ---\n");
    unsigned char test3[N];
    memcpy(test3, codeword, N);
    for (int i = 0; i < 9; i++)
    {
        test3[i * 4] ^= (0x11 * (i + 1));
    }
    result = rs_decode(test3, decoded);
    if (result >= 0)
    {
        printf("Errori corretti: %d\n", result);
        int ok = (memcmp(decoded, data, K) == 0);
        printf("Verifica: %s\n", ok ? "OK" : "ERRORE");
    }
    else
    {
        printf("Decodifica FALLITA\n");
    }

                 getch();
    return 0;
}
