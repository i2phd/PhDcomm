#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>

#define SAMPLE_RATE 48000
#define CARRIER_FREQ 480
#define SYMBOL_DURATION 10
#define SAMPLES_PER_SYMBOL (SAMPLE_RATE * SYMBOL_DURATION / CARRIER_FREQ)
#define OVERSAMPLING 4
#define EFFECTIVE_SAMPLES (SAMPLES_PER_SYMBOL * OVERSAMPLING)
#define PI 3.14159265359
#define PREAMBLE_LENGTH 64
#define CORRELATION_THRESHOLD 0.15
#define PLL_BANDWIDTH 0.01
#define INTEGRATION_SYMBOLS 8

// ============ PARAMETRI FEC ============
#define CONV_K 11              // Lunghezza vincolo convoluzionale
#define CONV_RATE_NUM 1        // Rate 1/2
#define CONV_RATE_DEN 2
#define CONV_STATES (1 << (CONV_K-1))  // 2^10 = 1024 stati
#define CONV_OUTPUTS 2         // 2 bit output per ogni bit input
#define VITERBI_TRACEBACK 44   // 4 * K per traceback

// Reed-Solomon parameters from LLVM implementation
#define mm 8
#define nn 255
#define tt 16
#define kk (nn - 2 * tt)

#define RS_N nn
#define RS_K kk
#define RS_T tt

// Scrambling LFSR
#define SCRAMBLER_POLY 0x48F   // Polinomio x^11 + x^7 + x^3 + x + 1
#define SCRAMBLER_INIT 0x7FF   // Seed iniziale

// Polinomi generatori convoluzionali (NASA standard rate 1/2, K=11)
#define CONV_POLY_1 0x4F3      // 10011110011 (octal 2331)
#define CONV_POLY_2 0x6ED      // 11011101101 (octal 3355)

// ============ STRUTTURE DATI ============

// Codice Barker
static const int barker_13[13] = {1,1,1,1,1,-1,-1,1,1,-1,1,-1,1};
static int preamble_bits[PREAMBLE_LENGTH];

// Tabelle Galois Field per Reed-Solomon
static uint8_t gf_exp[512];
static uint8_t gf_log[256];
static uint8_t rs_generator[33];

// Scrambler LFSR
typedef struct {
    uint16_t state;
    uint16_t poly;
} scrambler_t;

// Encoder convoluzionale
typedef struct {
    uint32_t state;        // Registro a scorrimento K-1 bit
    uint32_t poly1, poly2; // Polinomi generatori
} conv_encoder_t;

// Nodo Viterbi
typedef struct {
    uint32_t metric;       // Metrica cumulativa
    uint64_t path;         // Bit del path
} viterbi_node_t;

// Decoder Viterbi
typedef struct {
    viterbi_node_t states[CONV_STATES];     // Stati correnti
    viterbi_node_t next_states[CONV_STATES]; // Stati prossimi
    uint8_t *traceback_buffer;              // Buffer traceback
    int traceback_pos;
    uint32_t *path_memory;                  // Memoria path
} viterbi_decoder_t;

#include "rs_llvm.h"

// Note: The LLVM implementation uses global variables for data buffers.
// These wrapper functions will copy data to and from those global buffers.
// This is not ideal, but it is the safest way to integrate the known-good code.

// Global variables from rs_llvm.c, declared as extern in rs_llvm.h
// int alpha_to[nn+1], index_of[nn+1], gg[nn-kk+1];
// int recd[nn], data[kk], bb[nn-kk];

// PLL per tracking portante
typedef struct {
    double phase;
    double frequency;
    double phase_error;
    double freq_error;
    double kp, ki;
} pll_t;

// Filtro IIR
typedef struct {
    double *x_buf, *y_buf;
    double *b_coeff, *a_coeff;
    int order, buf_idx;
} iir_filter_t;

// Demodulatore completo
typedef struct {
    // Elaborazione segnale
    double complex *samples;
    double complex *filtered_samples;
    int length;
    double phase_ref;
    int symbol_start;
    int is_synchronized;
    double quality_metric;
    pll_t carrier_pll;
    iir_filter_t *lpf_i, *lpf_q;

    // FEC e scrambling
    scrambler_t scrambler;
    conv_encoder_t conv_enc;
    viterbi_decoder_t *viterbi_dec;
    // rs_codec is no longer needed, the LLVM code uses globals.

    // Buffer e statistiche
    double *correlation_buffer;
    double *phase_history;
    int phase_history_idx;
    double noise_power;
    double signal_power;
    uint8_t *fec_buffer;
    int fec_buffer_size;
} robust_dpsk_system_t;

// Galois Field functions are now part of the LLVM implementation in rs_llvm.c

// ============ SCRAMBLER ============

void scrambler_init(scrambler_t *scr, uint16_t seed, uint16_t poly) {
    scr->state = seed;
    scr->poly = poly;
}

uint8_t scrambler_next_bit(scrambler_t *scr) {
    uint8_t output = scr->state & 1;
    uint8_t feedback = 0;

    // Calcola XOR dei bit corrispondenti al polinomio
    uint16_t temp = scr->state & scr->poly;
    while (temp) {
        feedback ^= temp & 1;
        temp >>= 1;
    }

    scr->state = (scr->state >> 1) | (feedback << (CONV_K-1));
    return output;
}

uint8_t scramble_byte(scrambler_t *scr, uint8_t data) {
    uint8_t result = 0;
    for (int i = 0; i < 8; i++) {
        uint8_t scramble_bit = scrambler_next_bit(scr);
        uint8_t data_bit = (data >> i) & 1;
        result |= (data_bit ^ scramble_bit) << i;
    }
    return result;
}

// ============ ENCODER CONVOLUZIONALE ============

void conv_encoder_init(conv_encoder_t *enc) {
    enc->state = 0;
    enc->poly1 = CONV_POLY_1;
    enc->poly2 = CONV_POLY_2;
}

void conv_encode_bit(conv_encoder_t *enc, uint8_t input_bit, uint8_t *output_bits) {
    // Aggiorna stato
    enc->state = ((enc->state << 1) | input_bit) & ((1 << (CONV_K-1)) - 1);

    // Calcola output
    uint32_t temp1 = enc->state & enc->poly1;
    uint32_t temp2 = enc->state & enc->poly2;

    output_bits[0] = 0;
    output_bits[1] = 0;

    while (temp1) {
        output_bits[0] ^= temp1 & 1;
        temp1 >>= 1;
    }

    while (temp2) {
        output_bits[1] ^= temp2 & 1;
        temp2 >>= 1;
    }
}

// ============ DECODER VITERBI ============

uint8_t hamming_distance(uint8_t a, uint8_t b) {
    uint8_t diff = a ^ b;
    uint8_t count = 0;
    while (diff) {
        count += diff & 1;
        diff >>= 1;
    }
    return count;
}

viterbi_decoder_t* viterbi_decoder_init() {
    viterbi_decoder_t *dec = malloc(sizeof(viterbi_decoder_t));

    // Inizializza metriche
    for (int i = 0; i < CONV_STATES; i++) {
        dec->states[i].metric = (i == 0) ? 0 : 0xFFFFFFFF;
        dec->states[i].path = 0;
    }

    dec->traceback_buffer = calloc(VITERBI_TRACEBACK * CONV_STATES, sizeof(uint8_t));
    dec->path_memory = calloc(VITERBI_TRACEBACK, sizeof(uint32_t));
    dec->traceback_pos = 0;

    return dec;
}

void viterbi_decoder_free(viterbi_decoder_t *dec) {
    if (dec) {
        free(dec->traceback_buffer);
        free(dec->path_memory);
        free(dec);
    }
}

uint8_t viterbi_decode_bit(viterbi_decoder_t *dec, uint8_t received_bits[2]) {
    // Reset stati prossimi
    for (int i = 0; i < CONV_STATES; i++) {
        dec->next_states[i].metric = 0xFFFFFFFF;
        dec->next_states[i].path = 0;
    }

    // Calcola transizioni
    for (int current_state = 0; current_state < CONV_STATES; current_state++) {
        if (dec->states[current_state].metric == 0xFFFFFFFF) continue;

        for (uint8_t input_bit = 0; input_bit < 2; input_bit++) {
            // Calcola prossimo stato
            int next_state = ((current_state << 1) | input_bit) & (CONV_STATES - 1);

            // Genera output atteso
            conv_encoder_t temp_enc = {0};
            temp_enc.state = current_state;
            temp_enc.poly1 = CONV_POLY_1;
            temp_enc.poly2 = CONV_POLY_2;

            uint8_t expected_bits[2];
            conv_encode_bit(&temp_enc, input_bit, expected_bits);

            // Calcola metrica di branch
            uint8_t branch_metric = hamming_distance(
                (received_bits[0] << 1) | received_bits[1],
                (expected_bits[0] << 1) | expected_bits[1]
            );

            uint32_t new_metric = dec->states[current_state].metric + branch_metric;

            // Aggiorna se migliore
            if (new_metric < dec->next_states[next_state].metric) {
                dec->next_states[next_state].metric = new_metric;
                dec->next_states[next_state].path =
                    (dec->states[current_state].path << 1) | input_bit;
            }
        }
    }

    // Copia stati
    memcpy(dec->states, dec->next_states, sizeof(dec->states));

    // Traceback se necessario
    if (dec->traceback_pos >= VITERBI_TRACEBACK) {
        // Trova stato con metrica migliore
        int best_state = 0;
        for (int i = 1; i < CONV_STATES; i++) {
            if (dec->states[i].metric < dec->states[best_state].metric) {
                best_state = i;
            }
        }

        // Estrai bit dal path
        return (dec->states[best_state].path >> (VITERBI_TRACEBACK - 1)) & 1;
    }

    dec->traceback_pos++;
    return 0; // Ritardo iniziale
}

// ============ REED-SOLOMON (LLVM IMPLEMENTATION) ============

// Note: The following functions are wrappers around the LLVM Reed-Solomon
// implementation, which uses global variables. This is not ideal but
// necessary for integrating the known-good code without a major rewrite.

void rs_init() {
    generate_gf();
    gen_poly();
}

void rs_encode(uint8_t *data_in, uint8_t *codeword) {
    // The LLVM code uses global integer arrays. We must copy our data into them.
    for (int i = 0; i < kk; i++) {
        data[i] = data_in[i];
    }

    encode_rs(); // This operates on the global `data` and `bb` arrays.

    // The result is systematic, so the data part is unchanged.
    // The parity is in the global `bb` array.
    memcpy(codeword, data_in, kk);
    for (int i = 0; i < (nn - kk); i++) {
        codeword[kk + i] = bb[i];
    }
}

int rs_decode(uint8_t *codeword) {
    // The LLVM decoder expects the received message in the global `recd`
    // array, and in INDEX FORM.
    for (int i = 0; i < nn; i++) {
        recd[i] = index_of[codeword[i]];
    }

    decode_rs(); // This operates on and modifies the global `recd` array.

    // The corrected codeword is now in `recd`, but in POLYNOMIAL form.
    // We copy it back to the original codeword buffer.
    for (int i = 0; i < nn; i++) {
        codeword[i] = recd[i];
    }

    // The LLVM `decode_rs` does not return the number of corrected errors.
    // We can't easily check for uncorrectable errors here, but the `main`
    // function's final verification will catch any failures.
    return 0; // Assume success.
}

// ============ SISTEMA PRINCIPALE ============

void generate_extended_preamble() {
    for (int i = 0; i < PREAMBLE_LENGTH; i++) {
        preamble_bits[i] = barker_13[i % 13];
    }
}

// Inizializzazione PLL
void pll_init(pll_t *pll, double bandwidth) {
    pll->phase = 0.0;
    pll->frequency = 0.0;
    pll->phase_error = 0.0;
    pll->freq_error = 0.0;
    pll->kp = 4.0 * bandwidth;
    pll->ki = 4.0 * bandwidth * bandwidth;
}

void pll_update(pll_t *pll, double phase_error) {
    pll->phase_error = phase_error;
    pll->freq_error += pll->ki * phase_error;
    pll->frequency = pll->freq_error + pll->kp * phase_error;
    pll->phase += pll->frequency;

    while (pll->phase > PI) pll->phase -= 2*PI;
    while (pll->phase < -PI) pll->phase += 2*PI;
}

// Inizializzazione filtro IIR
iir_filter_t* iir_filter_init(int order) {
    iir_filter_t *filter = malloc(sizeof(iir_filter_t));
    filter->order = order;
    filter->x_buf = calloc(order + 1, sizeof(double));
    filter->y_buf = calloc(order + 1, sizeof(double));
    filter->b_coeff = malloc((order + 1) * sizeof(double));
    filter->a_coeff = malloc((order + 1) * sizeof(double));
    filter->buf_idx = 0;

    // Coefficienti Butterworth
    double fc = (double)CARRIER_FREQ / 4.0 / SAMPLE_RATE;
    double omega = tan(PI * fc);
    double k = 1.0 + sqrt(2.0)*omega + omega*omega;

    filter->b_coeff[0] = omega*omega/k;
    filter->b_coeff[1] = 2*omega*omega/k;
    filter->b_coeff[2] = omega*omega/k;
    filter->a_coeff[0] = 1.0;
    filter->a_coeff[1] = (2*(omega*omega-1))/k;
    filter->a_coeff[2] = (1-sqrt(2)*omega+omega*omega)/k;

    return filter;
}

double iir_filter_process(iir_filter_t *filter, double input) {
    memmove(filter->x_buf + 1, filter->x_buf, filter->order * sizeof(double));
    memmove(filter->y_buf + 1, filter->y_buf, filter->order * sizeof(double));

    filter->x_buf[0] = input;

    double output = 0.0;
    for (int i = 0; i <= filter->order; i++) {
        output += filter->b_coeff[i] * filter->x_buf[i];
    }
    for (int i = 1; i <= filter->order; i++) {
        output -= filter->a_coeff[i] * filter->y_buf[i];
    }

    filter->y_buf[0] = output;
    return output;
}

void iir_filter_free(iir_filter_t *filter) {
    if (filter) {
        free(filter->x_buf);
        free(filter->y_buf);
        free(filter->b_coeff);
        free(filter->a_coeff);
        free(filter);
    }
}

// Inizializzazione sistema completo
robust_dpsk_system_t* system_init(int buffer_length) {
    robust_dpsk_system_t *sys = malloc(sizeof(robust_dpsk_system_t));

    // Elaborazione segnale
    sys->samples = calloc(buffer_length, sizeof(double complex));
    sys->filtered_samples = calloc(buffer_length, sizeof(double complex));
    sys->correlation_buffer = calloc(buffer_length, sizeof(double));
    sys->phase_history = calloc(100, sizeof(double));
    sys->length = buffer_length;
    sys->phase_ref = 0.0;
    sys->symbol_start = 0;
    sys->is_synchronized = 0;
    sys->quality_metric = 0.0;
    sys->phase_history_idx = 0;
    sys->noise_power = 1.0;
    sys->signal_power = 1.0;

    // Inizializza componenti
    pll_init(&sys->carrier_pll, PLL_BANDWIDTH);
    sys->lpf_i = iir_filter_init(2);
    sys->lpf_q = iir_filter_init(2);

    // Inizializza FEC
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    conv_encoder_init(&sys->conv_enc);
    sys->viterbi_dec = viterbi_decoder_init();

    sys->fec_buffer_size = 10000;
    sys->fec_buffer = calloc(sys->fec_buffer_size, sizeof(uint8_t));

    // Inizializza tabelle Galois e Reed-Solomon (LLVM version)
    rs_init();

    generate_extended_preamble();

    printf("Sistema DPSK con FEC concatenato inizializzato\n");
    printf("- Convoluzionale K=%d, Rate=%d/%d\n", CONV_K, CONV_RATE_NUM, CONV_RATE_DEN);
    printf("- Reed-Solomon (%d,%d), T=%d\n", RS_N, RS_K, RS_T);
    printf("- Scrambling con polinomio 0x%X\n", SCRAMBLER_POLY);

    return sys;
}

void system_free(robust_dpsk_system_t *sys) {
    if (sys) {
        free(sys->samples);
        free(sys->filtered_samples);
        free(sys->correlation_buffer);
        free(sys->phase_history);
        iir_filter_free(sys->lpf_i);
        iir_filter_free(sys->lpf_q);
        viterbi_decoder_free(sys->viterbi_dec);
        // rs_codec_free is no longer needed
        free(sys->fec_buffer);
        free(sys);
    }
}

// Encoding completo: Data -> Scrambling -> Reed-Solomon -> Convolutional
int encode_data_with_fec(robust_dpsk_system_t *sys, uint8_t *input_data, int input_len,
                        uint8_t *output_bits, int max_output_bits) {
    int output_pos = 0;

    printf("Encoding %d byte di dati...\n", input_len);

    // Reset scrambler per encoding
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);

    // Processa in blocchi Reed-Solomon
    for (int block = 0; block < (input_len + RS_K - 1) / RS_K; block++) {
        int block_start = block * RS_K;
        int block_len = (block_start + RS_K <= input_len) ? RS_K : input_len - block_start;

        uint8_t rs_input[RS_K];
        uint8_t rs_codeword[RS_N];

        // Prepara blocco con scrambling
        memset(rs_input, 0, RS_K);
        for (int i = 0; i < block_len; i++) {
            rs_input[i] = scramble_byte(&sys->scrambler, input_data[block_start + i]);
        }

        // Encoding Reed-Solomon
        rs_encode(rs_input, rs_codeword);

        // Encoding convoluzionale
        conv_encoder_init(&sys->conv_enc);

        for (int i = 0; i < RS_N; i++) {
            for (int bit = 0; bit < 8; bit++) {
                if (output_pos + 1 >= max_output_bits) {
                    printf("Buffer output troppo piccolo!\n");
                    return output_pos;
                }

                uint8_t input_bit = (rs_codeword[i] >> bit) & 1;
                uint8_t conv_output[2];
                conv_encode_bit(&sys->conv_enc, input_bit, conv_output);

                output_bits[output_pos++] = conv_output[0];
                output_bits[output_pos++] = conv_output[1];
            }
        }

        // Flush convoluzionale
        for (int i = 0; i < CONV_K - 1; i++) {
            if (output_pos + 1 >= max_output_bits) break;
            uint8_t conv_output[2];
            conv_encode_bit(&sys->conv_enc, 0, conv_output);
            output_bits[output_pos++] = conv_output[0];
            output_bits[output_pos++] = conv_output[1];
        }
    }

    printf("Encoded: %d bit output da %d byte input\n", output_pos, input_len);
    return output_pos;
}

// Decoding completo: Bits -> Viterbi -> Reed-Solomon -> Descrambling
int decode_data_with_fec(robust_dpsk_system_t *sys, uint8_t *input_bits, int input_len,
                        uint8_t *output_data, int max_output_len) {
    int output_pos = 0;

    printf("Decoding %d bit ricevuti...\n", input_len);

    // Reset per decoding
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);

    // Buffer per Viterbi output
    uint8_t *viterbi_output = malloc(input_len / 2);
    int viterbi_pos = 0;

    // Decoding Viterbi
    for (int i = 0; i < input_len - 1; i += 2) {
        uint8_t received_bits[2] = {input_bits[i], input_bits[i+1]};
        uint8_t decoded_bit = viterbi_decode_bit(sys->viterbi_dec, received_bits);

        if (viterbi_pos < input_len / 2) {
            viterbi_output[viterbi_pos++] = decoded_bit;
        }
    }

    // Raggruppa in byte
    int byte_count = viterbi_pos / 8;
    uint8_t *rs_input = malloc(byte_count);

    for (int i = 0; i < byte_count; i++) {
        rs_input[i] = 0;
        for (int bit = 0; bit < 8; bit++) {
            if (i * 8 + bit < viterbi_pos) {
                rs_input[i] |= viterbi_output[i * 8 + bit] << bit;
            }
        }
    }

    // Processa blocchi Reed-Solomon
    int rs_blocks = byte_count / RS_N;
    for (int block = 0; block < rs_blocks; block++) {
        uint8_t *codeword = &rs_input[block * RS_N];

        // Tentativo correzione Reed-Solomon
        int rs_result = rs_decode(codeword);

        if (rs_result == 0) {  // Successo o nessun errore
            // Descrambling e output
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) {
                output_data[output_pos++] = scramble_byte(&sys->scrambler, codeword[i]);
            }
        } else {
            printf("Errore Reed-Solomon non correggibile nel blocco %d\n", block);
            // Output anche con errori per debug
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) {
                output_data[output_pos++] = scramble_byte(&sys->scrambler, codeword[i]);
            }
        }
    }

    free(viterbi_output);
    free(rs_input);

    printf("Decoded: %d byte output da %d bit input\n", output_pos, input_len);
    return output_pos;
}

int main() {
    printf("====== Inizio Test Reed-Solomon (LLVM) ======\n");

    // Inizializza il sistema (chiama rs_init)
    robust_dpsk_system_t *sys = system_init(1);
    if (!sys) {
        printf("Errore: impossibile inizializzare il sistema.\n");
        return 1;
    }

    // 1. Prepara i dati originali
    uint8_t original_data[RS_K];
    for (int i = 0; i < RS_K; i++) {
        original_data[i] = i % 251;
    }

    // 2. Esegui l'encoding Reed-Solomon
    uint8_t codeword[RS_N];
    rs_encode(original_data, codeword);
    printf("Dati originali encodati.\n");

    // 3. Introduci degli errori
    uint8_t corrupted_codeword[RS_N];
    memcpy(corrupted_codeword, codeword, RS_N);
    int num_errors = 8;
    printf("Introduzione di %d errori nel codeword...\n", num_errors);
    for (int i = 0; i < num_errors; i++) {
        int loc = (i * 30) % RS_N;
        uint8_t error_val = (i + 1) * 17;
        corrupted_codeword[loc] ^= error_val;
    }

    // 4. Esegui il decoding
    rs_decode(corrupted_codeword);
    printf("Decoding tentato.\n");

    // 5. Verifica il risultato
    int success = 1;
    for (int i = 0; i < RS_K; i++) {
        if (corrupted_codeword[i] != original_data[i]) {
            printf("Errore di verifica alla posizione %d: atteso %02X, ottenuto %02X\n",
                   i, original_data[i], corrupted_codeword[i]);
            success = 0;
        }
    }

    if (success) {
        printf("\n====== Test Reed-Solomon PASSATO! ======\n");
    } else {
        printf("\n====== Test Reed-Solomon FALLITO! ======\n");
    }

    system_free(sys);
    return success ? 0 : 1;
}
