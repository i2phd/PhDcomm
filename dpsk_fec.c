#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
#include "reedsolomon.h"
#include "fft-complex.h"

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
#define CONV_K 11
#define CONV_RATE_NUM 1
#define CONV_RATE_DEN 2
#define CONV_STATES (1 << (CONV_K-1))
#define CONV_OUTPUTS 2
#define VITERBI_TRACEBACK 64

// Reed-Solomon parameters from LLVM implementation
#define mm 8
#define nn 255
#define tt 16
#define kk (nn - 2 * tt)

#define RS_N nn
#define RS_K kk
#define RS_T tt

// Scrambling LFSR
#define SCRAMBLER_POLY 0x48F
#define SCRAMBLER_INIT 0x7FF

// Polinomi generatori convoluzionali (NASA standard rate 1/2, K=11)
#define CONV_POLY_1 0x4F3
#define CONV_POLY_2 0x6ED

// ============ STRUTTURE DATI ============
static const int barker_13[13] = {1,1,1,1,1,-1,-1,1,1,-1,1,-1,1};
static int preamble_bits[PREAMBLE_LENGTH];

typedef struct {
    uint16_t state;
    uint16_t poly;
} scrambler_t;

typedef struct {
    uint32_t state;
    uint32_t poly1, poly2;
} conv_encoder_t;

typedef struct {
    uint32_t metric;
} viterbi_node_t;

typedef struct {
    viterbi_node_t states[CONV_STATES];
    viterbi_node_t next_states[CONV_STATES];
} viterbi_decoder_t;

typedef struct {
    double phase;
    double frequency;
    double phase_error;
    double freq_error;
    double kp, ki;
} pll_t;

typedef struct {
    double *x_buf, *y_buf;
    double *b_coeff, *a_coeff;
    int order, buf_idx;
} iir_filter_t;

typedef struct {
    double complex *samples;
    double complex *filtered_samples;
    int length;
    double phase_ref;
    int symbol_start;
    int is_synchronized;
    double complex *preamble_samples;
    int preamble_len;
    int sync_pos;
    double quality_metric;
    pll_t carrier_pll;
    iir_filter_t *lpf_i, *lpf_q;
    scrambler_t scrambler;
    conv_encoder_t conv_enc;
    viterbi_decoder_t *viterbi_dec;
    double *correlation_buffer;
    double *phase_history;
    int phase_history_idx;
    double noise_power;
    double signal_power;
    uint8_t *fec_buffer;
    int fec_buffer_size;
} robust_dpsk_system_t;

// ============ SCRAMBLER ============
void scrambler_init(scrambler_t *scr, uint16_t seed, uint16_t poly) { scr->state = seed; scr->poly = poly; }
uint8_t scrambler_next_bit(scrambler_t *scr) {
    uint8_t output = scr->state & 1;
    uint8_t feedback = 0;
    uint16_t temp = scr->state & scr->poly;
    while (temp) { feedback ^= temp & 1; temp >>= 1; }
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
void conv_encoder_init(conv_encoder_t *enc) { enc->state = 0; enc->poly1 = CONV_POLY_1; enc->poly2 = CONV_POLY_2; }
void conv_encode_bit(conv_encoder_t *enc, uint8_t input_bit, uint8_t *output_bits) {
    enc->state = ((enc->state << 1) | input_bit) & ((1 << (CONV_K-1)) - 1);
    uint32_t temp1 = enc->state & enc->poly1;
    uint32_t temp2 = enc->state & enc->poly2;
    output_bits[0] = 0; output_bits[1] = 0;
    while (temp1) { output_bits[0] ^= temp1 & 1; temp1 >>= 1; }
    while (temp2) { output_bits[1] ^= temp2 & 1; temp2 >>= 1; }
}

// ============ DECODER VITERBI ============
uint8_t hamming_distance(uint8_t a, uint8_t b) {
    uint8_t diff = a ^ b; uint8_t count = 0;
    while (diff) { count += diff & 1; diff >>= 1; }
    return count;
}
viterbi_decoder_t* viterbi_decoder_init() {
    viterbi_decoder_t *dec = malloc(sizeof(viterbi_decoder_t));
    for (int i = 0; i < CONV_STATES; i++) { dec->states[i].metric = (i == 0) ? 0 : 0xFFFFFFFF; }
    return dec;
}
void viterbi_decoder_free(viterbi_decoder_t *dec) { if (dec) { free(dec); } }
void viterbi_update_trellis(viterbi_decoder_t *dec, uint8_t received_bits[2], uint16_t *trellis_column) {
    for (int i = 0; i < CONV_STATES; i++) { dec->next_states[i].metric = 0xFFFFFFFF; }

    for (int current_state = 0; current_state < CONV_STATES; current_state++) {
        if (dec->states[current_state].metric == 0xFFFFFFFF) continue;

        for (uint8_t input_bit = 0; input_bit < 2; input_bit++) {
            int next_state = ((current_state << 1) | input_bit) & (CONV_STATES - 1);

            conv_encoder_t temp_enc = { .state = current_state, .poly1 = CONV_POLY_1, .poly2 = CONV_POLY_2 };
            uint8_t expected_bits[2];
            conv_encode_bit(&temp_enc, input_bit, expected_bits);
            uint8_t branch_metric = hamming_distance((received_bits[0] << 1) | received_bits[1], (expected_bits[0] << 1) | expected_bits[1]);
            uint32_t new_metric = dec->states[current_state].metric + branch_metric;

            if (new_metric < dec->next_states[next_state].metric) {
                dec->next_states[next_state].metric = new_metric;
                trellis_column[next_state] = current_state;
            }
        }
    }
    memcpy(dec->states, dec->next_states, sizeof(dec->states));
}

void viterbi_traceback(viterbi_decoder_t *dec, uint16_t *trellis, int num_bits, uint8_t *decoded_bits) {
    int best_state = 0;
    for (int i = 1; i < CONV_STATES; i++) {
        if (dec->states[i].metric < dec->states[best_state].metric) {
            best_state = i;
        }
    }

    for (int i = num_bits; i > 0; i--) {
        decoded_bits[i-1] = best_state & 1;
        best_state = trellis[(i-1) * CONV_STATES + best_state];
    }
}

// ============ REED-SOLOMON ============
void rs_init() {
    init_galois_tables();
    initialize_ecc(NPAR);
}
void rs_encode(uint8_t *data_in, int nbytes, uint8_t *codeword) {
    encode_data(data_in, nbytes, codeword);
}
int rs_decode(uint8_t *codeword, int nbytes) {
    decode_data(codeword, nbytes);
    if (check_syndrome() != 0) {
        return correct_errors_erasures(codeword, nbytes, 0, NULL);
    }
    return 1; // No errors
}

// ============ SISTEMA PRINCIPALE ============
int correlate(double complex *signal, int signal_len, double complex *preamble, int preamble_len) {
    int n = 1;
    while (n < signal_len + preamble_len - 1) {
        n *= 2;
    }

    double complex *signal_padded = calloc(n, sizeof(double complex));
    double complex *preamble_padded = calloc(n, sizeof(double complex));
    memcpy(signal_padded, signal, signal_len * sizeof(double complex));
    memcpy(preamble_padded, preamble, preamble_len * sizeof(double complex));

    Fft_transform(signal_padded, n, false);
    Fft_transform(preamble_padded, n, false);

    for (int i = 0; i < n; i++) {
        signal_padded[i] *= conj(preamble_padded[i]);
    }

    Fft_transform(signal_padded, n, true);

    double max_corr = 0;
    int max_idx = 0;
    for (int i = 0; i < n; i++) {
        double corr = cabs(signal_padded[i]);
        if (corr > max_corr) {
            max_corr = corr;
            max_idx = i;
        }
    }

    free(signal_padded);
    free(preamble_padded);

    return max_idx;
}

void generate_m_sequence(uint8_t *seq, int len, uint16_t poly, int degree) {
    uint16_t lfsr = 1;
    for (int i = 0; i < len; i++) {
        seq[i] = lfsr & 1;
        uint16_t bit = 0;
        uint16_t temp = lfsr & poly;
        while (temp) {
            bit ^= temp & 1;
            temp >>= 1;
        }
        lfsr = (lfsr >> 1) | (bit << (degree - 1));
    }
}

void generate_extended_preamble() { for (int i = 0; i < PREAMBLE_LENGTH; i++) { preamble_bits[i] = barker_13[i % 13]; } }
void pll_init(pll_t *pll, double bandwidth) { *pll = (pll_t){ .kp = 4.0 * bandwidth, .ki = 4.0 * bandwidth * bandwidth }; }
void pll_update(pll_t *pll, double phase_error) {
    pll->phase_error = phase_error;
    pll->freq_error += pll->ki * phase_error;
    pll->frequency = pll->freq_error + pll->kp * phase_error;
    pll->phase += pll->frequency;
    while (pll->phase > PI) pll->phase -= 2*PI;
    while (pll->phase < -PI) pll->phase += 2*PI;
}
iir_filter_t* iir_filter_init(int order) {
    iir_filter_t *filter = malloc(sizeof(iir_filter_t));
    *filter = (iir_filter_t){ .order = order, .x_buf = calloc(order + 1, sizeof(double)), .y_buf = calloc(order + 1, sizeof(double)), .b_coeff = malloc((order + 1) * sizeof(double)), .a_coeff = malloc((order + 1) * sizeof(double)) };
    double fc = (double)CARRIER_FREQ / 4.0 / SAMPLE_RATE;
    double omega = tan(PI * fc);
    double k = 1.0 + sqrt(2.0)*omega + omega*omega;
    filter->b_coeff[0] = omega*omega/k; filter->b_coeff[1] = 2*omega*omega/k; filter->b_coeff[2] = omega*omega/k;
    filter->a_coeff[0] = 1.0; filter->a_coeff[1] = (2*(omega*omega-1))/k; filter->a_coeff[2] = (1-sqrt(2)*omega+omega*omega)/k;
    return filter;
}
double iir_filter_process(iir_filter_t *filter, double input) {
    memmove(filter->x_buf + 1, filter->x_buf, filter->order * sizeof(double));
    memmove(filter->y_buf + 1, filter->y_buf, filter->order * sizeof(double));
    filter->x_buf[0] = input;
    double output = 0.0;
    for (int i = 0; i <= filter->order; i++) { output += filter->b_coeff[i] * filter->x_buf[i]; }
    for (int i = 1; i <= filter->order; i++) { output -= filter->a_coeff[i] * filter->y_buf[i]; }
    filter->y_buf[0] = output;
    return output;
}
void iir_filter_free(iir_filter_t *filter) { if (filter) { free(filter->x_buf); free(filter->y_buf); free(filter->b_coeff); free(filter->a_coeff); free(filter); } }

robust_dpsk_system_t* system_init(int buffer_length) {
    robust_dpsk_system_t *sys = malloc(sizeof(robust_dpsk_system_t));
    *sys = (robust_dpsk_system_t){ .samples = calloc(buffer_length, sizeof(double complex)), .filtered_samples = calloc(buffer_length, sizeof(double complex)), .correlation_buffer = calloc(buffer_length, sizeof(double)), .phase_history = calloc(100, sizeof(double)), .length = buffer_length, .fec_buffer_size = 10000, .fec_buffer = calloc(10000, sizeof(uint8_t)) };
    pll_init(&sys->carrier_pll, PLL_BANDWIDTH);
    sys->lpf_i = iir_filter_init(2); sys->lpf_q = iir_filter_init(2);
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    conv_encoder_init(&sys->conv_enc);
    sys->viterbi_dec = viterbi_decoder_init();
    rs_init();

    sys->preamble_len = 1023;
    uint8_t m_seq[sys->preamble_len];
    generate_m_sequence(m_seq, sys->preamble_len, 0x409, 10);
    sys->preamble_samples = malloc(sys->preamble_len * SAMPLES_PER_SYMBOL * sizeof(double complex));
    for (int i=0; i<sys->preamble_len; i++) {
        double phase = m_seq[i] * PI;
        for (int j=0; j<SAMPLES_PER_SYMBOL; j++) {
            double t = (double)j / SAMPLE_RATE;
            sys->preamble_samples[i*SAMPLES_PER_SYMBOL + j] = cexp(I * (2*PI*CARRIER_FREQ*t + phase));
        }
    }

    printf("Sistema DPSK con FEC concatenato inizializzato\n- Convoluzionale K=%d, Rate=%d/%d\n- Reed-Solomon (%d,%d), T=%d\n- Scrambling con polinomio 0x%X\n", CONV_K, CONV_RATE_NUM, CONV_RATE_DEN, RS_N, RS_K, RS_T, SCRAMBLER_POLY);
    return sys;
}

void system_free(robust_dpsk_system_t *sys) {
    if (sys) {
        free(sys->samples); free(sys->filtered_samples); free(sys->correlation_buffer); free(sys->phase_history);
        iir_filter_free(sys->lpf_i); iir_filter_free(sys->lpf_q);
        viterbi_decoder_free(sys->viterbi_dec);
        free(sys->fec_buffer);
        free(sys->preamble_samples);
        free(sys);
    }
}

int encode_data_with_fec(robust_dpsk_system_t *sys, uint8_t *input_data, int input_len, uint8_t *output_bits, int max_output_bits) {
    int output_pos = 0;
    printf("Encoding %d byte di dati...\n", input_len);
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    for (int block = 0; block < (input_len + RS_K - 1) / RS_K; block++) {
        int block_start = block * RS_K;
        int block_len = (block_start + RS_K <= input_len) ? RS_K : input_len - block_start;
        uint8_t rs_input[RS_K];
        uint8_t rs_codeword[RS_N];
        memset(rs_input, 0, RS_K);
        memcpy(rs_input, &input_data[block_start], block_len);

        rs_encode(rs_input, RS_K, rs_codeword);

        conv_encoder_init(&sys->conv_enc);
        for (int i = 0; i < RS_N; i++) {
            for (int bit = 0; bit < 8; bit++) {
                if (output_pos + 1 >= max_output_bits) { printf("Buffer output troppo piccolo!\n"); return output_pos; }
                uint8_t input_bit = (rs_codeword[i] >> bit) & 1;
                uint8_t conv_output[2];
                conv_encode_bit(&sys->conv_enc, input_bit, conv_output);

                // Scramble the bits
                uint8_t scrambled_bits[2];
                scrambled_bits[0] = conv_output[0] ^ scrambler_next_bit(&sys->scrambler);
                scrambled_bits[1] = conv_output[1] ^ scrambler_next_bit(&sys->scrambler);

                output_bits[output_pos++] = scrambled_bits[0];
                output_bits[output_pos++] = scrambled_bits[1];
            }
        }
        // No flushing here, it's done in the convolutional encoder test
    }
    printf("Encoded: %d bit output da %d byte input\n", output_pos, input_len);
    return output_pos;
}

void synchronize(robust_dpsk_system_t *sys) {
    int preamble_samples_len = sys->preamble_len * SAMPLES_PER_SYMBOL;
    sys->sync_pos = correlate(sys->samples, sys->length, sys->preamble_samples, preamble_samples_len);
    sys->is_synchronized = 1;
    printf("Synchronization found at position %d\n", sys->sync_pos);
}

int decode_data_with_fec(robust_dpsk_system_t *sys, uint8_t *input_bits, int input_len, uint8_t *output_data, int max_output_len, int *corrected_blocks, int *uncorrectable_blocks) {
    if (!sys->is_synchronized) {
        synchronize(sys);
    }
    int output_pos = 0;
    printf("Decoding %d bit ricevuti...\n", input_len);
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    uint8_t descrambled_bits[input_len];
    for (int i=0; i<input_len; i++) {
        descrambled_bits[i] = input_bits[i] ^ scrambler_next_bit(&sys->scrambler);
    }

    int num_columns = input_len / 2;
    uint16_t *trellis = calloc(CONV_STATES * num_columns, sizeof(uint16_t));
    for (int i = 0; i < num_columns; i++) {
        viterbi_update_trellis(sys->viterbi_dec, &descrambled_bits[i * 2], &trellis[i * CONV_STATES]);
    }
    uint8_t *viterbi_output = malloc(num_columns);
    viterbi_traceback(sys->viterbi_dec, trellis, num_columns, viterbi_output);
    free(trellis);
    int viterbi_len = num_columns;
    if (viterbi_len < 0) viterbi_len = 0;
    int byte_count = viterbi_len / 8;
    uint8_t *rs_input = malloc(byte_count);
    for (int i = 0; i < byte_count; i++) {
        rs_input[i] = 0;
        for (int bit = 0; bit < 8; bit++) {
            if (i * 8 + bit < viterbi_len) { rs_input[i] |= viterbi_output[i * 8 + bit] << bit; }
        }
    }
    int rs_blocks = byte_count / RS_N;
    *corrected_blocks = 0;
    *uncorrectable_blocks = 0;
    for (int block = 0; block < rs_blocks; block++) {
        uint8_t *codeword = &rs_input[block * RS_N];
        int rs_result = rs_decode(codeword, RS_N);
        if (rs_result >= 0) { // 0 for no errors, 1 for corrected errors
            if (rs_result == 1) {
                printf("Blocco %d: errori correggibili trovati e corretti.\n", block);
                (*corrected_blocks)++;
            }
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) {
                output_data[output_pos++] = codeword[i];
            }
        } else { // -1 for uncorrectable errors
            printf("Errore Reed-Solomon non correggibile nel blocco %d\n", block);
            (*uncorrectable_blocks)++;
            // Even with uncorrectable errors, we might still want to output the data
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) {
                output_data[output_pos++] = codeword[i];
            }
        }
    }
    free(viterbi_output);
    free(rs_input);
    printf("Decoded: %d byte output da %d bit input\n", output_pos, input_len);
    return output_pos;
}

int test_scrambler() {
    printf("\n====== Inizio Test Scrambler ======\n");
    scrambler_t scrambler, descrambler;
    scrambler_init(&scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    scrambler_init(&descrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);

    uint8_t original_data[100];
    uint8_t scrambled_data[100];
    uint8_t descrambled_data[100];

    for (int i = 0; i < 100; i++) {
        original_data[i] = i;
    }

    for (int i = 0; i < 100; i++) {
        scrambled_data[i] = scramble_byte(&scrambler, original_data[i]);
    }

    for (int i = 0; i < 100; i++) {
        descrambled_data[i] = scramble_byte(&descrambler, scrambled_data[i]);
    }

    int success = 1;
    for (int i = 0; i < 100; i++) {
        if (original_data[i] != descrambled_data[i]) {
            printf("Errore di verifica alla posizione %d: atteso %02X, ottenuto %02X\n", i, original_data[i], descrambled_data[i]);
            success = 0;
        }
    }

    if (success) {
        printf("====== Test Scrambler PASSATO! ======\n");
    } else {
        printf("====== Test Scrambler FALLITO! ======\n");
    }
    return success ? 0 : 1;
}

int test_convolutional_coder() {
    printf("\n====== Inizio Test Coder Convoluzionale ======\n");
    conv_encoder_t encoder;
    conv_encoder_init(&encoder);
    viterbi_decoder_t *decoder = viterbi_decoder_init();

    int num_bits = 1000;
    int num_columns = num_bits + VITERBI_TRACEBACK - 1;
    uint8_t original_bits[num_bits];
    uint8_t encoded_bits[num_columns * 2];
    uint8_t decoded_bits[num_columns];

    for (int i = 0; i < num_bits; i++) {
        original_bits[i] = rand() % 2;
    }

    for (int i = 0; i < num_bits; i++) {
        conv_encode_bit(&encoder, original_bits[i], &encoded_bits[i * 2]);
    }
    // Flush the encoder
    for (int i = 0; i < VITERBI_TRACEBACK - 1; i++) {
        conv_encode_bit(&encoder, 0, &encoded_bits[(num_bits + i) * 2]);
    }

    // Introduce some errors
    int num_errors = 10;
    printf("Introdotti %d errori di bit.\n", num_errors);
    for (int i = 0; i < num_errors; i++) {
        int loc = rand() % (num_bits * 2);
        encoded_bits[loc] ^= 1;
    }

    uint16_t *trellis = calloc(CONV_STATES * num_columns, sizeof(uint16_t));

    for (int i = 0; i < num_columns; i++) {
        viterbi_update_trellis(decoder, &encoded_bits[i * 2], &trellis[i * CONV_STATES]);
    }

    viterbi_traceback(decoder, trellis, num_columns, decoded_bits);

    int success = 1;
    for (int i = 0; i < num_bits; i++) {
        if (original_bits[i] != decoded_bits[i]) {
            printf("Errore di verifica al bit %d: atteso %d, ottenuto %d\n", i, original_bits[i], decoded_bits[i]);
            success = 0;
        }
    }

    if (success) {
        printf("====== Test Coder Convoluzionale PASSATO! ======\n");
    } else {
        printf("====== Test Coder Convoluzionale FALLITO! ======\n");
    }

    free(trellis);
    viterbi_decoder_free(decoder);
    return success ? 0 : 1;
}

int test_custom_fec_chain() {
    printf("\n====== Inizio Test Catena FEC Custom ======\n");
    int num_bytes = 64; // 512 bits
    uint8_t original_data[num_bytes];
    for (int i = 0; i < num_bytes; i++) {
        original_data[i] = rand() % 256;
    }

    // 1. RS encode
    uint8_t rs_codeword[num_bytes + NPAR];
    rs_encode(original_data, num_bytes, rs_codeword);

    // 2. Convolutional encode
    conv_encoder_t conv_enc;
    conv_encoder_init(&conv_enc);
    int conv_input_len_bits = (num_bytes + NPAR) * 8;
    int conv_output_len_bits_unpadded = (conv_input_len_bits + VITERBI_TRACEBACK - 1) * 2;
    int conv_output_len_bits = (conv_output_len_bits_unpadded + 7) & ~7; // Align to 8 bits
    uint8_t conv_output_bits[conv_output_len_bits];
    memset(conv_output_bits, 0, conv_output_len_bits);

    for (int i = 0; i < (num_bytes + NPAR); i++) {
        for (int bit = 0; bit < 8; bit++) {
            uint8_t input_bit = (rs_codeword[i] >> bit) & 1;
            conv_encode_bit(&conv_enc, input_bit, &conv_output_bits[(i * 8 + bit) * 2]);
        }
    }
    for (int i = 0; i < VITERBI_TRACEBACK - 1; i++) {
        conv_encode_bit(&conv_enc, 0, &conv_output_bits[(conv_input_len_bits + i) * 2]);
    }

    // 3. Scramble
    scrambler_t scrambler;
    scrambler_init(&scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    int scrambled_len_bytes = conv_output_len_bits / 8;
    uint8_t scrambled_data[scrambled_len_bytes];
    for (int i = 0; i < scrambled_len_bytes; i++) {
        uint8_t byte = 0;
        for (int bit = 0; bit < 8; bit++) {
            byte |= conv_output_bits[i * 8 + bit] << bit;
        }
        scrambled_data[i] = scramble_byte(&scrambler, byte);
    }

    // 4. Introduce errors
    int num_errors = 20;
    printf("Introdotti %d errori di bit.\n", num_errors);
    for (int i = 0; i < num_errors; i++) {
        int byte_loc = rand() % scrambled_len_bytes;
        int bit_loc = rand() % 8;
        scrambled_data[byte_loc] ^= (1 << bit_loc);
    }

    // 5. Descramble
    scrambler_init(&scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    uint8_t descrambled_data[scrambled_len_bytes];
    for (int i = 0; i < scrambled_len_bytes; i++) {
        descrambled_data[i] = scramble_byte(&scrambler, scrambled_data[i]);
    }

    // 6. Convolutional decode
    viterbi_decoder_t *viterbi_dec = viterbi_decoder_init();
    int num_columns = conv_output_len_bits / 2;
    uint16_t *trellis = calloc(CONV_STATES * num_columns, sizeof(uint16_t));
    uint8_t conv_decoded_bits[num_columns];
    for (int i = 0; i < num_columns; i++) {
        uint8_t received_bits[2];
        int bit_idx0 = i * 2;
        int bit_idx1 = i * 2 + 1;
        received_bits[0] = (descrambled_data[bit_idx0 / 8] >> (bit_idx0 % 8)) & 1;
        received_bits[1] = (descrambled_data[bit_idx1 / 8] >> (bit_idx1 % 8)) & 1;

        viterbi_update_trellis(viterbi_dec, received_bits, &trellis[i * CONV_STATES]);
    }
    viterbi_traceback(viterbi_dec, trellis, num_columns, conv_decoded_bits);
    free(trellis);
    viterbi_decoder_free(viterbi_dec);

    // 7. RS decode
    uint8_t rs_decoded_codeword[num_bytes + NPAR];
    for (int i = 0; i < (num_bytes + NPAR); i++) {
        uint8_t byte = 0;
        for (int bit = 0; bit < 8; bit++) {
            byte |= conv_decoded_bits[i * 8 + bit] << bit;
        }
        rs_decoded_codeword[i] = byte;
    }
    rs_decode(rs_decoded_codeword, num_bytes + NPAR);

    // 8. Compare
    int success = 1;
    for (int i = 0; i < num_bytes; i++) {
        if (original_data[i] != rs_decoded_codeword[i]) {
            printf("Errore di verifica al byte %d: atteso %02X, ottenuto %02X\n", i, original_data[i], rs_decoded_codeword[i]);
            success = 0;
        }
    }

    if (success) {
        printf("====== Test Catena FEC Custom PASSATO! ======\n");
    } else {
        printf("====== Test Catena FEC Custom FALLITO! ======\n");
    }
    return success ? 0 : 1;
}

int test_synchronization() {
    printf("\n====== Inizio Test Sincronizzazione ======\n");
    robust_dpsk_system_t *sys = system_init(1);

    // 1. Generate signal
    int data_len = 1000;
    int preamble_len_samples = sys->preamble_len * SAMPLES_PER_SYMBOL;
    int signal_len = preamble_len_samples + data_len;
    sys->samples = realloc(sys->samples, signal_len * sizeof(double complex));
    sys->length = signal_len;
    memcpy(sys->samples, sys->preamble_samples, preamble_len_samples * sizeof(double complex));
    // Fill the rest with random data
    for (int i=preamble_len_samples; i<signal_len; i++) {
        sys->samples[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }

    // 2. Add noise
    double signal_power = 0;
    for (int i=0; i<signal_len; i++) {
        signal_power += cabs(sys->samples[i]) * cabs(sys->samples[i]);
    }
    signal_power /= signal_len;
    double noise_power = signal_power / pow(10, -20.0/10.0);
    double noise_stddev = sqrt(noise_power / 2.0);
    for (int i=0; i<signal_len; i++) {
        sys->samples[i] += ( (double)rand()/RAND_MAX * 2.0 - 1.0) * noise_stddev + I * ((double)rand()/RAND_MAX * 2.0 - 1.0) * noise_stddev;
    }

    // 3. Add offset
    int offset = rand() % 100;
    printf("Offset: %d\n", offset);
    memmove(&sys->samples[offset], sys->samples, (signal_len - offset) * sizeof(double complex));
    for (int i=0; i<offset; i++) {
        sys->samples[i] = 0;
    }

    // 4. Synchronize
    synchronize(sys);

    // 5. Verify
    int success = (abs(sys->sync_pos - offset) < 5); // Allow for some tolerance

    if (success) {
        printf("====== Test Sincronizzazione PASSATO! ======\n");
    } else {
        printf("====== Test Sincronizzazione FALLITO! ======\n");
    }

    system_free(sys);
    return success ? 0 : 1;
}

int main() {
    int fec_test_result, scrambler_test_result, conv_test_result;

    scrambler_test_result = test_scrambler();
    conv_test_result = test_convolutional_coder();

    int custom_fec_test_result = test_custom_fec_chain();

    int sync_test_result = test_synchronization();

    printf("\n====== Inizio Test Completo FEC (Scrambling, RS, Conv) ======\n");
    robust_dpsk_system_t *sys = system_init(1);
    if (!sys) { printf("Errore: impossibile inizializzare il sistema.\n"); return 1; }
    uint8_t original_data[RS_K];
    for (int i = 0; i < RS_K; i++) { original_data[i] = i % 251; }
    printf("Preparati %d byte di dati originali.\n", RS_K);
    int max_output_bits = RS_N * 8 * 2 * 2;
    uint8_t *encoded_bits = malloc(max_output_bits);
    int encoded_len = encode_data_with_fec(sys, original_data, RS_K, encoded_bits, max_output_bits);
    if (encoded_len == 0) { printf("Errore durante l'encoding!\n"); free(encoded_bits); system_free(sys); return 1; }
    int num_bit_errors = 5;
    printf("Introduzione di %d errori di bit nel bitstream encodato...\n", num_bit_errors);
    for (int i = 0; i < num_bit_errors; i++) { int loc = rand() % encoded_len; encoded_bits[loc] ^= 1; }
    uint8_t decoded_data[RS_K * 2];
    memset(decoded_data, 0, sizeof(decoded_data));
    int corrected_blocks, uncorrectable_blocks;
    int decoded_len = decode_data_with_fec(sys, encoded_bits, encoded_len, decoded_data, sizeof(decoded_data), &corrected_blocks, &uncorrectable_blocks);
    if (decoded_len == 0 && RS_K > 0) { printf("Errore durante il decoding: nessun byte restituito.\n"); }

    printf("\n--- Statistiche decodifica ---\n");
    printf("Blocchi corretti: %d\n", corrected_blocks);
    printf("Blocchi non correggibili: %d\n", uncorrectable_blocks);
    int success = 1;
    if (decoded_len != RS_K) {
        printf("Errore: lunghezza dei dati decodificati non corretta. Atteso %d, ottenuto %d\n", RS_K, decoded_len);
        success = 0;
    }
    for (int i = 0; i < RS_K; i++) {
        if (decoded_data[i] != original_data[i]) {
            printf("Errore di verifica alla posizione %d: atteso %02X, ottenuto %02X\n", i, original_data[i], decoded_data[i]);
            success = 0;
        }
    }
    if (success) { printf("\n====== Test FEC Completo PASSATO! ======\n"); }
    else { printf("\n====== Test FEC Completo FALLITO! ======\n"); }
    free(encoded_bits);
    system_free(sys);
    fec_test_result = success ? 0 : 1;

    return (fec_test_result || scrambler_test_result || conv_test_result || custom_fec_test_result || sync_test_result);
}

// Implementations of missing complex arithmetic functions
double complex __muldc3(double a, double b, double c, double d) {
    return (a * c - b * d) + I * (a * d + b * c);
}

double complex __divdc3(double a, double b, double c, double d) {
    double denom = c * c + d * d;
    return ((a * c + b * d) / denom) + I * ((b * c - a * d) / denom);
}
