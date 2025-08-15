#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>
#include "rs_llvm.h"

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
#define VITERBI_TRACEBACK 44

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
    uint64_t path;
} viterbi_node_t;

typedef struct {
    viterbi_node_t states[CONV_STATES];
    viterbi_node_t next_states[CONV_STATES];
    uint8_t *traceback_buffer;
    int traceback_pos;
    uint32_t *path_memory;
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
    for (int i = 0; i < CONV_STATES; i++) { dec->states[i].metric = (i == 0) ? 0 : 0xFFFFFFFF; dec->states[i].path = 0; }
    dec->traceback_buffer = calloc(VITERBI_TRACEBACK * CONV_STATES, sizeof(uint8_t));
    dec->path_memory = calloc(VITERBI_TRACEBACK, sizeof(uint32_t));
    dec->traceback_pos = 0;
    return dec;
}
void viterbi_decoder_free(viterbi_decoder_t *dec) { if (dec) { free(dec->traceback_buffer); free(dec->path_memory); free(dec); } }
uint8_t viterbi_decode_bit(viterbi_decoder_t *dec, uint8_t received_bits[2]) {
    for (int i = 0; i < CONV_STATES; i++) { dec->next_states[i].metric = 0xFFFFFFFF; dec->next_states[i].path = 0; }
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
                dec->next_states[next_state].path = (dec->states[current_state].path << 1) | input_bit;
            }
        }
    }
    memcpy(dec->states, dec->next_states, sizeof(dec->states));
    if (dec->traceback_pos >= VITERBI_TRACEBACK) {
        int best_state = 0;
        for (int i = 1; i < CONV_STATES; i++) { if (dec->states[i].metric < dec->states[best_state].metric) best_state = i; }
        return (dec->states[best_state].path >> (VITERBI_TRACEBACK - 1)) & 1;
    }
    dec->traceback_pos++;
    return 0;
}

// ============ REED-SOLOMON (LLVM IMPLEMENTATION) ============
void rs_init() { generate_gf(); gen_poly(); }
void rs_encode(uint8_t *data_in, uint8_t *codeword) {
    for (int i = 0; i < kk; i++) { data[i] = data_in[i]; }
    encode_rs();
    memcpy(codeword, data_in, kk);
    for (int i = 0; i < (nn - kk); i++) { codeword[kk + i] = bb[i]; }
}
int rs_decode(uint8_t *codeword) {
    for (int i = 0; i < nn; i++) { recd[i] = index_of[codeword[i]]; }
    decode_rs();
    for (int i = 0; i < nn; i++) { codeword[i] = recd[i]; }
    return 0;
}

// ============ SISTEMA PRINCIPALE ============
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
    generate_extended_preamble();
    printf("Sistema DPSK con FEC concatenato inizializzato\n- Convoluzionale K=%d, Rate=%d/%d\n- Reed-Solomon (%d,%d), T=%d\n- Scrambling con polinomio 0x%X\n", CONV_K, CONV_RATE_NUM, CONV_RATE_DEN, RS_N, RS_K, RS_T, SCRAMBLER_POLY);
    return sys;
}

void system_free(robust_dpsk_system_t *sys) {
    if (sys) {
        free(sys->samples); free(sys->filtered_samples); free(sys->correlation_buffer); free(sys->phase_history);
        iir_filter_free(sys->lpf_i); iir_filter_free(sys->lpf_q);
        viterbi_decoder_free(sys->viterbi_dec);
        free(sys->fec_buffer); free(sys);
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
        for (int i = 0; i < block_len; i++) { rs_input[i] = scramble_byte(&sys->scrambler, input_data[block_start + i]); }
        rs_encode(rs_input, rs_codeword);
        conv_encoder_init(&sys->conv_enc);
        for (int i = 0; i < RS_N; i++) {
            for (int bit = 0; bit < 8; bit++) {
                if (output_pos + 1 >= max_output_bits) { printf("Buffer output troppo piccolo!\n"); return output_pos; }
                uint8_t input_bit = (rs_codeword[i] >> bit) & 1;
                uint8_t conv_output[2];
                conv_encode_bit(&sys->conv_enc, input_bit, conv_output);
                output_bits[output_pos++] = conv_output[0];
                output_bits[output_pos++] = conv_output[1];
            }
        }
        // Flush convoluzionale with enough bits for Viterbi traceback
        for (int i = 0; i < VITERBI_TRACEBACK - 1; i++) {
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

int decode_data_with_fec(robust_dpsk_system_t *sys, uint8_t *input_bits, int input_len, uint8_t *output_data, int max_output_len) {
    int output_pos = 0;
    printf("Decoding %d bit ricevuti...\n", input_len);
    scrambler_init(&sys->scrambler, SCRAMBLER_INIT, SCRAMBLER_POLY);
    uint8_t *raw_viterbi_output = malloc(input_len / 2);
    int viterbi_pos = 0;
    for (int i = 0; i < input_len - 1; i += 2) {
        uint8_t received_bits[2] = {input_bits[i], input_bits[i+1]};
        uint8_t decoded_bit = viterbi_decode_bit(sys->viterbi_dec, received_bits);
        if (viterbi_pos < input_len / 2) { raw_viterbi_output[viterbi_pos++] = decoded_bit; }
    }
    uint8_t *viterbi_output = &raw_viterbi_output[VITERBI_TRACEBACK - 1];
    int viterbi_len = viterbi_pos - (VITERBI_TRACEBACK - 1);
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
    for (int block = 0; block < rs_blocks; block++) {
        uint8_t *codeword = &rs_input[block * RS_N];
        int rs_result = rs_decode(codeword);
        if (rs_result == 0) {
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) { output_data[output_pos++] = scramble_byte(&sys->scrambler, codeword[i]); }
        } else {
            printf("Errore Reed-Solomon non correggibile nel blocco %d\n", block);
            for (int i = 0; i < RS_K && output_pos < max_output_len; i++) { output_data[output_pos++] = scramble_byte(&sys->scrambler, codeword[i]); }
        }
    }
    free(raw_viterbi_output);
    free(rs_input);
    printf("Decoded: %d byte output da %d bit input\n", output_pos, input_len);
    return output_pos;
}

int main() {
    printf("====== Inizio Test Completo FEC (Scrambling, RS, Conv) ======\n");
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
    int decoded_len = decode_data_with_fec(sys, encoded_bits, encoded_len, decoded_data, sizeof(decoded_data));
    if (decoded_len == 0 && RS_K > 0) { printf("Errore durante il decoding: nessun byte restituito.\n"); }
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
    return success ? 0 : 1;
}
