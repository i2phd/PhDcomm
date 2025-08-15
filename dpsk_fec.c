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

// Reed-Solomon parametri (255,223) - 32 byte correzione
#define RS_N 255               // Lunghezza codeword
#define RS_K 223               // Lunghezza messaggio
#define RS_T 16                // Capacità correzione (32 byte errori)

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

// Reed-Solomon encoder/decoder
typedef struct {
    uint8_t *generator;
    uint8_t *syndromes;
    int *error_positions;
    uint8_t *error_values;
} rs_codec_t;

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
    rs_codec_t *rs_codec;

    // Buffer e statistiche
    double *correlation_buffer;
    double *phase_history;
    int phase_history_idx;
    double noise_power;
    double signal_power;
    uint8_t *fec_buffer;
    int fec_buffer_size;
} robust_dpsk_system_t;

// ============ FUNZIONI GALOIS FIELD ============

uint8_t gf_add(uint8_t a, uint8_t b) {
    return a ^ b;
}

uint8_t gf_mult(uint8_t a, uint8_t b) {
    if (a == 0 || b == 0) return 0;
    return gf_exp[(gf_log[a] + gf_log[b]) % 255];
}

uint8_t gf_div(uint8_t a, uint8_t b) {
    if (a == 0) return 0;
    if (b == 0) return 0; // Errore divisione per zero
    return gf_exp[(gf_log[a] - gf_log[b] + 255) % 255];
}

void init_galois_field() {
    int poly = 0x11D; // x^8 + x^4 + x^3 + x^2 + 1

    // Inizializza tabelle exp e log
    gf_exp[0] = 1;
    gf_log[0] = 0; // Non definito, ma serve per evitare errori

    for (int i = 1; i < 255; i++) {
        gf_exp[i] = gf_exp[i-1] * 2;
        if (gf_exp[i] >= 256) {
            gf_exp[i] ^= poly;
        }
        gf_log[gf_exp[i]] = i;
    }

    // Estendi tabella exp per evitare modulo
    for (int i = 255; i < 512; i++) {
        gf_exp[i] = gf_exp[i - 255];
    }

    gf_log[1] = 0;
}

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

// ============ REED-SOLOMON ============

void rs_generate_polynomial() {
    // Genera polinomio generatore per RS(255,223)
    rs_generator[0] = 1;

    for (int i = 1; i <= 2 * RS_T; i++) {
        rs_generator[i] = 1;
        for (int j = i; j > 0; j--) {
            rs_generator[j] = gf_mult(rs_generator[j], gf_exp[i]) ^ rs_generator[j-1];
        }
        rs_generator[0] = gf_mult(rs_generator[0], gf_exp[i]);
    }
}

rs_codec_t* rs_codec_init() {
    rs_codec_t *codec = malloc(sizeof(rs_codec_t));
    codec->generator = malloc(33 * sizeof(uint8_t));
    codec->syndromes = malloc(33 * sizeof(uint8_t));
    codec->error_positions = malloc(RS_T * sizeof(int));
    codec->error_values = malloc(RS_T * sizeof(uint8_t));

    memcpy(codec->generator, rs_generator, 33);
    return codec;
}

void rs_codec_free(rs_codec_t *codec) {
    if (codec) {
        free(codec->generator);
        free(codec->syndromes);
        free(codec->error_positions);
        free(codec->error_values);
        free(codec);
    }
}

void rs_encode(rs_codec_t *codec, uint8_t *data, uint8_t *codeword) {
    // Copia dati
    memcpy(codeword, data, RS_K);
    memset(codeword + RS_K, 0, RS_N - RS_K);

    // Divisione sintetica
    for (int i = 0; i < RS_K; i++) {
        uint8_t feedback = codeword[i];
        if (feedback != 0) {
            for (int j = 1; j <= RS_N - RS_K; j++) {
                codeword[i + j] ^= gf_mult(codec->generator[j], feedback);
            }
        }
    }

    // Ripristina dati originali
    memcpy(codeword, data, RS_K);
}

int rs_calculate_syndromes(rs_codec_t *codec, uint8_t *codeword) {
    int error_count = 0;

    int fcr = 1; // First consecutive root is alpha^1
    for (int i = 0; i < 2 * RS_T; i++) {
        uint8_t s = 0;
        for (int j = 0; j < RS_N; j++) {
            s ^= gf_mult(codeword[j], gf_exp[((i + fcr) * j) % 255]);
        }
        codec->syndromes[i] = s;
        if (s != 0) error_count++;
    }

    return error_count;
}

#ifndef min
#define min(a,b) ((a) < (b) ? (a) : (b))
#endif

int rs_decode(rs_codec_t* codec, uint8_t* codeword)
{
    int nroots = 2 * RS_T;

    if (rs_calculate_syndromes(codec, codeword) == 0) {
        return 0; // No errors
    }
    printf("Syndromes: ");
    for(int i=0; i<nroots; i++) printf("%02X ", codec->syndromes[i]);
    printf("\n");

    // Berlekamp-Massey Algorithm (from Schifra)
    uint8_t lambda[nroots + 1];
    uint8_t prev_lambda[nroots + 1];

    memset(lambda, 0, sizeof(lambda));
    memset(prev_lambda, 0, sizeof(prev_lambda));
    lambda[0] = 1;
    prev_lambda[0] = 1;

    int l = 0;
    int i = -1;

    for (int r = 0; r < nroots; ++r) {
        uint8_t discrepancy = codec->syndromes[r];
        for (int j = 1; j <= l; ++j) {
            discrepancy ^= gf_mult(lambda[j], codec->syndromes[r - j]);
        }

        if (discrepancy != 0) {
            uint8_t tau[nroots + 1];
            memset(tau, 0, sizeof(tau));

            uint8_t d_inv = gf_div(1, discrepancy);

            for(int k=0; k<=nroots; ++k) {
                tau[k] = lambda[k] ^ gf_mult(discrepancy, prev_lambda[k]);
            }

            if (l < (r - i)) {
                int tmp_l = r - i;
                i = r - l;
                l = tmp_l;

                for(int k=0; k<=nroots; ++k) {
                    prev_lambda[k] = gf_mult(lambda[k], d_inv);
                }
            }

            memcpy(lambda, tau, sizeof(lambda));
        }

        // Shift prev_lambda
        memmove(&prev_lambda[1], prev_lambda, nroots * sizeof(uint8_t));
        prev_lambda[0] = 0;
    }

    int deg_lambda = 0;
    for (int j = 0; j <= nroots; ++j) {
        if (lambda[j] != 0) {
            deg_lambda = j;
        }
    }

    if (deg_lambda > RS_T) {
        printf("Error: deg_lambda > RS_T (%d > %d)\n", deg_lambda, RS_T);
        return -1;
    }
    printf("Lambda (deg %d): ", deg_lambda);
    for(int i=0; i<=deg_lambda; i++) printf("%02X ", lambda[i]);
    printf("\n");

    // Chien Search
    int error_locs[RS_T];
    int error_count = 0;
    for (int j = 0; j < RS_N; ++j) {
        uint8_t q = 1;
        for (int k = 1; k <= deg_lambda; ++k) {
            q ^= gf_mult(lambda[k], gf_exp[(j * k) % 255]);
        }
        if (q == 0) {
            if (error_count < RS_T) {
                error_locs[error_count] = 255 - j;
            }
            error_count++;
        }
    }

    printf("Found %d roots. Error locations: ", error_count);
    for(int i=0; i<error_count; i++) printf("%d ", error_locs[i]);
    printf("\n");
    if (error_count != deg_lambda) {
        printf("Error: error_count != deg_lambda (%d != %d)\n", error_count, deg_lambda);
        return -1;
    }

    // Forney's Algorithm
    uint8_t omega[nroots + 1];
    memset(omega, 0, sizeof(omega));
    for (int j = 0; j < deg_lambda; ++j) {
        uint8_t term = 0;
        for (int k = 0; k <= j; ++k) {
            term ^= gf_mult(lambda[k], codec->syndromes[j - k]);
        }
        omega[j] = term;
    }

    for (int j = 0; j < error_count; ++j) {
        int loc = error_locs[j];
        if (loc >= RS_N) continue;

        uint8_t root = gf_exp[255 - loc];

        uint8_t omega_val = 0;
        for (int k = 0; k < deg_lambda; ++k) {
            omega_val ^= gf_mult(omega[k], gf_exp[(gf_log[root] * k) % 255]);
        }

        uint8_t lambda_prime_val = 0;
        for (int k = 1; k <= deg_lambda; k += 2) {
            lambda_prime_val ^= gf_mult(lambda[k], gf_exp[(gf_log[root] * (k - 1)) % 255]);
        }

        uint8_t error_val = gf_div(omega_val, lambda_prime_val);
        codeword[loc] ^= error_val;
    }

    return error_count;
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
    sys->rs_codec = rs_codec_init();

    sys->fec_buffer_size = 10000;
    sys->fec_buffer = calloc(sys->fec_buffer_size, sizeof(uint8_t));

    // Inizializza tabelle Galois e Reed-Solomon
    init_galois_field();
    rs_generate_polynomial();

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
        rs_codec_free(sys->rs_codec);
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
        rs_encode(sys->rs_codec, rs_input, rs_codeword);

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
        int rs_result = rs_decode(sys->rs_codec, codeword);

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
    printf("====== Inizio Test Reed-Solomon ======\n");

    // Inizializza il sistema (questo chiama init_galois_field, etc.)
    robust_dpsk_system_t *sys = system_init(1);
    if (!sys) {
        printf("Errore: impossibile inizializzare il sistema.\n");
        return 1;
    }

    // 1. Prepara i dati originali
    uint8_t original_data[RS_K];
    for (int i = 0; i < RS_K; i++) {
        original_data[i] = i % 251; // Un pattern di dati
    }

    // 2. Esegui l'encoding Reed-Solomon
    uint8_t codeword[RS_N];
    rs_encode(sys->rs_codec, original_data, codeword);
    printf("Dati originali encodati.\n");

    // 3. Introduci degli errori nel codeword
    uint8_t corrupted_codeword[RS_N];
    memcpy(corrupted_codeword, codeword, RS_N);

    int num_errors = 8;
    printf("Introduzione di %d errori nel codeword...\n", num_errors);
    for (int i = 0; i < num_errors; i++) {
        int loc = (i * 30) % RS_N; // Posizioni degli errori sparse
        uint8_t error_val = (i + 1) * 17;
        corrupted_codeword[loc] ^= error_val;
        printf("  - Errore in posizione %d, valore xor %02X\n", loc, error_val);
    }

    // 4. Esegui il decoding
    int corrected_errors = rs_decode(sys->rs_codec, corrupted_codeword);

    if (corrected_errors < 0) {
        printf("Errore: il decoder ha fallito (errore non correggibile).\n");
    } else {
        printf("Decoder ha corretto %d errori.\n", corrected_errors);
    }

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
        printf("I dati sono stati ripristinati correttamente dopo la correzione degli errori.\n");
    } else {
        printf("\n====== Test Reed-Solomon FALLITO! ======\n");
    }

    system_free(sys);
    return success ? 0 : 1;
}
