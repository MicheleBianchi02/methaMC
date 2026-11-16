#include <string.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <time.h>
#include <math.h>
#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif
#include "./pcg_basic.h"


// TODO: Cose da aggiungere:
// Aggiungere voltage offset (altrimenti problema quando t=0 e noise negativo, ADC non misura)
// Aggiungere la possibilità di salvare anche il numero di particelle in ogni
// integrazione.
// Rolling sha256 della sezione dati + index.json (come fa l'applicazione)


#define SAVE_FILE_NAME "./simulation_results.csv"
#define RAMP_FILE_NAME_CH1 "simulation_ramp_ch1.bin"
#define RAMP_FILE_NAME_CH2 "simulation_ramp_ch2.bin"

#define MAX_NUM_BIN_SPEC 10000   // Max number of bins for the idx_array of spectrum function
#define MAX_NUM_BIN_RATE 100000  // Max number of bins for the idx_array of rate function
// Macros used to remove the compiler warning and automatically check that the input 
// parameters match the expected one.
#define XSCANF(expected_matches, fmt, ...) do { \
    if (scanf(fmt, __VA_ARGS__) != (expected_matches)) { \
        fprintf(stderr, "Fatal: Input matching failed for format '%s'\n", fmt); \
        exit(EXIT_FAILURE); \
    } \
} while(0)

// ============ Data structures ===============
typedef struct DataPoint {
    double x;
    double y;
} DataPoint;

typedef struct SpectrumBuffer {
    DataPoint *data;
    size_t *idx_array;  // Array used for faster access
    size_t sz;          // Total allocated space 
    size_t ptr;         // Current number of elements stored (the pointer position)
} SpectrumBuffer;


typedef struct RateBuffer {
    DataPoint *data;
    size_t *idx_array;     // Array used for faster access
    size_t sz;             // Total allocated space 
    size_t ptr;            // Current number of elements stored (the pointer position)
    double inv_bin_width;  // inverse of widht of each bin for the idx_array
    double period;         // Period of the given rate function
    double t_min;          // Time offset (first point of the given rate function)
    double rate_gain;      // Multiplicative factor for the rate
    double max_rate;       // Maximum value of the rate
} RateBuffer;

typedef struct SimulationParameters {
    SpectrumBuffer *spectrum; 
    RateBuffer *rate_buff;
    pcg32_random_t *rng;       // Random number generator
    double t_int;              // integration time [s]
    double t_rst;              // reset time [s]
    double rate;               // expressed in particle / s
    double ramp_dt;            // time between two ramp samples [s]
    double int_capacitor;      // Integrating capacitor [pF]
    double ADC_min_range_V;    // minimum value the adc input can reach (in V)
    double ADC_max_range_V;    // maximum value the adc input can reach (in V)
    double electronic_noise;   // electronic noise (from user input) expressed in V
                               // sampled from a gaussian distribution
    double sinusoid_omega;     // = 2*pi/period
    double sinusoid_aplitude;  // expressed in particle / s
    double leakage_current;    // expressed in pA 
    double relative_gain;      // Gain relative to the spectrum data (simple multiplicative 
                               // factor accounting for different gain between the spectrum
                               // and the simulation)
    int ADC_bit;               // ADC number of bit
    int N_points;              // Number of simulated integrations
    int output_in_charge;      // If == 1 the output is in charge, otherwise in ADC channels
    int add_sinusoid;          // If == 1 a sinusoid is present
    int add_noise;             // If == 1 the noise of the ADC and input noise are included
    int sample_ramp;           // If == 1 the integration ramp is continuously sampled
    int custom_rate;           // If == 1 a custom rate function is used to sample the rate
    const char *ramp_dir;      // Output directory of the ramp .bin files (NULL = CWD)
} SimulationParameters;

typedef struct Simulation {
    SimulationParameters *params;
    double *result_ch1;       // buffer containing the results for channel 1
    double *result_ch2;       // buffer containing the results for channel 2
    uint32_t *ramp_out_ch1;   // full-run quantized ramp for channel 1 (API
                              // output; NULL = keep only the files)
    uint32_t *ramp_out_ch2;   // full-run quantized ramp for channel 2
    uint64_t ramp_samples_ch1;// Exact number of stored ramp samples for channel 1
    uint64_t ramp_samples_ch2;// Exact number of stored ramp samples for channel 2
} Simulation;


/* Ramp file header, byte-compatible with the application format "AQR1 v1"
 * (services/ramp_storage.py).  Fixed 4096-byte header followed by raw
 * little-endian uint32 ADC samples; the offsets below are asserted to match
 * the Python struct layout exactly (same asserts as the application tests).
 * Endianness: the fields and the samples are written as-is, so a
 * little-endian host is assumed (x86-64, ARM LE). */
#pragma pack(push, 1)
typedef struct RampFileHeader {
    char     magic[4];             // "AQR1"
    uint16_t header_version;       // 1
    uint16_t header_size;          // 4096
    uint32_t flags;                // 0
    uint32_t sample_dtype;         // 0 = little-endian uint32
    uint8_t  channel_index;        // 0 = ch1, 1 = ch2
    uint8_t  state;                // 0 running, 1 finalized, 2 aborted
    uint16_t timing_mode;          // 1 = ramp
    uint32_t adc_bits;
    uint64_t samples_per_window;   // N = round(t_int/ramp_dt) + 1
    uint64_t samples_written;      // watermark, patched at finalize
    uint64_t data_offset;          // 4096
    double   ramp_dt_s;
    double   period_s;             // t_int + t_rst
    double   integration_time_s;
    double   reset_time_s;
    double   sample_rate_hz;       // 1 / ramp_dt
    double   range_min_v;
    double   range_max_v;
    double   start_epoch_s;
    uint8_t  data_sha256[32];      // zeros: rolling hash not implemented yet
    uint8_t  config_sha256[32];    // zeros
    char     measurement_uuid[32]; // 32 ASCII hex chars, NOT NUL terminated
    char     start_iso[64];
    char     channel_label[32];    // "ch1" / "ch2"
    char     driver_code[32];      // "simulation"
    char     firmware_version[32];
    char     app_version[16];
    uint8_t  reserved[64];
} RampFileHeader;
#pragma pack(pop)

#define RAMP_HEADER_SIZE 4096
#define RAMP_STATE_RUNNING  0
#define RAMP_STATE_FINALIZED 1

_Static_assert(sizeof(RampFileHeader) == 448, "AQR1 header size drifted from services/ramp_storage.py");
_Static_assert(offsetof(RampFileHeader, samples_per_window) == 0x18, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, samples_written) == 0x20, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, data_offset) == 0x28, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, ramp_dt_s) == 0x30, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, data_sha256) == 0x70, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, measurement_uuid) == 0xB0, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, start_iso) == 0xD0, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, app_version) == 0x170, "AQR1 layout mismatch");
_Static_assert(offsetof(RampFileHeader, reserved) == 0x180, "AQR1 layout mismatch");


// ====================== Helper functions ================

/* Wrapper for malloc. Check if the memory allocation was successull */
void *xmalloc(size_t size) {
    void *ptr = malloc(size);
    
    // If malloc returns NULL, the allocation failed
    if (ptr == NULL) {
        fprintf(stderr, "Fatal error: memory allocation failed for %zu bytes.\n", size);
        exit(EXIT_FAILURE);
    }
    
    return ptr;
}



/* Wrapper for realloc. Check if the memory allocation was successull */
void *xrealloc(void *buff, size_t size) {
    void *ptr = realloc(buff, size);

    if (ptr == NULL) {
        fprintf(stderr, "Fatal: Memory reallocation failed.\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

// ================ AQR1 ramp file helpers ==============

/* Best-effort creation of the ramp output directory (ignored if it exists) */
static void EnsureRampDirectory(const char *path) {
    if (path == NULL || path[0] == '\0') return;
#ifdef _WIN32
    _mkdir(path);
#else
    mkdir(path, 0777);
#endif
}

/* Cheap deterministic uuid look-alike: 32 hex chars, unique per run */
static uint64_t SplitMix64(uint64_t x) {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

static void GenerateMeasurementUuid(char out[32], uint64_t seed) {
    char buf[33];
    uint64_t a = SplitMix64(seed);
    uint64_t b = SplitMix64(a);
    snprintf(buf, sizeof(buf), "%016" PRIx64 "%016" PRIx64, a, b);
    memcpy(out, buf, 32);
}

/* Local time as ISO 8601 with offset, e.g. 2026-10-20T14:23:45+02:00
 * (the application writes the same shape).  tm_gmtoff is POSIX only. */
static void FormatStartIso(char out[64]) {
    time_t now = time(NULL);
    struct tm tm_local;
    long offset_s = 0;
#ifdef _WIN32
    localtime_s(&tm_local, &now);
#else
    localtime_r(&now, &tm_local);
    offset_s = tm_local.tm_gmtoff;
#endif
    char sign = (offset_s < 0) ? '-' : '+';
    long abs_s = (offset_s < 0) ? -offset_s : offset_s;
    snprintf(out, 64, "%04d-%02d-%02dT%02d:%02d:%02d%c%02ld:%02ld",
             tm_local.tm_year + 1900, tm_local.tm_mon + 1, tm_local.tm_mday,
             tm_local.tm_hour, tm_local.tm_min, tm_local.tm_sec,
             sign, abs_s / 3600, (abs_s % 3600) / 60);
}

static void RampHeaderInit(RampFileHeader *h, const SimulationParameters *p,
                           int channel_index, uint64_t samples_per_window,
                           double start_epoch_s, const char *start_iso,
                           const char uuid[32]) {
    memset(h, 0, sizeof(*h));
    memcpy(h->magic, "AQR1", 4);
    h->header_version = 1;
    h->header_size = RAMP_HEADER_SIZE;
    h->sample_dtype = 0;
    h->channel_index = (uint8_t)channel_index;
    h->state = RAMP_STATE_RUNNING;
    h->timing_mode = 1;   // ramp
    h->adc_bits = (uint32_t)p->ADC_bit;
    h->samples_per_window = samples_per_window;
    h->samples_written = 0;
    h->data_offset = RAMP_HEADER_SIZE;
    h->ramp_dt_s = p->ramp_dt;
    h->period_s = p->t_int + p->t_rst;
    h->integration_time_s = p->t_int;
    h->reset_time_s = p->t_rst;
    h->sample_rate_hz = (p->ramp_dt > 0.0) ? 1.0 / p->ramp_dt : 0.0;
    h->range_min_v = p->ADC_min_range_V;
    h->range_max_v = p->ADC_max_range_V;
    h->start_epoch_s = start_epoch_s;
    memcpy(h->measurement_uuid, uuid, 32);
    snprintf(h->start_iso, sizeof(h->start_iso), "%s", start_iso);
    snprintf(h->channel_label, sizeof(h->channel_label), "ch%d", channel_index + 1);
    snprintf(h->driver_code, sizeof(h->driver_code), "simulation");
    snprintf(h->app_version, sizeof(h->app_version), "spectrum-mc");
}

/* (Re)write the whole 4096-byte header block at offset 0 */
static void WriteRampHeader(FILE *f, const RampFileHeader *h) {
    static const uint8_t pad[RAMP_HEADER_SIZE - sizeof(RampFileHeader)] = {0};
    if (fseek(f, 0, SEEK_SET) != 0 ||
        fwrite(h, sizeof(*h), 1, f) != 1 ||
        fwrite(pad, sizeof(pad), 1, f) != 1) {
        fprintf(stderr, "Fatal: could not write the ramp file header.\n");
        exit(EXIT_FAILURE);
    }
}

/* Validate the ramp configuration and return the single source of truth for
 * the window size: N = round(t_int/ramp_dt) + 1, the same geometry the
 * application chunks the stream with (it refuses N < 2). */
static uint64_t RampSamplesPerWindow(const SimulationParameters *p) {
    if (p->ramp_dt <= 0.0) {
        fprintf(stderr, "Fatal: sample_ramp requires ramp_dt > 0.\n");
        exit(EXIT_FAILURE);
    }
    if (p->ramp_dt > p->t_int) {
        fprintf(stderr, "Fatal: ramp_dt > t_int: a window needs at least 2 "
                        "samples (application requirement).\n");
        exit(EXIT_FAILURE);
    }
    if (p->ADC_bit < 1 || p->ADC_bit > 32) {
        fprintf(stderr, "Fatal: ramp files need 1 <= ADC_bit <= 32.\n");
        exit(EXIT_FAILURE);
    }
    if (p->int_capacitor <= 0.0) {
        fprintf(stderr, "Fatal: ramp files need int_capacitor > 0 "
                        "(charge must be convertible to volts).\n");
        exit(EXIT_FAILURE);
    }
    if (!(p->ADC_max_range_V > p->ADC_min_range_V)) {
        fprintf(stderr, "Fatal: ramp files need ADC_max_range_V > ADC_min_range_V.\n");
        exit(EXIT_FAILURE);
    }
    return (uint64_t)((p->t_int / p->ramp_dt) + 0.5) + 1;
}

// ================== Spectrum buffer functions ===========

int InitSpectrumBuffer(SpectrumBuffer *buff, size_t initial_capacity) {
    buff->data = xmalloc(initial_capacity * sizeof(DataPoint));

    buff->sz = initial_capacity;
    buff->ptr = 0;
    return 1;
}

void AppendSpectrumPoint(SpectrumBuffer *buff, DataPoint point) {
    if (buff->ptr >= buff->sz) {
        buff->sz *= 2;
        buff->data = (DataPoint *)xrealloc(buff->data, buff->sz * sizeof(DataPoint));
    }

    buff->data[buff->ptr] = point;
    buff->ptr++;
}

void FreeSpectrumBuffer(SpectrumBuffer *buff) {
    free(buff->data);
    if (buff->idx_array) free(buff->idx_array);
    buff->sz = 0;
    buff->ptr = 0;
}

// Calculate the Comulative Distribution Function of the charge spectrum. 
void CalculateCDF(SpectrumBuffer *buffer) {

    if (buffer->ptr == 0) return;

    for (size_t i=0; i<buffer->ptr - 1; i++) {
        buffer->data[i+1].y += buffer->data[i].y;
    }

    // Normalize to have the frequency spectrum
    double total_sum = buffer->data[buffer->ptr - 1].y;

    if (total_sum == 0.0) return;

    // Multiplication is much faster for the CPU than division
    double inv_sum = 1.0 / total_sum;

    for (size_t i = 0; i < buffer->ptr; i++) {
        buffer->data[i].y *= inv_sum;
    }
}


void ComputeSpectrumIdxArray(SpectrumBuffer *buff) {
    size_t num_bins = MAX_NUM_BIN_SPEC;
    buff->idx_array = xmalloc(num_bins * sizeof(size_t));
    double bin_width = 1.0 / num_bins;

    size_t idx = 0;
    for (size_t b = 0; b < num_bins; b++) {
        double r_val = b * bin_width;
        // Find the first index where the CDF (y) is greater than or equal to r_val
        while (idx < buff->ptr - 1 && buff->data[idx].y < r_val) {
            idx++;
        }
        // Step back by one to ensure we don't overshoot the lower bound during interpolation
        buff->idx_array[b] = (idx > 0) ? idx - 1 : 0; 
    }
}



// ================= Rate buffer functions =============

int InitRateBuffer(RateBuffer *buff, size_t initial_capacity) {
    buff->data = xmalloc(initial_capacity * sizeof(DataPoint));

    buff->sz = initial_capacity;
    buff->ptr = 0;
    return 1;
}

void AppendRatePoint(RateBuffer *buff, DataPoint point) {
    if (buff->ptr >= buff->sz) {
        buff->sz *= 2;
        buff->data = (DataPoint *)realloc(buff->data, buff->sz * sizeof(DataPoint));
    }

    buff->data[buff->ptr] = point;
    buff->ptr++;
}

void FreeRateBuffer(RateBuffer *buff) {
    free(buff->data);
    buff->sz = 0;
    buff->ptr = 0;
}


void ComputeRateIdxArray(RateBuffer *buff) {
    size_t num_bins = MAX_NUM_BIN_RATE;
    buff->idx_array = xmalloc(num_bins * sizeof(size_t));
    
    buff->inv_bin_width = 1 / (buff->period / num_bins);

    size_t idx = 0;
    for (size_t b = 0; b < num_bins; b++) {
        double t = (double)b / num_bins * buff->period;

        // Find the index i where t_i <= t
        // The + 1 is needed because we are incrementing inside the while loop. Without 
        // the +1 we would increment the index at the required t_i.
        while (idx < buff->ptr - 1 && buff->data[idx + 1].x <= t) {
            idx++;
        }

        buff->idx_array[b] = idx;
    }
}




// ================= File functions ===================

// Import spectrum from CSV file
void ReadSpectrumCSV(const char *filename, SpectrumBuffer *buffer) {

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening the file");
        exit(EXIT_FAILURE);
    }

    char line[256];
    
    // skip header
    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error: File is empty or cannot be read.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    DataPoint point;
    while (fgets(line, sizeof(line), file) != NULL) {
        if (sscanf(line, "%lf,%lf", &point.x, &point.y) == 2) {
            AppendSpectrumPoint(buffer, point);
        } else {
            fprintf(stderr, "Warning: Skipped malformed line -> %s", line);
        }
    }

    fclose(file);
}



// Import Rate function from CSV file. It is assumed that the x axis (first column) is 
// sorted and no duplicate exist.
void ReadRateCSV(const char *filename, RateBuffer *buffer) {

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening the file");
        exit(EXIT_FAILURE);
    }

    char line[256];
    
    // skip header
    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error: File is empty or cannot be read.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    DataPoint point;

    // First point to define t_min
    int ii = 0;
    while (fgets(line, sizeof(line), file) != NULL && ii == 0) {
        if (sscanf(line, "%lf,%lf", &point.x, &point.y) == 2) {
            point.y *= buffer->rate_gain;
            AppendRatePoint(buffer, point);
            buffer->t_min = point.x;
            ii = 1;
        } else {
            fprintf(stderr, "Warning: Skipped malformed line -> %s", line);
        }
    }

    double max_rate = point.y;
    while (fgets(line, sizeof(line), file) != NULL) {
        if (sscanf(line, "%lf,%lf", &point.x, &point.y) == 2) {
            point.y *= buffer->rate_gain;
            if (max_rate < point.y) max_rate = point.y;
            AppendRatePoint(buffer, point);
        } else {
            fprintf(stderr, "Warning: Skipped malformed line -> %s", line);
        }
    }

    // Last point to define the period
    buffer->period = point.x - buffer->t_min;
    buffer->max_rate = max_rate;

    if (max_rate <= 0) {
        fprintf(stderr, "Fatal: The provided custom rate function is faulty.\n");
        exit(EXIT_FAILURE);
    }

    fclose(file);
}


// Save result of simulation into a CSV file
void SaveDataToCSV(const char *filename,const char *header, Simulation *simul) {

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening the file");
        return;
    }

    fprintf(file, "%s\n", header);

    int n_points = simul->params->N_points;
    double t_step = simul->params->t_int + simul->params->t_rst;
    int is_charge = (simul->params->output_in_charge == 1);

    for (int i = 0; i < n_points; i++) {
        double timestamp = (i + 1) * t_step; 
        
        if (is_charge) {
            fprintf(file, "%f, %lf, %lf\n", timestamp, simul->result_ch1[i], 
                                                    simul->result_ch2[i]);
        } else {

            fprintf(file, "%f, %d, %d\n", timestamp, (int)simul->result_ch1[i], 
                                                (int)simul->result_ch2[i]);
        }
    }
    fclose(file);
}




// ========= Monte Carlo simulation ==========

// Generate a random double between 0.0 (inclusive) and 1.0 (exclusive)
static inline double GetUniformRandom(pcg32_random_t *rng) {
    // We multiply by 2^-32 to get a perfect [0.0, 1.0) distribution.
    // 0x1.0p-32 is the C99 standard way to write exactly 1 / 4294967296.
    double random_double = pcg32_random_r(rng) * 0x1.0p-32;
    return random_double;
}


double SampleData(SpectrumBuffer *spectrum_buffer, pcg32_random_t *rng, double rel_gain) {

    if (spectrum_buffer->ptr == 0) return 0.0;

    double r = GetUniformRandom(rng);
    
    size_t num_bin = MAX_NUM_BIN_SPEC;
    size_t bin = (size_t)(r * num_bin); 
    if (bin >= num_bin) bin = num_bin - 1;       // Safety catch
    
    size_t idx = spectrum_buffer->idx_array[bin];

    // Tiny linear search to find the exact segment (usually 0 to 2 iterations)
    while (idx < spectrum_buffer->ptr - 1 && spectrum_buffer->data[idx + 1].y < r) {
        idx++;
    }

    // Interpolate
    double x0 = spectrum_buffer->data[idx].x;
    double y0 = spectrum_buffer->data[idx].y;
    double x1 = spectrum_buffer->data[idx+1].x;
    double y1 = spectrum_buffer->data[idx+1].y;

    double y_diff = y1 - y0;
    if (y_diff < 1e-12) return x1 * rel_gain;

    return (x0 + (r - y0) * (x1 - x0) / y_diff) * rel_gain;
}



static inline double GetCustomRate(SimulationParameters *params, double t_wrapped ) {
    RateBuffer *rb = params->rate_buff;

    // this find the nominal bin: bin = t / period * num_bin
    // Used inverse of the bin_width because multiplication are faster than divisions
    size_t bin = t_wrapped * rb->inv_bin_width;
    if (bin >= MAX_NUM_BIN_RATE) bin = MAX_NUM_BIN_RATE - 1;

    size_t idx = rb->idx_array[bin];

    while (idx < rb->ptr - 1 && rb->data[idx + 1].x <= t_wrapped) {
        idx++;
    }

    // Safety catch to prevent out-of-bounds array access during interpolation
    if (idx >= rb->ptr - 1) {
        return rb->data[rb->ptr - 1].y;
    }

    // At this point we have that t_wrapped is inside rate_buff->data[idx] and rate_buff->data[idx+1]
    // We can do a linear interpolation to find the estimated rate
    return rb->data[idx].y + (rb->data[idx+1].y - rb->data[idx].y) * 
        (t_wrapped - rb->data[idx].x) / (rb->data[idx+1].x - rb->data[idx].x);

}


// Extract from an exponential distribution the time between two particle interactions
// The rate parameter is expressed in particles (or interactions) per second
static inline double GetArrivalTime(SimulationParameters *params, double current_time) {
    
    // Handle the stationary beam case
    if (params->add_sinusoid == 0 && params->custom_rate == 0) {
        double u = GetUniformRandom(params->rng);
        // If u is exactly 0.0, add the smallest possible non-zero increment
        if (u == 0.0) u = 0x1.0p-32;
        return -log(u) / params->rate;
    } else {
     
        // If the rate r = r0 + R * sin(omega * t), we have a
        // Non-Homogeneous Poisson Process. In this case we use the 
        // Thinning Algorithm (Lewis-Shedler Method)
        double max_rate;
        double t = current_time; // 't' acts as our running timeline
        double t_wrapped = 0;
        if (params->add_sinusoid == 1) {
            max_rate = params->rate + params->sinusoid_aplitude;
        } else {
            max_rate = params->rate_buff->max_rate;

            /* current_time = t: 
             * Considering that the rate function is periodic with period T we have that 
             * t = a*T + b with a = int and b < T. 
             * The fmod(t, T) function return the value t - n*T where n = floor(t/T). So, 
             * n = a and fmod(t, T) = b as expected. If t < T a = n = 0. 
             * So fmod can be used to wrap around the time value. */
            t_wrapped = fmod(current_time, params->rate_buff->period);

            // fmod can return negative numbers in C if t is negative. 
            // This forces it into the positive [0, T) range:
            if (t_wrapped < 0.0) t_wrapped += params->rate_buff->period;
        }
        
        while(1) {
            double u = GetUniformRandom(params->rng);
            if (u == 0.0) u = 0x1.0p-32;

            // Calculate the time step to the next proposed particle
            double dt = -log(u) / max_rate;
            
            // Advance the running clock. Time passes whether we accept or reject
            t += dt;


            double actual_rate;
            if (params->add_sinusoid == 1) {
                actual_rate =  params->rate + params->sinusoid_aplitude * 
                                            sin(params->sinusoid_omega * t);
            } else {
                t_wrapped += dt;
                while (t_wrapped >= params->rate_buff->period) {
                    t_wrapped -= params->rate_buff->period;
                }
                actual_rate = GetCustomRate(params, t_wrapped);
            }
            
            double u2 = GetUniformRandom(params->rng);
            if (u2 <= (actual_rate / max_rate)) {
                break;
            }

        }

        return t - current_time;
    }
}


/* Calculate the value at the end of the integration while sampling the ramp */
static inline double CalculateIntegration(SimulationParameters *params, 
                                   double *current_time, 
                                   double next_sample, 
                                   double end_of_integration, 
                                   double end_of_cycle,
                                   double relative_gain,
                                   double *ramp_buffer,
                                   size_t *ramp_idx) {
    

    double integrated_yval = 0.0;
    const double ramp_dt = params->ramp_dt;
    // next_sample is initially passed as the start of the cycle. Save it to measure elapsed time.
    const double cycle_start_time = next_sample;


    if (params->sample_ramp) {
        // Using integer comparison because comparing floatig numbers is prone to 
        // small comulated errors

        // Adding 0.5 before casting to size_t acts as a safe rounding mechanism
        size_t expected_samples = (size_t)((params->t_int / ramp_dt) + 0.5) + 1;
        size_t samples_taken = 0;

        /* current_time_ch1 is incremented beyond the beginning of the next 
         * integration time. Since when we sample the ramp we do it always after the 
         * same amount of time we must consider the case where current_time_ch1 is 
         * after the moment the ramp has to be sampled. */
        while (samples_taken < expected_samples && next_sample <= *current_time) {
            double current_leakage = params->leakage_current * (next_sample - cycle_start_time);
            ramp_buffer[(*ramp_idx)++] = integrated_yval + current_leakage;
            samples_taken++;
            // Multiplication prevents floating point drift from accumulating!
            next_sample = cycle_start_time + (samples_taken * ramp_dt); 
        }


        // Keep accumulating particles as long as they hit before the window closes
        while (*current_time < end_of_integration) {

            while (samples_taken < expected_samples && next_sample <= *current_time) {
                double current_leakage = params->leakage_current * (next_sample - cycle_start_time);
                ramp_buffer[(*ramp_idx)++] = integrated_yval + current_leakage;
                samples_taken++;
                next_sample = cycle_start_time + (samples_taken * ramp_dt);
            }
            
            integrated_yval += SampleData(params->spectrum, params->rng, relative_gain);
            
            // Generate the next particle. Time keeps increasing 
            *current_time += GetArrivalTime(params, *current_time);
        }

        /* Enter here if the next particle will arrive after end_of_interval but some 
         * ramp sample still need to be taken. */
        while (samples_taken < expected_samples) {
            double current_leakage = params->leakage_current * (next_sample - cycle_start_time);
            ramp_buffer[(*ramp_idx)++] = integrated_yval + current_leakage;
            samples_taken++;
            next_sample = cycle_start_time + (samples_taken * ramp_dt);
        }

    } else {
        
        while (*current_time < end_of_integration) {
            integrated_yval += SampleData(params->spectrum, params->rng, relative_gain);
            *current_time += GetArrivalTime(params, *current_time);
        }
    }

    while (*current_time < end_of_cycle) {
        *current_time += GetArrivalTime(params, *current_time);
    }



    // Once current_time overshoots the window, we save the result.
    // We do not reset current_time. That "overshoot" particle becomes 
    // the first candidate for the next integration window.
    return integrated_yval + params->leakage_current * params->t_int;
}






// Sample elctronic noise from gaussian distribution. Same unit of measure of 
// params->electronic_noise
double SampleNoise(SimulationParameters *params) {
    double u1 = GetUniformRandom(params->rng);
    double u2 = GetUniformRandom(params->rng);
    if (u1 == 0.0) u1 = 0x1.0p-32;
    if (u2 == 0.0) u2 = 0x1.0p-32;
    
    // The Box-Muller Transform is used to sample from a uniform gaussian distribution. 
    // since the mean is zero, given that the rms value is the standard deviation, the 
    // noise distribution is simply Z * sigma
    double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return z0 * params->electronic_noise;
}



// Quantize one cycle's ramp buffer (pC) into uint32 ADC channels.
//
// The application converts back with V = V_min + counts/(2^bits - 1) * span
// (RampChannel.to_volts), so the forward conversion MUST use 2^bits - 1 as
// the full-scale code - using 2^bits (the old behaviour) left a systematic
// (2^bits - 1)/2^bits gain mismatch of ~1 LSB on every sample.
//
// The ramp files always store channels: the AQR1 format has a single dtype
// code, exactly like the real ADC stream.  Charge output stays in the
// results CSV / API arrays only.
void ConvertRampBufferToChannels(const double * restrict ramp_buffer, size_t count,
                                        SimulationParameters *params,
                                        uint32_t * restrict out) {
    const double num_bin = (double)((1ULL << params->ADC_bit) - 1);
    const double adc_range = params->ADC_max_range_V - params->ADC_min_range_V;
    const double adc_max_ch = num_bin;
    const double adc_min_ch = 0.0;

    const double conversion_pC_ch = num_bin / (params->int_capacitor * adc_range);
    const double conversion_V_ch = num_bin / adc_range;
    const double ch_offset = (0.0 - params->ADC_min_range_V) * conversion_V_ch;

    if (params->add_noise == 0) {
        for (size_t i = 0; i < count; i++) {
            double val = round((ramp_buffer[i] * conversion_pC_ch) + ch_offset);
            out[i] = (uint32_t)fmax(adc_min_ch, fmin(val, adc_max_ch));
        }
    } else {
        for (size_t i = 0; i < count; i++) {
            double noise = SampleNoise(params);
            double val = round((ramp_buffer[i] * conversion_pC_ch) + ch_offset + (noise * conversion_V_ch));
            out[i] = (uint32_t)fmax(adc_min_ch, fmin(val, adc_max_ch));
        }
    }
}



void ConvertIntoChannel(Simulation *simul) {
    const int N_points = simul->params->N_points;
    // Same 2^bits - 1 full-scale convention as ConvertRampBufferToChannels:
    // the CSV results and the ramp files must agree bit for bit.
    const double num_bin = (double)((1ULL << simul->params->ADC_bit) - 1);
    
    double adc_range = simul->params->ADC_max_range_V - simul->params->ADC_min_range_V;
    
    const double adc_max_ch = num_bin; 
    const double adc_min_ch = 0.0;
    
    // Calculate slopes
    const double conversion_pC_ch = num_bin / (simul->params->int_capacitor * adc_range);
    const double conversion_V_ch = num_bin / adc_range;
    
    // Calculate offset if the minimum voltage is not 0V.
    const double ch_offset = (0.0 - simul->params->ADC_min_range_V) * conversion_V_ch;

    // restrict is used to enable SIMD (vectorization) to increase performances
    double * restrict ch1 = simul->result_ch1;
    double * restrict ch2 = simul->result_ch2;

    /* --- CONVERT FINAL INTEGRATION RESULTS --- */
    if (simul->params->add_noise == 0) {
        for (int i = 0; i < N_points; i++) {
            double val1 = round((ch1[i] * conversion_pC_ch) + ch_offset);
            double val2 = round((ch2[i] * conversion_pC_ch) + ch_offset);

            ch1[i] = fmax(adc_min_ch, fmin(val1, adc_max_ch));
            ch2[i] = fmax(adc_min_ch, fmin(val2, adc_max_ch));
        }
    } else {
        for (int i = 0; i < N_points; i++) {
            double noise1 = SampleNoise(simul->params);
            double noise2 = SampleNoise(simul->params);

            double val1 = round((ch1[i] * conversion_pC_ch) + ch_offset + (noise1 * conversion_V_ch));
            double val2 = round((ch2[i] * conversion_pC_ch) + ch_offset + (noise2 * conversion_V_ch));
            
            ch1[i] = fmax(adc_min_ch, fmin(val1, adc_max_ch));
            ch2[i] = fmax(adc_min_ch, fmin(val2, adc_max_ch));
        }
    }

}





void RunSimulation(Simulation *simul, const int verbose) {
    
    // Calculate the absolute maximum time the simulation will reach
    double total_sim_time = simul->params->N_points * (simul->params->t_int + simul->params->t_rst);
    
    // Find the maximum possible rate
    double peak_rate = simul->params->rate;
    if (simul->params->add_sinusoid == 1) {
        peak_rate += simul->params->sinusoid_aplitude;
    } else if (simul->params->custom_rate == 1) {
        peak_rate = simul->params->rate_buff->max_rate;
    }
    double estimated_particles = total_sim_time * peak_rate;

    if (verbose) {
        // Calculate the maximum possible rate depending on the mode
        double peak_rate = simul->params->rate;
        if (simul->params->add_sinusoid == 1) {
            peak_rate += simul->params->sinusoid_aplitude;
        } else if (simul->params->custom_rate == 1) {
            peak_rate = simul->params->rate_buff->max_rate;
        }

        // Calculate total time and estimated particles

        printf("\n--- Simulation Setup ---\n");
        printf("Iterations: %d\n", simul->params->N_points);
        printf("Total simulated time: %.4f seconds\n", total_sim_time);
        printf("Peak particle rate: %.0f p/s\n", peak_rate);
        // Using %.0f or %e to prevent messy decimal printing for massive numbers
        if (estimated_particles > 1e6) {
            printf("Estimated total particles: %.2e\n", estimated_particles);
        } else {
            printf("Estimated total particles: %.0f\n", estimated_particles);
        }
        printf("------------------------\n\n");

    }
    // The Danger Threshold (1 Trillion particles)
    if (estimated_particles > 1e12) {
        printf("\n==================== WARNING ====================\n");
        printf("Estimated total particles: %.2e\n", estimated_particles);
        printf("This simulation is approaching the physical precision limits\n");
        printf("of 64-bit computer memory. Time calculations (t += dt) may\n");
        printf("experience floating-point truncation, skewing your statistics.\n");
        printf("Consider breaking this into smaller N_points runs.\n");
        printf("=================================================\n\n");
    }



    int N = 0;
    double current_time_ch1 = 0.0; 
    double current_time_ch2 = 0.0; 
    double relative_gain = simul->params->relative_gain;


    FILE *f_ramp1 = NULL;
    FILE *f_ramp2 = NULL;
    RampFileHeader ramp_hdr1, ramp_hdr2;
    double *ramp_scratch1 = NULL;   // per-cycle pC ramps (reused every cycle)
    double *ramp_scratch2 = NULL;
    uint32_t *ramp_stage1 = NULL;   // per-cycle quantized channels
    uint32_t *ramp_stage2 = NULL;
    char ramp_path1[1024];
    char ramp_path2[1024];

    if (simul->params->sample_ramp) {

        uint64_t samples_per_window = RampSamplesPerWindow(simul->params);

        ramp_scratch1 = xmalloc(samples_per_window * sizeof(double));
        ramp_scratch2 = xmalloc(samples_per_window * sizeof(double));
        ramp_stage1 = xmalloc(samples_per_window * sizeof(uint32_t));
        ramp_stage2 = xmalloc(samples_per_window * sizeof(uint32_t));

        const char *ramp_dir = (simul->params->ramp_dir != NULL) ?
                               simul->params->ramp_dir : ".";
        EnsureRampDirectory(ramp_dir);
        snprintf(ramp_path1, sizeof(ramp_path1), "%s/%s", ramp_dir, RAMP_FILE_NAME_CH1);
        snprintf(ramp_path2, sizeof(ramp_path2), "%s/%s", ramp_dir, RAMP_FILE_NAME_CH2);

        f_ramp1 = fopen(ramp_path1, "wb");
        f_ramp2 = fopen(ramp_path2, "wb");
        if (f_ramp1 == NULL || f_ramp2 == NULL) {
            perror("Fatal: could not open the ramp output files");
            exit(EXIT_FAILURE);
        }
        setvbuf(f_ramp1, NULL, _IOFBF, 16 * 1024 * 1024);
        setvbuf(f_ramp2, NULL, _IOFBF, 16 * 1024 * 1024);

        // Both channels of one run share uuid and start time (the
        // application's index sidecar relies on that)
        char measurement_uuid[32];
        char start_iso[64];
        double start_epoch_s = (double)time(NULL);
        GenerateMeasurementUuid(measurement_uuid,
                                ((uint64_t)time(NULL) << 20) ^
                                (uint64_t)(uintptr_t)simul);
        FormatStartIso(start_iso);

        RampHeaderInit(&ramp_hdr1, simul->params, 0, samples_per_window,
                       start_epoch_s, start_iso, measurement_uuid);
        RampHeaderInit(&ramp_hdr2, simul->params, 1, samples_per_window,
                       start_epoch_s, start_iso, measurement_uuid);

        // Running headers (state = RUNNING, watermark = 0)
        WriteRampHeader(f_ramp1, &ramp_hdr1);
        WriteRampHeader(f_ramp2, &ramp_hdr2);

        if (verbose) {
            printf("Ramp files: %s, %s (%llu samples/window)\n",
                   ramp_path1, ramp_path2,
                   (unsigned long long)samples_per_window);
        }
    }

    uint64_t total_ramp_samples_ch1 = 0;
    uint64_t total_ramp_samples_ch2 = 0;

    current_time_ch1 += GetArrivalTime(simul->params, current_time_ch1);
    current_time_ch2 += GetArrivalTime(simul->params, current_time_ch2);



    const int update_interval = simul->params->N_points / 10;

    while(N < simul->params->N_points) {
        double cycle_start_time = (double)N * (simul->params->t_int + simul->params->t_rst);
        double end_of_integration = cycle_start_time + simul->params->t_int;
        double end_of_cycle = end_of_integration + simul->params->t_rst;

        // Reset buffer indices back to 0 for the new cycle
        size_t cycle_idx_ch1 = 0;
        size_t cycle_idx_ch2 = 0;

        simul->result_ch1[N] = CalculateIntegration(simul->params, &current_time_ch1, 
                                cycle_start_time, end_of_integration, end_of_cycle, 
                                relative_gain, ramp_scratch1, &cycle_idx_ch1);

        simul->result_ch2[N] = CalculateIntegration(simul->params, &current_time_ch2, 
                                cycle_start_time, end_of_integration, end_of_cycle, 
                                relative_gain, ramp_scratch2, &cycle_idx_ch2);

        // Convert and stream to disk immediately.  The conversion order
        // (ch1 then ch2, after both integrations) keeps the RNG stream
        // identical to the previous file format for a given seed.
        if (simul->params->sample_ramp) {
            ConvertRampBufferToChannels(ramp_scratch1, cycle_idx_ch1, simul->params, ramp_stage1);
            ConvertRampBufferToChannels(ramp_scratch2, cycle_idx_ch2, simul->params, ramp_stage2);

            fwrite(ramp_stage1, sizeof(uint32_t), cycle_idx_ch1, f_ramp1);
            fwrite(ramp_stage2, sizeof(uint32_t), cycle_idx_ch2, f_ramp2);

            if (simul->ramp_out_ch1 != NULL) {
                memcpy(simul->ramp_out_ch1 + total_ramp_samples_ch1,
                       ramp_stage1, cycle_idx_ch1 * sizeof(uint32_t));
            }
            if (simul->ramp_out_ch2 != NULL) {
                memcpy(simul->ramp_out_ch2 + total_ramp_samples_ch2,
                       ramp_stage2, cycle_idx_ch2 * sizeof(uint32_t));
            }

            total_ramp_samples_ch1 += cycle_idx_ch1;
            total_ramp_samples_ch2 += cycle_idx_ch2;
        }

        N += 1;
        if (verbose && update_interval > 0 && N % update_interval == 0) {
            printf("Progress: %lld%% complete\n", (N * 100LL) / simul->params->N_points);
            fflush(stdout);
        }
    }

    // Rewrite the headers with the real watermark and finalize
    if (simul->params->sample_ramp) {
        ramp_hdr1.samples_written = total_ramp_samples_ch1;
        ramp_hdr1.state = RAMP_STATE_FINALIZED;
        ramp_hdr2.samples_written = total_ramp_samples_ch2;
        ramp_hdr2.state = RAMP_STATE_FINALIZED;

        WriteRampHeader(f_ramp1, &ramp_hdr1);
        WriteRampHeader(f_ramp2, &ramp_hdr2);
        fclose(f_ramp1);
        fclose(f_ramp2);

        free(ramp_scratch1);
        free(ramp_scratch2);
        free(ramp_stage1);
        free(ramp_stage2);

        simul->ramp_samples_ch1 = total_ramp_samples_ch1;
        simul->ramp_samples_ch2 = total_ramp_samples_ch2;
    }

    if (verbose) printf("Simulation complete!\n");
}






// ===========================================================
// The API function intended for Python
// ===========================================================
#ifdef _WIN32
__declspec(dllexport)
#else
__attribute__((visibility("default"))) 
#endif
void run_simulation_api(SimulationParameters *params, 
                        const char *spectrum_filename,
                        const char *rate_filename, 
                        double custom_rate_gain,
                        double *out_ch1, double *out_ch2, 
                        uint32_t *ramp_ch1, uint32_t *ramp_ch2,
                        uint64_t random_seed) {
    
    SpectrumBuffer spectrum_buffer;
    InitSpectrumBuffer(&spectrum_buffer, 100);
    ReadSpectrumCSV(spectrum_filename, &spectrum_buffer);
    CalculateCDF(&spectrum_buffer);
    ComputeSpectrumIdxArray(&spectrum_buffer);
    
    params->spectrum = &spectrum_buffer;

    RateBuffer rate_buff;
    if (params->custom_rate == 1) {
        if (rate_filename == NULL || strlen(rate_filename) == 0) {
            fprintf(stderr, "Fatal: custom_rate is 1 but no rate_filename was provided.\n");
            exit(EXIT_FAILURE);
        }
        memset(&rate_buff, 0, sizeof(RateBuffer));
        InitRateBuffer(&rate_buff, 100);
        
        rate_buff.rate_gain = custom_rate_gain; 
        
        ReadRateCSV(rate_filename, &rate_buff);
        ComputeRateIdxArray(&rate_buff);
        params->rate_buff = &rate_buff;
    } else {
        params->rate_buff = NULL;
    }

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, random_seed, (intptr_t)&rng);
    params->rng = &rng;

    Simulation simul = { params, out_ch1, out_ch2, ramp_ch1, ramp_ch2, 0, 0 };
    // TODO: add verbosity
    RunSimulation(&simul, 0);
    if(params->output_in_charge == 0) {
        ConvertIntoChannel(&simul);
    }
    
    FreeSpectrumBuffer(&spectrum_buffer);
    if (params->custom_rate == 1) {
        free(rate_buff.idx_array);
        FreeRateBuffer(&rate_buff);
    }
}


// ===========================================================
// CLI Executable Code (Ignored when building Python library)
// ===========================================================
#ifndef BUILD_AS_LIB


void ParseConfig(const char* filename, SimulationParameters *params) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Could not open settings file '%s'\n", filename);
        exit(EXIT_FAILURE);
    }

    char key[64];
    char str_val[256];
    double double_val;
    char line[256];
    double sinusoid_period = 1;

    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '#' || line[0] == '\n') continue;

        // String-valued key: the output directory of the ramp files
        if (sscanf(line, "%63s %255s", key, str_val) == 2 &&
            strcmp(key, "ramp_dir") == 0) {
            free((void *)params->ramp_dir);
            params->ramp_dir = strdup(str_val);
            continue;
        }

        if (sscanf(line, "%63s %lf", key, &double_val) == 2) {
            if (strcmp(key, "t_int") == 0) params->t_int = double_val / 1000.0;
            else if (strcmp(key, "t_rst") == 0) params->t_rst = double_val / 1000.0;
            else if (strcmp(key, "rate") == 0) params->rate = double_val;
            else if (strcmp(key, "ramp_dt") == 0) params->ramp_dt = double_val / 1000000;
            else if (strcmp(key, "int_capacitor") == 0) params->int_capacitor = double_val;
            else if (strcmp(key, "ADC_min_range_V") == 0) params->ADC_min_range_V = double_val;
            else if (strcmp(key, "ADC_max_range_V") == 0) params->ADC_max_range_V = double_val;
            else if (strcmp(key, "electronic_noise") == 0) params->electronic_noise = double_val;
            else if (strcmp(key, "sinusoid_period") == 0) sinusoid_period = double_val / 1000;
            else if (strcmp(key, "sinusoid_amplitude") == 0) params->sinusoid_aplitude = double_val;
            else if (strcmp(key, "leakage_current") == 0) params->leakage_current = double_val;
            else if (strcmp(key, "relative_gain") == 0) params->relative_gain = double_val;
            else if (strcmp(key, "ADC_bit") == 0) params->ADC_bit = (int)double_val;
            else if (strcmp(key, "N_points") == 0) params->N_points = (int)double_val;
            else if (strcmp(key, "output_in_charge") == 0) params->output_in_charge = (int)double_val;
            else if (strcmp(key, "add_sinusoid") == 0) params->add_sinusoid = (int)double_val;
            else if (strcmp(key, "add_noise") == 0) params->add_noise = (int)double_val;
            else if (strcmp(key, "sample_ramp") == 0) params->sample_ramp = (int)double_val;
            else if (strcmp(key, "custom_rate") == 0) params->custom_rate = (int)double_val;
            else if (strcmp(key, "custom_rate_gain") == 0) params->rate_buff->rate_gain = double_val;

            params->sinusoid_omega = 2 * M_PI / sinusoid_period;

            if (params->custom_rate) {
                params->add_sinusoid = 0;
            }

        }
    }
    fclose(file);
    printf("Loaded settings from %s\n", filename);
}

void UserInput(SimulationParameters *params) {
    printf("Integration time [ms]: ");
    XSCANF(1, " %lf", &params->t_int);

    params->t_int /= 1000; // converted in seconds

    printf("Reset time [ms]: ");
    XSCANF(1, " %lf", &params->t_rst);


    params->t_rst /= 1000; // converted in seconds

    
    char answer;
    
    printf("Do you want to use a custom rate function? [y/n]: ");
    XSCANF(1, " %c", &answer); 
    if(answer == 'y' || answer == 'Y') {
        params->custom_rate = 1;
        printf("Relative rate gain: ");
        XSCANF(1, " %lf", &params->rate_buff->rate_gain);

    } else if(answer == 'n' || answer == 'N') {
        params->custom_rate = 0;
        printf("Rate [p/s]: ");
        XSCANF(1, " %lf", &params->rate);
    }

    printf("Number of iterations: ");
    XSCANF(1, " %d", &params->N_points);

    printf("Leakage current [pA]: ");
    XSCANF(1, " %lf", &params->leakage_current);
    
    printf("Relative gain: ");
    XSCANF(1, " %lf", &params->relative_gain);
    
    if (params->custom_rate == 0) {
        printf("Add beam fluctuations (sinusoid)? [y/n]: ");
        XSCANF(1, " %c", &answer); 
        
        if(answer == 'y' || answer == 'Y') {
            params->add_sinusoid = 1;
        }
        else if(answer == 'n' || answer == 'N') {
            params->add_sinusoid = 0;
        }
        else {
            fprintf(stderr, "Wrong answer\n");
            exit(EXIT_FAILURE);
        }
        

        if(params->add_sinusoid == 1) {
            double period;
            printf("Period [ms]: ");
            XSCANF(1, " %lf", &period);

            if(period <= 0) {
                fprintf(stderr, "The period can't be lower or equal to zero\n");
                exit(EXIT_FAILURE);
            }

            period /= 1000;  // convert into s
            params->sinusoid_omega = 2 * M_PI / period;
            
            printf("Amplitude [p/s]: ");
            XSCANF(1, " %lf", &params->sinusoid_aplitude);
        }
    }



    printf("Save data in Charge [pC] (1) or ADC channels [ch] (2)? [1/2]: ");
    XSCANF(1, " %c", &answer); 
    
    if(answer == '1' ) {
        params->output_in_charge = 1;
    }
    else if(answer == '2') {
        params->output_in_charge = 0;

        int adc_bit;
        printf("Number of ADC bit: ");
        XSCANF(1, "%d", &adc_bit);

        double adc_max;
        printf("ADC maximum voltage [V]: ");
        XSCANF(1, "%lf", &adc_max);

        double adc_min;
        printf("ADC minimum voltage [V]: ");
        XSCANF(1, "%lf", &adc_min);

        double int_capacitor;
        printf("Integration capacitance [pF]: ");
        XSCANF(1, "%lf", &int_capacitor);

        
        params->ADC_bit = adc_bit;
        params->ADC_max_range_V = adc_max;
        params->ADC_min_range_V = adc_min;
        params->int_capacitor = int_capacitor;

        printf("Do you want to add electronic noise? [y/n]: ");
        XSCANF(1, " %c", &answer); 
        if(answer == 'y' || answer == 'Y') {
            params->add_noise = 1;

            // The ADC noise is already considered when converting into channels

            double noise;
            printf("Electronic rms noise [mV]: ");
            XSCANF(1, " %lf", &noise);

            params->electronic_noise = noise / 1000;
            
            
        }
        else if(answer == 'n' || answer == 'N') {
            params->add_noise = 0;
        }
        else {
            fprintf(stderr, "Wrong answer\n");
            exit(EXIT_FAILURE);
        }


    }

    else {
        fprintf(stderr, "Wrong answer\n");
        exit(EXIT_FAILURE);
    }


    printf("Do you want to sample the ramp? [y/n]: ");
    XSCANF(1, " %c", &answer); 

    if(answer == 'y' || answer == 'Y') {
        params->sample_ramp = 1;

        double dt;
        printf("Time between ramp samples [us]: ");
        XSCANF(1, " %lf", &dt);

        params->ramp_dt = dt / 1000000;

        /* The ramp .bin files always store ADC channels (application AQR1
         * format): in charge mode the ADC parameters were never asked, but
         * the quantization needs them. */
        if (params->output_in_charge == 1) {
            printf("Ramp files store ADC channels: ADC parameters needed.\n");

            printf("Number of ADC bit: ");
            XSCANF(1, "%d", &params->ADC_bit);

            printf("ADC maximum voltage [V]: ");
            XSCANF(1, " %lf", &params->ADC_max_range_V);

            printf("ADC minimum voltage [V]: ");
            XSCANF(1, " %lf", &params->ADC_min_range_V);

            printf("Integration capacitance [pF]: ");
            XSCANF(1, " %lf", &params->int_capacitor);
        }
    }
    else if(answer == 'n' || answer == 'N') {
        params->sample_ramp = 0;
    }
    else {
        fprintf(stderr, "Wrong answer\n");
        exit(EXIT_FAILURE);
    }
    
}



int main(int argc, char *argv[]) {
    char *spectrum_filename = NULL;
    char *settings_filename = NULL;
    char *rate_filename = NULL;
    int use_settings = 0;
    int use_custom_rate = 0;
    int verbose = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--settings") == 0) {
            if (i + 1 < argc) {
                settings_filename = argv[++i];
                use_settings = 1;
            } else {
                fprintf(stderr, "Error: -s requires a filename argument.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--custom-rate") == 0) {
            if (i + 1 < argc) {
                rate_filename = argv[++i];
                use_custom_rate = 1;
            } else {
                fprintf(stderr, "Error: -r requires a filename argument.\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--verbose") == 0) {
            verbose = 1;
        } else {
            spectrum_filename = argv[i];
        }
    }

    if (!spectrum_filename) {
        printf("Usage: %s <spectrum_filename.csv> [-s settings.txt] [-r rate.csv]\n", argv[0]);
        return 1; 
    }

    SimulationParameters params; 
    memset(&params, 0, sizeof(SimulationParameters)); // Zero out the struct

    RateBuffer rate_buff;
    memset(&rate_buff, 0, sizeof(RateBuffer)); 
    InitRateBuffer(&rate_buff, 100);

    params.rate_buff = &rate_buff;



    if (use_settings) {
        ParseConfig(settings_filename, &params);
    } else {
        printf("No settings file provided. Falling back to manual input...\n");
        UserInput(&params);
    }
    
    SpectrumBuffer spectrum_buffer;
    InitSpectrumBuffer(&spectrum_buffer, 100);
    ReadSpectrumCSV(spectrum_filename, &spectrum_buffer);
    printf("Successfully read %zu data points from spectrum function.\n", spectrum_buffer.ptr + 1);
    CalculateCDF(&spectrum_buffer);

    ComputeSpectrumIdxArray(&spectrum_buffer);

    params.spectrum = &spectrum_buffer;

    pcg32_random_t rng;
    pcg32_srandom_r(&rng, time(NULL), (intptr_t)&rng);
    params.rng = &rng;


    if (params.custom_rate == 1) {

        if (use_custom_rate == 0) {
            fprintf(stderr, "Error: Provide a file for the custom rate.\n");
            exit(EXIT_FAILURE);
        }
        InitRateBuffer(&rate_buff, 100);
        ReadRateCSV(rate_filename, &rate_buff);
        printf("Successfully read %zu data points from custom rate function.\n\n", rate_buff.ptr +1);
        ComputeRateIdxArray(&rate_buff);
    }
    
    Simulation simul;
    memset(&simul, 0, sizeof(Simulation));
    simul.params = &params;

    simul.result_ch1 = xmalloc(simul.params->N_points * sizeof(*simul.result_ch1));
    simul.result_ch2 = xmalloc(simul.params->N_points * sizeof(*simul.result_ch2));

    // The per-cycle ramp scratch buffers are owned by RunSimulation; only
    // the ramp file paths/watermark live on afterwards (simul.ramp_out_*
    // stay NULL: the CLI consumes the files, not the arrays).

    RunSimulation(&simul, verbose);


    char header[256]; 

    snprintf(header, sizeof(header),
            "integration time [s]: %.15f, reset time [s]: %.15f, "
            "sample frequency [Hz]: %.15f, starting time : 2000-01-01 00:00:00.000",
            simul.params->t_int, 
            simul.params->t_rst,
            1.0 / (simul.params->t_int + simul.params->t_rst));


    if(params.output_in_charge == 0) {
        ConvertIntoChannel(&simul);
    }
    
    SaveDataToCSV(SAVE_FILE_NAME, header, &simul);

    free(simul.result_ch1);
    free(simul.result_ch2);
    free((void *)params.ramp_dir);
    if (params.custom_rate == 1) {
        free(rate_buff.idx_array);
        FreeRateBuffer(&rate_buff);
    }
    FreeSpectrumBuffer(&spectrum_buffer);

    return 0;
}



#endif // BUILD_AS_LIB
