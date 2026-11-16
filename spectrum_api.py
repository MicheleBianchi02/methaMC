import ctypes
import os

import numpy as np

# np.set_printoptions(formatter={"float": lambda x: "{0:0.5f}".format(x)})


class SimulationParameters(ctypes.Structure):
    _fields_ = [
        ("spectrum", ctypes.c_void_p),
        ("rate_buff", ctypes.c_void_p),
        ("rng", ctypes.c_void_p),
        ("t_int", ctypes.c_double),
        ("t_rst", ctypes.c_double),
        ("rate", ctypes.c_double),
        ("ramp_dt", ctypes.c_double),
        ("int_capacitor", ctypes.c_double),
        ("ADC_min_range_V", ctypes.c_double),
        ("ADC_max_range_V", ctypes.c_double),
        ("electronic_noise", ctypes.c_double),
        ("sinusoid_omega", ctypes.c_double),
        ("sinusoid_aplitude", ctypes.c_double),
        ("leakage_current", ctypes.c_double),
        ("relative_gain", ctypes.c_double),
        ("ADC_bit", ctypes.c_int),
        ("N_points", ctypes.c_int),
        ("output_in_charge", ctypes.c_int),
        ("add_sinusoid", ctypes.c_int),
        ("add_noise", ctypes.c_int),
        ("sample_ramp", ctypes.c_int),
        ("custom_rate", ctypes.c_int),
        ("ramp_dir", ctypes.c_char_p),
    ]

_SO_DIR = os.path.dirname(os.path.abspath(__file__))
_SO_PATH = os.path.join(_SO_DIR, "libspectrum.so")
lib = ctypes.CDLL(_SO_PATH)
lib.run_simulation_api.argtypes = [
    ctypes.POINTER(SimulationParameters),
    ctypes.c_char_p,
    ctypes.c_char_p,
    ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_uint32),
    ctypes.POINTER(ctypes.c_uint32),
    ctypes.c_uint64,
]
lib.run_simulation_api.restype = None


def run_monte_carlo(
    config_dict,
    spectrum_file,
    seed,
    rate_file=None,
    custom_rate_gain=1.0,
    ramp_dir=None,
):
    """Run the Monte Carlo through the C library.

    With ``sample_ramp = 1`` the ramp of every integration window is stored
    twice, in the same quantized uint32 ADC channels:

    - on disk, as AQR1 files (application format, ``services/ramp_storage.py``)
      written to ``ramp_dir`` (default: the current directory) as
      ``simulation_ramp_ch1.bin`` / ``simulation_ramp_ch2.bin``;
    - in memory, returned as the two trailing numpy uint32 arrays (the full
      run, window after window — identical to the file contents).

    Returns ``(ch1, ch2)`` integration results (float64; charge or channels
    depending on ``output_in_charge``) plus ``(ramp_ch1, ramp_ch2)`` when the
    ramp is sampled.
    """
    clean_config = config_dict.copy()
    clean_config.pop("custom_rate_gain", None)
    if ramp_dir is None:
        ramp_dir = clean_config.pop("ramp_dir", None)

    params = SimulationParameters(**clean_config)

    ch1_out = np.zeros(params.N_points, dtype=np.float64)
    ch2_out = np.zeros(params.N_points, dtype=np.float64)

    # Handle conditional ramp array allocation
    ramp_ch1_ptr = ctypes.cast(None, ctypes.POINTER(ctypes.c_uint32))
    ramp_ch2_ptr = ctypes.cast(None, ctypes.POINTER(ctypes.c_uint32))
    ramp_ch1_out = None
    ramp_ch2_out = None

    if getattr(params, "sample_ramp", 0) == 1:
        if getattr(params, "ramp_dt", 0.0) <= 0:
            raise ValueError("sample_ramp=1 requires ramp_dt > 0")
        if params.ramp_dt > params.t_int:
            raise ValueError(
                "ramp_dt > t_int: every window needs at least 2 samples "
                "(application requirement)"
            )
        if not 1 <= params.ADC_bit <= 32:
            raise ValueError("sample_ramp=1 requires 1 <= ADC_bit <= 32")

        # Same formula as RampSamplesPerWindow() in spectrum.c: every window
        # stores exactly this many samples, so the run total is exact.
        samples_per_window = int(params.t_int / params.ramp_dt + 0.5) + 1
        total_samples = params.N_points * samples_per_window

        ramp_ch1_out = np.zeros(total_samples, dtype=np.uint32)
        ramp_ch2_out = np.zeros(total_samples, dtype=np.uint32)

        ramp_ch1_ptr = ramp_ch1_out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))
        ramp_ch2_ptr = ramp_ch2_out.ctypes.data_as(ctypes.POINTER(ctypes.c_uint32))

    # Output directory of the AQR1 ramp files (None = current directory).
    # Created here (recursively, like the application does) because the C
    # side only mkdirs a single level as a CLI convenience.
    if ramp_dir:
        os.makedirs(ramp_dir, exist_ok=True)
    params.ramp_dir = os.fsencode(ramp_dir) if ramp_dir else None

    # Encode filenames. Handles None safely.
    c_spectrum_file = spectrum_file.encode("utf-8")
    c_rate_file = rate_file.encode("utf-8") if rate_file is not None else None

    lib.run_simulation_api(
        ctypes.byref(params),
        c_spectrum_file,
        c_rate_file,
        ctypes.c_double(custom_rate_gain),
        ch1_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ch2_out.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        ramp_ch1_ptr,
        ramp_ch2_ptr,
        ctypes.c_uint64(seed),
    )

    if getattr(params, "sample_ramp", 0) == 1:
        return ch1_out, ch2_out, ramp_ch1_out, ramp_ch2_out

    return ch1_out, ch2_out


# TODO:
# - plottare con log y gli errori (nel report.tex) : calcolare pendenza, dovrebbe essere
#    -1/2
# - simulare con tante N diverse e misurare yd e incertezza e comparare con valore
#    teorico
# - calcolare incertezza dell'incertezza analiticamente
# - per var-cov runnare tante volte con stesso N e misure yd e incertezza. graficare
#    spettro frequenza di yd e incertezza e misurare se std di yd (con diverse simulazioni)
#    corrisponde all'incertezza misurata. da spettro frequenza incertezza mostrare
#    che la stima dell'incertezza è accurata (misurare incertezza dell'incertezza = std
#    a diversi N)
# - fare tutto con spettri gaussiani (tint > 100ms per rateo = 8000)
# - salvare i dati in un file csv e fare i plot con latex invece che con python
# - simulare con stessi parametri ma cambiando il block size per vedere qual'è quello ottimale
