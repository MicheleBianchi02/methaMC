# Microdosimetric Multi-Event Monte Carlo Simulator

methaMC (Multi-Event Tool for HAdrontherapy - Monte Carlo) is a custom Monte
Carlo code written in C (with Python bindings and analysis tools) designed to
simulate charge collection on a gated-integrator readout circuit coupled to a
microdosimetric detector (such as a mini-TEPC). 

Developed within the framework of the **INFN AToMiQA** project, this simulator
models stochastic single-event charge depositions under both stationary and
time-varying particle beams. It is used to benchmark multi-event
variance-covariance estimators of the dose-mean lineal energy $\bar{y}_D$ in
the time and frequency domains.

For a detailed theoretical and algorithmic review of the simulator
architecture, flowcharts, and post-processing stages, please refer to
**[`Architecture.pdf`](.docs/Architecture.pdf)** included in this repository.

---

## Repository Structure

- **`spectrum.c`** — Core Monte Carlo simulation engine in C. Simulates two independent detector channels (`ch1` and `ch2`), continuous leakage current, optional high-frequency ramp sampling (saved in binary format), Gaussian electronic noise, and ADC quantization. Can be compiled as a standalone CLI executable or as a shared library.

- **`Makefile`** — Build configuration to compile both the standalone CLI executable (`spectrum`) and the shared C library (`libspectrum.so`).

- **`config.txt`** — Example configuration file specifying timing parameters ($`t_{\text{int}}`$, $`t_{\text{rst}}`$, `ramp_dt`), circuit characteristics ($`C_{\text{int}}`$, leakage current, relative gain), ADC settings, and beam modulation flags.

- **`spectrum_api.py`** — Python `ctypes` wrapper (`run_monte_carlo`) around `libspectrum.so` for running simulations and retrieving integrated charge/channel arrays and ramp data directly as NumPy arrays.

- **`generate_rate.py`** — PyQt6 GUI and Python API (`generate_rate_file`) to generate custom time-varying relative beam rate profiles ($`\lambda(t)`$) in CSV format, supporting sinusoidal rates, band-limited Gaussian noise, pulse trains, and sudden beam drops.

- **`api/calc_stat.py`** — Post-processing and statistical analysis module implementing multi-event estimators and circular block bootstrap uncertainty estimation.

- **`docs/Architecture.pdf`** — In-depth documentation of the Monte Carlo simulator architecture and physics loop.

- **`docs/output_file_formats.md`** — Documentation regarding the output files.

---

## Prerequisites & Dependencies

### 1. External C Dependency: PCG Random Number Generator
This simulator relies on the minimal C implementation of the **PCG32** random number generator (`pcg_basic.c` and `pcg_basic.h`), which is not bundled directly in this repository.

Before compiling, download `pcg_basic.c` and `pcg_basic.h` from the official [PCG C Basic repository](https://github.com/imneme/pcg-c-basic) (or [pcg-random.org](https://www.pcg-random.org/download.html)) and place them in the root directory of this project:

```bash
curl -O https://raw.githubusercontent.com/imneme/pcg-c-basic/master/pcg_basic.c
curl -O https://raw.githubusercontent.com/imneme/pcg-c-basic/master/pcg_basic.h
```

### 2. Python Environment
To use the Python API, rate generator GUI, and statistical analysis routines, install the required packages:

```bash
pip install numpy scipy arch PyQt6 pyqtgraph
```

---

## Compilation

Once `pcg_basic.c` and `pcg_basic.h` are in the project directory, compile both the CLI executable and the Python shared library using `make`:

```bash
make
```

This produces:
* `spectrum` (standalone command-line executable)
* `libspectrum.so` (shared library loaded by `spectrum_api.py`)

To clean build artifacts, run `make clean`.

While the simulator support both Linux and windows machines, `make` can only be 
used on Linux. 



---

## Usage

### 1. Command-Line Interface (CLI)
Run the simulator by passing a single-event charge/lineal-energy spectrum CSV file (`x, y` columns with header), along with optional configuration and custom rate files:

```bash
./spectrum <spectrum_file.csv> [-s config.txt] [-r custom_rate.csv] [-v]
```

* `-s`, `--settings`: Path to the configuration file (e.g., `config.txt`). If omitted, an interactive prompt will ask for the parameters.
* `-r`, `--custom-rate`: Path to a custom time-varying rate CSV file (required if `custom_rate` is set to `1`).
* `-v`, `--verbose`: Print simulation setup and progress updates.

The output is saved to `simulation_results.csv` (and `simulation_ramp_ch1.bin` / `simulation_ramp_ch2.bin` if `sample_ramp` is enabled).
The format is specified [here](./docs/output_file_formats.md).


### 2. Custom Beam Rate Generation
Launch the interactive GUI to design, inspect (Time Domain, FFT, Sinc Integrator Filter, and Autocorrelation), and export a `custom_rate.csv` profile:

```bash
python generate_rate.py
```

Alternatively, call `generate_rate_file()` directly from your Python scripts.

### 3. Python API & Statistical Analysis
You can run the simulation and evaluate the multi-event estimators directly in Python:

```python
import numpy as np
from spectrum_api import run_monte_carlo
from calc_stat import calculate_results

config = {
    "t_int": 0.488e-3,          # [s]
    "t_rst": 0.427e-3,          # [s]
    "rate": 9.6e6,              # [p/s]
    "ramp_dt": 1e-6,            # [s]
    "int_capacitor": 22.0,      # [pF]
    "ADC_min_range_V": 0.0,     # [V]
    "ADC_max_range_V": 5.0,     # [V]
    "electronic_noise": 0.0003, # [V]
    "sinusoid_omega": 2 * np.pi / 1.0,
    "sinusoid_aplitude": 3000.0,
    "leakage_current": 0.0,     # [pA]
    "relative_gain": 0.0057,
    "ADC_bit": 12,
    "N_points": 20000,
    "output_in_charge": 1,
    "add_sinusoid": 0,
    "add_noise": 0,
    "sample_ramp": 0,
    "custom_rate": 0,
}

# Run Monte Carlo simulation
ch1, ch2 = run_monte_carlo(config, spectrum_file="spectrum.csv", seed=42)

# Compute Variance-Covariance and FBA estimators
data = np.column_stack((ch1, ch2))
fs = 1.0 / (config["t_int"] + config["t_rst"])
results = calculate_results(data, cov_both_ch=True, calc_fft=True, sample_freq_hz=fs)
```

---

## Citation & Credits

* **Author & Thesis Reference**:  
  Michele Bianchi, *"First implementation of the Variance-Covariance Method applied to Microdosimetry of ion beams. The INFN AToMiQA experiment"*, M.Sc. Thesis in Nuclear Engineering, Politecnico di Milano, Academic Year 2025–2026.
* **PCG Random Number Generator**:  
  Melissa E. O'Neill, *"PCG: A Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number Generation"*, Harvey Mudd College Computer Science Department Technical Report HMC-CS-2014-0905, 2014. Code available at [github.com/imneme/pcg-c-basic](https://github.com/imneme/pcg-c-basic) under the Apache License 2.0.
