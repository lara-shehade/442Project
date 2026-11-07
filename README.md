# EECE 442 Course Project, Part I: Source Encoder and Decoder

MATLAB code for Part I of the EECE 442 (Communication Systems) course project at the American University of Beirut, Fall 2026.

Part I builds and tests the source side of a digital communication chain:

```
m(t) -> sampling -> quantization -> source coding -> bits -> decoding -> reconstruction -> m^(t)
```

The message signal is the opening of Beethoven's *Für Elise* (E5, D#5, E5, D#5, E5). It is modelled as a sum of 22 damped sinusoids, 4 to 5 harmonics per note, measured from a real piano recording. Because the model is a formula, `m(t)` can be evaluated at any time instant and used as an exact reference.

**Team:** Lara Chehade, Joey Jamous, Jonathan Achkouti

## Requirements

- MATLAB R2022a or newer (tested on R2025b). No toolbox is needed for the main results.
- Audio Toolbox (R2024a or newer), optional. It is only used for the ViSQOL perceived-quality scores in `part2.m`. Without it, that block is skipped and everything else runs.

## Quick start

1. Put all the files in one folder and make it MATLAB's current folder.
2. In the Command Window, run:

   ```matlab
   run_all_figs
   ```

This runs the whole project in order, in about 3 to 5 minutes. It produces:

| Output | Content |
|---|---|
| `figs_1_1/`, `figs_1_2/`, `figs_part2/`, `figs_part3/` | every figure as a PNG, one folder per section |
| `all_figs.zip` | the four figure folders, zipped |
| `run_log.txt` | everything the scripts print (tables, parameters, MSE, bit rates) |
| `quantized.mat` | the quantized sequences, written by `part2.m` and read by `part3.m` |
| `figs/`, `figs.zip` | also written by `fourierApp1_2.m` itself (the same figures as `figs_1_2/`) |

`part2.m` plays the original message, the sampled version and every quantized version, about 15 s of sound in total.

Each script starts with `clc`, so after `run_all_figs` the Command Window only shows the output of the last script. The full output of every script is in `run_log.txt`.

## Project structure

The project follows the three modules of the Part I specification.

### Section 1: Signal representation through sampling

| File | Type | What it does |
|---|---|---|
| `compare_1.m` | script | 1.1: samples `cos(2*pi*250*t)` below Nyquist (250 Hz) and above Nyquist (1000 Hz), reconstructs both and prints the MSE |
| `fourierApp1_2.m` | script | 1.2: builds the message, finds its bandwidth and sampling rate, runs the Fourier series experiments (a) to (e), checks the reconstruction of the message, and applies the DTFT on frequency bands |
| `piano_signal.m` | function | the message signal `m = piano_signal(t)` |
| `sample.m` | function | `[t_sample, x_sample] = sample(t, xt, fs)` |
| `reconstruct.m` | function | `xrcon = reconstruct(t, x_sample, fs)`, ideal sinc interpolation |
| `ffs.m` | function | `[xhat, ck] = ffs(xt, t, n, T)`, finite Fourier series computed from the integral (no built-in transform) |

Main values: bandwidth `B = 2516 Hz` (99.99% of the energy), `fN = 2B = 5032 Hz`, `fs = 2fN = 10064 Hz`, giving `N = 12480` samples. The Fourier series uses a 50 ms chunk starting at 0.60 s, with `T = 50 ms`.

### Section 2: Quantization

| File | Type | What it does |
|---|---|---|
| `part2.m` | script | two-level quantizer (optimal, M = 2), uniform quantizers (M = 4, 16, 32), Lloyd-Max quantizers from the empirical distribution (M = 4, 16, 32), the end-to-end chain with MSE of the reconstructed message, ViSQOL scores, and playback |
| `quan.m` | function | `xq = quan(x, thr, lvl)` |

It prints the thresholds and levels of every quantizer and the MSE table, and it saves `quantized.mat` for Section 3.

### Section 3: Source coding and decoding

| File | Type | What it does |
|---|---|---|
| `part3.m` | script | fixed-length coding, entropy, symbol-wise Huffman, block Huffman (k = 2, 3, 4), total bits, throughput and latency, the distortion-versus-rate comparison, and LZ77 (bonus) with the maximum match length tuned |
| `huff.m` | function | `codes = huff(p)`, Huffman code for the probabilities `p` |
| `huff_dec.m` | function | `s = huff_dec(bits, codes)`, decoder for a prefix-free code |
| `lz77.m` | function | `tok = lz77(s, W, L)`, LZ77 encoder; each token is (offset, length, next symbol) |
| `lz77dec.m` | function | `s = lz77dec(tok)`, LZ77 decoder |

Every code (Huffman, block Huffman and LZ77) is decoded and compared with the original symbol stream, and the script prints whether the reconstruction is lossless.

### Helpers

| File | What it does |
|---|---|
| `run_all_figs.m` | runs all four scripts in order, saves the figures and the printed output, and zips the figures |
| `save_figs.m` | `save_figs(folder)` saves every open figure as a PNG with a white background |

## Running one part at a time

Each script can also be run on its own:

```matlab
compare_1        % Section 1.1
fourierApp1_2    % Section 1.2
part2            % Section 2 (writes quantized.mat)
part3            % Section 3 (needs quantized.mat, so run part2 first)
```

To save the figures of a script, call `save_figs` right after it, for example `part2; save_figs('figs_part2')`.

## Notes

- Continuous time is represented by a 44.1 kHz grid (a 48 kHz grid for ViSQOL).
- Reconstruction errors are measured without the first and last 50 ms of the message, where the sinc pulses of samples outside the message are missing.
- Nothing is random, so every run gives the same numbers.
