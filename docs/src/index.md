# Introduction

## Overview

AutocorrelationShell is package for using autocorrelation wavelets in [julia](https://github.com/JuliaLang/julia).

Signal representations using autocorrelation wavelets are redundant and non-orthogonal. Some desirable properties of autocorrelation wavelet transforms are symmetry without losing vanishing moments, edge detection and characterization capabilities, and shift invariance. Autocorrelation wavelets can be used as a tool for data analysis such as time series analysis and image analysis.

The Autocorrelation Shell is heavily inspired by the paper [Wavelets, their autocorrelation functions, and multiresolution representation of signals](https://www.math.ucdavis.edu/~saito/publications/saito_acs_spie.pdf).

## Installation

You can install the package at the Pkg REPL-mode with:
````julia
(v1.4) >] add AutocorrelationShell
````

## Highlights

- 1D autocorrelation wavelet transform.
- 2D autocorrelation wavelet transform.
- 1D autocorrelation wavelet packet transform.
- A variety of utility functions for image de-noising such as functions for Signal to Noise Ratio(SNR) and Structural Similarity Index(SSIM).
- The package is completely documented.

## Authors

This package was authored by [Rishi Subramanian](https://www.linkedin.com/in/rishi-subramanian-50550b104/), [Christina Chang](https://www.linkedin.com/in/christina-l-chang/), and [Shozen Dan](https://www.linkedin.com/in/shozendan/) under the supervision of Professor [Naoki Saito](https://www.math.ucdavis.edu/~saito/) at University of California, Davis.
