# TelecomUtils
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://disberd.github.io/TelecomUtils.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://disberd.github.io/TelecomUtils.jl/dev)

This package provides some basic utilities for helping with simulations of Telecommunication Satellites. The main contribution comes in the form of the `ReferenceView` structure defined in the [`refview_struct.jl`](notebooks/refview_struct.jl) notebook.

## ReferenceView 
The `ReferenceView` structure, together with the associated methods, allow to compute visibility, distance and pointing angles between points in space or on ground. For each of the functions it is possible to specify a reference _face_ to use for computing visibilities. Each object is considered to be a cube, so each of its 6 faces can be selected as a reference.

See the [`refview_struct.jl`](notebooks/refview_struct.jl) or the docstring for `ReferenceView` and the docstrings of other cuntions referenced for some example usage.

## Spectral efficiency computation
Some basic functionality to map between Signal to Noise Ratio (SNR) and Spectral efficiency for both DVB-S2x and 5G NR air interfaces is provided in the [`snr2speff.jl`](notebooks/snr2speff.jl) notebook.

Check the notebook directly or look at the docstring of [`speff2snr`](@ref) and the docstring of other functions referenced therein for some example usage.