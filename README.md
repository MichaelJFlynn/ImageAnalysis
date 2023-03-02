# ImageAnalysis

This repository contains code used in analyzing fluorescence microscopy images of neurons with HCR dots,
including segmentations of cells and dot calling.

This code is a personal library and has not been optimized for readability in any way.
However, it is important for scientific reproduciblity to make code public, so I have uploaded it here. 

The code currently requires Windows-specific concurrency libraries `ppl.h` and `ppltasks.h` and will not run on
Mac or Linux.

The code requires the installation of `OpenCV4` and `libTiff` via `vcpkg` in order to work,
and `cmake` to install(see `CMakeLists.txt`).

The library is contained within the `ImageAnalysis.hpp` file. For example usage, check out `example.cpp`.

