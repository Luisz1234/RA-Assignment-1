# Galton Box Simulation

This repository contains the C++ source code for Assignment 1: Galton Box Simulation. The program simulates the Galton Box experiment, analyzes the statistical distributions, and generates visual output.

## Program Overview

The `galtonbox.cpp` program performs the following:
1.  **Simulation:** Simulates `N` balls falling through an `n`-level Galton board.
2.  **Distribution Calculation:** Computes empirical frequencies, the theoretical binomial PMF $P(X=k) = \binom n k (0.5)^n$, and its normal approximation $\mathcal{N}(n/2, n/4)$.
3.  **Error Analysis:** Calculates Mean Quadratic Error (MQE) for:
    *   Binomial PMF vs. Normal PDF
    *   Experimental Frequencies vs. Binomial PMF
    *   Experimental Frequencies vs. Normal PDF
4.  **Automated Plotting:** Generates `.png` figures for visual comparison of the distributions using Gnuplot.

## Prerequisites

To compile and run the program, the following software must be installed and accessible in your system's PATH:

1.  **C++ Compiler:** A C++11 (or newer) compatible compiler (e.g., `g++` or Clang).
2.  **Gnuplot:** The plotting utility.

## Instructions

1.  **Compilation:**
    Navigate to the directory containing `galton_simulator.cpp` and compile using a C++ compiler:
    
    g++ -std=c++11 -O2 galton_simulator.cpp -o galton_simulator -lm
    

2.  **Execution:**
    Run the compiled executable. The program will prompt for the number of levels (`n`) and the number of balls (`N`).
    
    ./galton_simulator
    
    Upon execution, the program will ask for n and N values anf then print tabular results (frequencies and probabilities), MQE values, and generate data files (`.txt`) and plot images (`.png`) in the current directory.

## Generated Output

Running the program produces:
*   Standard output to the console, including tabular data and MQE values.
*   `.txt` files: Data points used for plotting.
*   `.png` files: Visual comparisons of the simulated, binomial, and normal distributions.
