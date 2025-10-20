#include <iostream>  // For input/output operations (cin, cout, cerr)
#include <vector>    // For dynamic arrays (std::vector)
#include <cmath>     // For mathematical functions (sqrt, exp, pow)
#include <iomanip>   // For output formatting (setw, setprecision, fixed)
#include <random>    // For random number generation (mt19937, uniform_real_distribution)
#include <numeric>   // Not directly used but good for general includes
#include <algorithm> // For std::max
#include <fstream>   // For file I/O operations (ofstream)
#include <string>    // For std::string
#include <cstdio>    // For std::remove

int n;
long long N;


/**
 * @brief Calculates "n choose r" (nCr) combination.
 * Uses an iterative method to avoid large factorial calculations and potential overflow.
 * Handles nCr(n, r) == nCr(n, n-r) optimization.
 * @param n Total number of items.
 * @param r Number of items to choose.
 * @return The combination value as a long double for precision.
 */
long double nCr(int n, int r) {
    if (r < 0 || r > n) {
        return 0; // Invalid input for combinations
    }
    if (r == 0 || r == n) {
        return 1; // Base cases: choose 0 or all items
    }
    if (r > n / 2) {
        r = n - r; // Optimize calculation: nCr(n, r) == nCr(n, n-r)
    }

    long double res = 1;
    for (int i = 1; i <= r; ++i) {
        // Calculation: C(n,r) = (n * (n-1) * ... * (n-r+1)) / (r * (r-1) * ... * 1)
        res = res * (n - i + 1) / i;
    }
    return res;
}

/**
 * @brief Calculates the Probability Mass Function (PMF) for a Binomial distribution.
 * Specifically for p = 0.5, as required by the Galton board.
 * P(X=k) = C(n, k) * p^k * (1-p)^(n-k)
 * With p=0.5, it simplifies to C(n, k) * (0.5)^n
 * @param k The number of "successes" (right turns in this case).
 * @param n The total number of trials (levels in this case).
 * @return The probability P(X=k) as a double.
 */
double binomialPMF(int k, int n) {
    if (k < 0 || k > n) {
        return 0.0; // Invalid k for binomial distribution
    }
    long double combinations = nCr(n, k);
    // (0.5)^k * (0.5)^(n-k) simplifies to (0.5)^n
    double prob_power = std::pow(0.5, n);
    return static_cast<double>(combinations * prob_power);
}

/**
 * @brief Calculates the Probability Density Function (PDF) for a Normal distribution.
 * f(x) = (1 / (stdDev * sqrt(2 * PI))) * exp(-0.5 * ((x - mean) / stdDev)^2)
 * @param x The value at which to evaluate the PDF.
 * @param mean The mean (mu) of the normal distribution.
 * @param stdDev The standard deviation (sigma) of the normal distribution.
 * @return The density f(x) as a double.
 */
double normalPDF(double x, double mean, double stdDev) {
    if (stdDev <= 0) {
        return 0.0; // Standard deviation must be positive
    }
    double exponent = -0.5 * std::pow((x - mean) / stdDev, 2);
    return (1.0 / (stdDev * std::sqrt(2.0 * M_PI))) * std::exp(exponent);
}

/**
 * @brief Simulates a single ball dropping through the Galton board.
 * At each level, the ball has a 50% chance to go left or right.
 * @param n The number of levels in the Galton board.
 * @param rng A random number generator engine (e.g., std::mt19937).
 * @return The final slot index (number of right turns) where the ball lands.
 */
int simulateBall(int n, std::mt19937& rng) {
    // Uniform distribution to simulate a 50/50 chance (0.0 to 1.0)
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    int right_turns = 0; // Counter for right turns (final slot index)
    for (int i = 0; i < n; ++i) {
        if (dist(rng) >= 0.5) { // If random number is 0.5 or greater, assume it goes right
            right_turns++;
        }
    }
    return right_turns;
}

/**
 * @brief Saves the experimental frequencies, binomial probabilities, and normal densities
 *        to a text file for plotting with Gnuplot.
 * @param n The number of levels.
 * @param experimental_frequencies Vector of experimental frequencies.
 * @param binomial_probabilities Vector of binomial probabilities.
 * @param normal_densities Vector of normal densities.
 * @param filename The name of the file to save data to.
 * @return true if data was saved successfully, false otherwise.
 */
bool saveDataForPlotting(int n,
                         const std::vector<double>& experimental_frequencies,
                         const std::vector<double>& binomial_probabilities,
                         const std::vector<double>& normal_densities,
                         const std::string& filename) {
    std::ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open data file " << filename << " for writing.\n";
        return false;
    }

    // Write header
    outputFile << "Slot ExperimentalFreq BinomialPMF NormalPDF\n";
    outputFile << std::fixed << std::setprecision(10); // High precision for data points

    for (int i = 0; i <= n; ++i) {
        outputFile << i << " "
                   << experimental_frequencies[i] << " "
                   << binomial_probabilities[i] << " "
                   << normal_densities[i] << "\n";
    }
    outputFile.close();
    std::cout << "\nData saved to " << filename << " for plotting.\n";
    return true;
}

/**
 * @brief Generates and executes a Gnuplot script to create a PNG image.
 * This version generates a plot similar to Figure 2 (Binomial PMF vs Normal PDF).
 * @param dataFilename The name of the data file.
 * @param outputPngFilename The desired name for the output PNG image.
 * @param plotTitle The title of the plot.
 * @return true if Gnuplot command was executed, false otherwise.
 */
bool generateBinomialVsNormalPlot(const std::string& dataFilename,
                                     const std::string& outputPngFilename,
                                     const std::string& plotTitle) {
    std::string scriptFilename = "plot_binomial_normal_n" + std::to_string(n) + "_N" + std::to_string(N) + ".gp";
    std::ofstream scriptFile(scriptFilename);
    if (!scriptFile.is_open()) {
        std::cerr << "Error: Could not open Gnuplot script file " << scriptFilename << " for writing.\n";
        return false;
    }

    scriptFile << "set terminal pngcairo enhanced font 'Arial,10' size 800,600\n";
    scriptFile << "set output '" << outputPngFilename << "'\n";
    scriptFile << "set title '" << plotTitle << "'\n";
    scriptFile << "set xlabel 'k (Slot)'\n";
    scriptFile << "set ylabel 'P(X=k) / Probability Density'\n";
    scriptFile << "set yrange [0:*]\n"; // Ensure y-axis starts at 0
    scriptFile << "set grid\n";
    scriptFile << "set style data histograms\n";
    scriptFile << "set style fill solid 0.5 border lc rgb 'black'\n"; // Filled bars with black border
    scriptFile << "set boxwidth 0.9\n";
    scriptFile << "set xtics nomirror\n"; // Ensure x-axis ticks are consistent

    // Plotting Binomial PMF as boxes and Normal PDF as a line
    // Binomial PMF is in column 3, Normal PDF in column 4
    scriptFile << "plot '" << dataFilename << "' using 1:3 with boxes title 'Binomial PMF' lc rgb '#ADD8E6', \\\n"; // Light blue boxes
    scriptFile << "     '' using 1:4 with lines title 'Normal PDF' lw 2 lc rgb 'black'\n"; // Solid black line
    
    scriptFile.close();
    std::cout << "Gnuplot script saved to " << scriptFilename << ".\n";

    std::string gnuplotCommand = "gnuplot " + scriptFilename;
    int result = system(gnuplotCommand.c_str());

    if (result == 0) {
        std::cout << "Plot generated successfully: " << outputPngFilename << "\n";
    } else {
        std::cerr << "Error: Gnuplot command failed. Make sure Gnuplot is installed and in your PATH.\n";
    }
    std::remove(scriptFilename.c_str()); // Clean up script file
    return (result == 0);
}

/**
 * @brief Generates and executes a Gnuplot script to create a PNG image.
 * This version compares Experimental Frequencies (bars) with Binomial PMF (line).
 * @param dataFilename The name of the data file.
 * @param outputPngFilename The desired name for the output PNG image.
 * @param plotTitle The title of the plot.
 * @return true if Gnuplot command was executed, false otherwise.
 */
bool generateExperimentalVsBinomialPlot(const std::string& dataFilename,
                                        const std::string& outputPngFilename,
                                        const std::string& plotTitle) {
    std::string scriptFilename = "plot_experimental_binomial_n" + std::to_string(n) + "_N" + std::to_string(N) + ".gp";
    std::ofstream scriptFile(scriptFilename);
    if (!scriptFile.is_open()) {
        std::cerr << "Error: Could not open Gnuplot script file " << scriptFilename << " for writing.\n";
        return false;
    }

    scriptFile << "set terminal pngcairo enhanced font 'Arial,10' size 800,600\n";
    scriptFile << "set output '" << outputPngFilename << "'\n";
    scriptFile << "set title '" << plotTitle << "'\n";
    scriptFile << "set xlabel 'k (Slot)'\n";
    scriptFile << "set ylabel 'P(X=k) / Frequency'\n";
    scriptFile << "set yrange [0:*]\n"; // Ensure y-axis starts at 0
    scriptFile << "set grid\n";
    scriptFile << "set style data histograms\n";
    scriptFile << "set style fill solid 0.7 border lc rgb 'black'\n"; // Filled bars with black border
    scriptFile << "set boxwidth 0.9\n";
    scriptFile << "set xtics nomirror\n"; // Ensure x-axis ticks are consistent

    // Plotting Experimental Frequencies as boxes and Binomial PMF as points/line for comparison
    // Experimental Freq is in column 2, Binomial PMF in column 3
    scriptFile << "plot '" << dataFilename << "' using 1:2 with boxes title 'Experimental Freq' lc rgb '#B0E0E6', \\\n"; // Powder blue boxes
    scriptFile << "     '' using 1:3 with linespoints title 'Binomial PMF' lw 1 pt 7 lc rgb 'red'\n"; // Red line with points
    
    scriptFile.close();
    std::cout << "Gnuplot script saved to " << scriptFilename << ".\n";

    std::string gnuplotCommand = "gnuplot " + scriptFilename;
    int result = system(gnuplotCommand.c_str());

    if (result == 0) {
        std::cout << "Plot generated successfully: " << outputPngFilename << "\n";
    } else {
        std::cerr << "Error: Gnuplot command failed. Make sure Gnuplot is installed and in your PATH.\n";
    }
    std::remove(scriptFilename.c_str()); // Clean up script file
    return (result == 0);
}


int main() {
    std::cout << "Welcome to the Galton Board Simulator!\n";

    // Get input for number of levels
    std::cout << "Enter the number of levels (n) [e.g., 10, 20]: ";
    std::cin >> n;
    if (n <= 0) {
        std::cerr << "Error: Number of levels must be positive.\n";
        return 1;
    }

    // Get input for number of balls
    std::cout << "Enter the number of balls (N) [e.g., 10000, 100000]: ";
    std::cin >> N;
    if (N <= 0) {
        std::cerr << "Error: Number of balls must be positive.\n";
        return 1;
    }

    // Initialize a high-quality random number generator
    std::random_device rd;  // Seeds the random number generator
    std::mt19937 rng(rd()); // Mersenne Twister engine

    // --- Simulation ---
    std::cout << "\nSimulating " << N << " balls through a " << n << "-level Galton board...\n";
    // Vector to store counts for each landing slot (0 to n)
    std::vector<int> experimental_counts(n + 1, 0);

    for (long long i = 0; i < N; ++i) {
        int slot = simulateBall(n, rng);
        // Ensure the slot is within the valid range [0, n]
        if (slot >= 0 && slot <= n) {
            experimental_counts[slot]++;
        }
    }

    // --- Calculate Experimental Frequencies ---
    std::vector<double> experimental_frequencies(n + 1);
    for (int i = 0; i <= n; ++i) {
        experimental_frequencies[i] = static_cast<double>(experimental_counts[i]) / N;
    }

    // --- Calculate Theoretical Binomial PMF ---
    std::vector<double> binomial_probabilities(n + 1);
    for (int i = 0; i <= n; ++i) {
        binomial_probabilities[i] = binomialPMF(i, n);
    }

    // --- Calculate Theoretical Normal PDF ---
    // For Bin(n, p) with p=0.5:
    // Mean (mu) = n * p = n * 0.5 = n/2
    // Variance (sigma^2) = n * p * (1-p) = n * 0.5 * 0.5 = n/4
    // Standard Deviation (sigma) = sqrt(n/4) = sqrt(n)/2
    double mean_norm = n / 2.0;
    double stddev_norm = std::sqrt(n) / 2.0;
    std::vector<double> normal_densities(n + 1);
    // Evaluate the normal PDF at each integer slot position
    for (int i = 0; i <= n; ++i) {
        normal_densities[i] = normalPDF(static_cast<double>(i), mean_norm, stddev_norm);
    }

    // --- Print Tabular Results ---
    std::cout << "\n--- Distribution Comparison ---\n";
    std::cout << std::fixed << std::setprecision(8); // Set precision for double output
    std::cout << "Slot | Experimental Freq | Binomial PMF     | Normal PDF      \n";
    std::cout << "--------------------------------------------------------------\n";
    for (int i = 0; i <= n; ++i) {
        std::cout << std::setw(4) << i << " | "
                  << std::setw(17) << experimental_frequencies[i] << " | "
                  << std::setw(16) << binomial_probabilities[i] << " | "
                  << std::setw(15) << normal_densities[i] << "\n";
    }

    // --- Calculate Mean Quadratic Error (MQE) ---
    double mqe_binomial_normal = 0.0;
    double mqe_exp_binomial = 0.0;
    double mqe_exp_normal = 0.0;

    for (int i = 0; i <= n; ++i) {
        // MQE between Binomial PMF and Normal PDF
        double diff_bn = binomial_probabilities[i] - normal_densities[i];
        mqe_binomial_normal += diff_bn * diff_bn;

        // MQE between Experimental Frequencies and Binomial PMF
        double diff_eb = experimental_frequencies[i] - binomial_probabilities[i];
        mqe_exp_binomial += diff_eb * diff_eb;

        // MQE between Experimental Frequencies and Normal PDF
        double diff_en = experimental_frequencies[i] - normal_densities[i];
        mqe_exp_normal += diff_en * diff_en;
    }

    // Divide by the number of points (n+1 slots) to get the mean
    mqe_binomial_normal /= (n + 1);
    mqe_exp_binomial /= (n + 1);
    mqe_exp_normal /= (n + 1);

    std::cout << "\n--- Mean Quadratic Errors (MQE) ---\n";
    std::cout << std::setprecision(10); // Higher precision for MQE values
    std::cout << "MQE (Binomial vs Normal Approximation): " << mqe_binomial_normal << "\n";
    std::cout << "MQE (Experimental vs Binomial):         " << mqe_exp_binomial << "\n";
    std::cout << "MQE (Experimental vs Normal):           " << mqe_exp_normal << "\n";

    // --- Generate PNG Plots ---
    std::string plotDataFilename = "galton_data_n" + std::to_string(n) + "_N" + std::to_string(N) + ".txt";
    std::string binomialNormalPng = "galton_binomial_vs_normal_n" + std::to_string(n) + "_N" + std::to_string(N) + ".png";
    std::string experimentalBinomialPng = "galton_experimental_vs_binomial_n" + std::to_string(n) + "_N" + std::to_string(N) + ".png";

    if (saveDataForPlotting(n, experimental_frequencies, binomial_probabilities, normal_densities, plotDataFilename)) {
        std::string title1 = "Binomial PMF vs Normal PDF (n=" + std::to_string(n) + ", N=" + std::to_string(N) + ")";
        generateBinomialVsNormalPlot(plotDataFilename, binomialNormalPng, title1);

        std::string title2 = "Experimental Freq vs Binomial PMF (n=" + std::to_string(n) + ", N=" + std::to_string(N) + ")";
        generateExperimentalVsBinomialPlot(plotDataFilename, experimentalBinomialPng, title2);
        
        // Optional: Clean up the temporary data file after plotting
        // std::remove(plotDataFilename.c_str());
        // std::cout << "Cleaned up data file: " << plotDataFilename << "\n";
    }
    
    std::cout << "\nSimulation complete. Press Enter to exit.";
    std::cin.ignore(); // Consume any leftover newline character from previous input
    std::cin.get();    // Wait for user to press Enter

    return 0; // Indicate successful execution
}