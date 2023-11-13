#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

class AsianCallOption {
private:
    double spotPrice;
    double strikePrice;
    double riskFreeRate;
    double volatility;
    double timeToMaturity;
    int numPaths;
    int numSteps;

public:
    AsianCallOption(double S, double K, double r, double sigma, double T, int paths, int steps)
        : spotPrice(S), strikePrice(K), riskFreeRate(r), volatility(sigma), timeToMaturity(T), numPaths(paths), numSteps(steps) {}

    double rand01() {
        return rand() / (double)RAND_MAX;
    }
    // Function to generate a vector of random numbers
    std::vector<double> generateRandomNumbers(int count) {
        std::vector<double> randomNumbers;
        for (int i = 0; i < count; ++i) {
            randomNumbers.push_back(rand01());
        }
        return randomNumbers;
    }

    // Monte Carlo simulation for Asian call option pricing using antithetic variance reduction
    double price() {
        double dt = timeToMaturity / numSteps;
        double sumPayoff = 0.0;
        srand(static_cast<unsigned>(time(0)));  // Seed for random number generation

        for (int i = 0; i < numPaths; ++i) {
            // Generate two sets of independent random numbers using antithetic variates
            std::vector<double> epsilon1 = generateRandomNumbers(numSteps);
            std::vector<double> epsilon2 = epsilon1;  // Antithetic variates
            for (double& e : epsilon2) {
                e = 1.0 - e;
            }

            double spotPriceCopy = spotPrice;
            double sumSpotPrice1 = spotPrice;
            double sumSpotPrice2 = spotPrice;

            for (int j = 1; j < numSteps; ++j) {
                double diffusion1 = volatility * sqrt(dt) * epsilon1[j - 1];
                double diffusion2 = volatility * sqrt(dt) * epsilon2[j - 1];

                // Apply antithetic variates to the diffusion term
                spotPriceCopy *= exp((riskFreeRate - 0.5 * volatility * volatility) * dt + 0.5 * (diffusion1 + diffusion2));
                sumSpotPrice1 += spotPriceCopy;
                sumSpotPrice2 += spotPriceCopy;
            }

            double averageSpotPrice1 = sumSpotPrice1 / numSteps;
            double averageSpotPrice2 = sumSpotPrice2 / numSteps;

            double payoff1 = fmax(averageSpotPrice1 - strikePrice, 0.0);
            double payoff2 = fmax(averageSpotPrice2 - strikePrice, 0.0);

            // Apply antithetic variates to the payoff
            double payoff = 0.5 * (payoff1 + payoff2);

            sumPayoff += payoff;
        }

        return exp(-riskFreeRate * timeToMaturity) * (sumPayoff / numPaths);
    }
};

int main() {
    // Example parameters
    double spotPrice = 200.0;
    double strikePrice = 220.0;
    double riskFreeRate = 0.06;
    double volatility = 0.1;
    double timeToMaturity = 2.0;
    int numPaths = 200000;
    int numSteps = 252;

    // Create an instance of the AsianCallOption class
    AsianCallOption asianOption(spotPrice, strikePrice, riskFreeRate, volatility, timeToMaturity, numPaths, numSteps);

    // Calculate Asian call option price using the price() method
    double optionPrice = asianOption.price();

    // Display the result
    std::cout << "Asian Call Option Price: " << optionPrice << std::endl;

    return 0;
}
