#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

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

    // Function to generate a random number between 0 and 1
    double rand01() {
        return rand() / (double)RAND_MAX;
    }

    // Monte Carlo simulation for Asian call option pricing
    double price() {
        double dt = timeToMaturity / numSteps;
        double sumPayoff = 0.0;
        srand(static_cast<unsigned>(time(0)));  // Seed for random number generation

        for (int i = 0; i < numPaths; ++i) {
            double spotPriceCopy = spotPrice;
            double sumSpotPrice = spotPrice;

            for (int j = 1; j < numSteps; ++j) {
                double epsilon = rand01();
                double drift = riskFreeRate - 0.5 * volatility * volatility;
                double diffusion = volatility * sqrt(dt) * epsilon;

                spotPriceCopy *= exp(drift * dt + diffusion);
                sumSpotPrice += spotPriceCopy;
            }

            double averageSpotPrice = sumSpotPrice / numSteps;
            double payoff = fmax(averageSpotPrice - strikePrice, 0.0);  // Asian call option payoff

            sumPayoff += payoff;
        }

        return exp(-riskFreeRate * timeToMaturity) * (sumPayoff / numPaths);
    }
};

int main() {
    // Example parameters
    double spotPrice = 100.0;
    double strikePrice = 100.0;
    double riskFreeRate = 0.05;
    double volatility = 0.2;
    double timeToMaturity = 1.0;
    int numPaths = 100000;
    int numSteps = 252;

    // Create an instance of the AsianCallOption class
    AsianCallOption asianOption(spotPrice, strikePrice, riskFreeRate, volatility, timeToMaturity, numPaths, numSteps);

    // Calculate Asian call option price using the price() method
    double optionPrice = asianOption.price();

    // Display the result
    std::cout << "Asian Call Option Price: " << optionPrice << std::endl;

    return 0;
}
