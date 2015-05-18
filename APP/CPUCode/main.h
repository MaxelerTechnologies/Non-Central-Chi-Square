//*********************************************************************//
// Non-Central Chi Square Random Number Generator
//
// Author:  Sébastien Racanière
//
// Date:    26 Jun 2012
//
//*********************************************************************//

#ifndef MAIN_H_
#define MAIN_H_

#include <vector>
#include <ctime>
#include <stdint.h>


/**
 * \brief Run CPU version of the program
 * \param [in] degree k, degree of the distribution
 * \param [in] outputCount Number of values to sample from distribution
 * \param [in] lambda descentralizing constant
 * \param [out] time Time taken to compute result, in seconds
 * \return vector of resulting random samples, caller has the responsability
 *                to delete() it
 */
double* runCPU(unsigned int degree,
        unsigned int outputCount,
        double lambda,
        double& runtime);


/**
 * \brief Accelerated, DFE version of \ref runCPU
 * \param [in] degree k, degree of the distribution
 * \param [in] outputCount Number of values to sample from distribution
 * \param [in] lambda descentralizing constant
 * \param [out] time Time taken to compute result, in seconds
 * \return vector of resulting random samples, caller has the responsability
 *                to delete() it
 */
double* runDFE(unsigned int degree,
        unsigned int outputCount,
        double lambda,
        double& runtime);


/**
 * \brief Mersenne Twister pseudo-random number generator.
 * Generate 32 bits numbers with uniform distribution.
 */
class MersenneTwister {
public:
    /**
     * \brief Constructor.
     * \param [in] seed initial value defining the sequence of values generated
     */
    MersenneTwister(uint32_t seed);

    /**
     * \return tempered pseudo-random number based on the index-th value,
     *         invokes updateSeeds() every DEGREES calls
     */
    uint32_t get();
private:
    /** 624 degrees of recurrence for this implementation */
    static const size_t DEGREES = 624;
    /** array containing the state of the generator */
    uint32_t MT[DEGREES];
    /** current state index */
    int index;

    /**
     * Initialize internal state array.
     * @param seed basis for internal state initial values
     */
    void initSeeds(uint32_t seed);

    /**
     * Update internal state array.
     */
    void updateSeeds();
};

/**
 * \brief Gaussian pseudo-random number generator.
 * Generate double-precision floating-point numbers with normal distribution.
 */
class Gaussian {
public:
    /**
     * \brief Constructor
     * \param [in] mean mu, the mean value of the normal distribution
     * \param [in] stdDeviation sigma, standard deviation of the distribution
     * \param [in] seed initial value defining the sequence of values generated
     */
    Gaussian(double mean, double stdDeviation, uint32_t seed);

    /**
     * @return a pseudo-random number
     */
    double get();

private:
    /** mean */
    double mu;
    /** standard deviation */
    double sigma;
    /**
     * cached value -- values are always generated in pairs, one is returned
     * directly by get(), the other is stored for the next call
     */
    double z;
    /** flag indicating whether cached value z is valid */
    bool hasCache;
    /** internal uniform-distribution pseudo-random number generator */
    MersenneTwister mt;
};

/**
 * \brief Central Chi-Square pseudo-random number generator.
 * Generate double-precision floating-point numbers following the distribution:
 *   X(k) = sum_over_k(Zi^2)
 * Where Z follows the Gaussian distribution.
 */
class CentralChiSquare {
public:
    /**
     * \brief Constructor
     * \param [in] k degree of the distribution
     * \param [in] seed initial value defining the sequence of values generated
     */
    CentralChiSquare(unsigned int k, uint32_t seed);

    /**
     * \return a pseudo-random number
     */
    double get();

private:
    /** whether the number of degrees of the distribution is odd */
    bool isOdd;
    /** internal Gaussian pseudo-random number generator */
    Gaussian g;
    /** internal array of uniform-distribution pseudo-random number generators */
    std::vector<MersenneTwister> mt;
};

/**
 * \brief Non-central Chi-Square pseudo-random number generator.
 * Generate double-precision floating-point numbers following the distribution:
 *   X(k, l) = sum_over_k-1(Zi^2) + (Z + l^0.5)^2
 * Where Z follows the Gaussian distribution.
 */
class NonCentralChiSquareCPU {
public:
    /**
     * \brief Constructor
     * \param [in] k degree of the distribution
     * \param [in] lambda descentralizing constant
     * \param [in] seed initial value defining the sequence of values generated
     */
    NonCentralChiSquareCPU(unsigned int k, double lambda, uint32_t seed);
    double get();

private:
    /** degree of the distribution */
    unsigned int k;
    /** square root of lambda */
    double sqrtLambda;
    /** internal Central Chi-Square pseudo-random number generator */
    CentralChiSquare ccs;
    /** internal Gaussian pseudo-random number generator */
    Gaussian g;
};

/**
 * \brief Utility class to measure time
 */
class Clock {
public:
    Clock();
    void start();
    void stop();
    long long getNanoSeconds();

private:
    struct timespec ts;
    long long total;
};

/**
 * \brief Output to stdout a report comparing results.
 * \param [in] resultCount number of values in each result array
 * \param [in] cpuResults results yielded by the CPU version
 * \param [in] dfeResults results yielded by the DFE version
 */
void printResultComparisonReport(size_t resultCount,
        const double* cpuResults,
        const double* dfeResults);

/**
 * \brief Output results to file.
 * \param [in] fileName path to output file
 * \param [in] resultCount number of values in each result array
 * \param [in] cpuResults results yielded by the CPU version
 * \param [in] dfeResults results yielded by the DFE version
 */
void outputResults(const char* fileName,
        size_t resultCount,
        const double* cpuResults,
        const double* dfeResults);

/**
 * \brief Utility function to determine the least count for an array of items
 *        that is byte-aligned for PCIE communication with the DFE.
 */
template<typename T>
static unsigned int getPciePaddedCount(unsigned int count);


/**
 * \brief Warm-up function to guarantee correct timing measurements.
 */
void warmup();


#endif /* MAIN_H_ */
