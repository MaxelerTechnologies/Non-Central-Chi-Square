//*********************************************************************//
// Non-Central Chi Square Random Number Generator
//
// Author:  S�bastien Racani�re
//
// Date:    26 Jun 2012
//
//*********************************************************************//

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <limits>
#include <stdio.h>

#include "cmdline.h"
#include "main.h"

#include "MaxSLiCInterface.h"
#include "NCChiSquare.h"



int main(int argc, char **argv) {
    struct gengetopt_args_info argsInfo;

    if (cmdline_parser(argc, argv, &argsInfo) != 0) {
        std::cout << "Try `" << argv[0] <<
                " --help` for more information." << std::endl;
        exit(EXIT_FAILURE);
    }
    if(argsInfo.degree_arg <= 0 ||
            argsInfo.degree_arg > NCChiSquare_maximumDegree) {
        std::cout << "degree must be an integer in (0," <<
                NCChiSquare_maximumDegree << "], not " << argsInfo.degree_arg << ".\nTry `" << argv[0] <<
                " --help` for more information." << std::endl;
        exit(EXIT_FAILURE);
    }

    double* cpuRandomValues = NULL;
    double* dfeRandomValues = NULL;

    // CPU version
    double cpuTime = 0.0;
    if(!argsInfo.dfe_given) {
        std::cout << "\nRunning CPU version..." << std::endl;
        cpuRandomValues = runCPU(argsInfo.degree_arg,
                argsInfo.output_count_arg,
                argsInfo.lambda_arg,
                cpuTime);
        std::cout << "Runtime: " << (cpuTime*1e3) << "ms" << std::endl;
    }

    // DFE version
    double dfeTime = 0.0;
    if(!argsInfo.cpu_given) {
        std::cout << "\nRunning accelerated DFE version..." << std::endl;
        warmup();
        dfeRandomValues = runDFE(argsInfo.degree_arg,
                argsInfo.output_count_arg,
                argsInfo.lambda_arg,
                dfeTime);
        std::cout << "Runtime: " << (dfeTime*1e3) << "ms" << std::endl;
    }

    if(!argsInfo.dfe_given && !argsInfo.cpu_given) {
        std::cout << "\nMeasured speedup: " << (cpuTime/dfeTime) << std::endl;
        printResultComparisonReport(argsInfo.output_count_arg,
                cpuRandomValues, dfeRandomValues);
    }

    if(argsInfo.file_output_given) {
        outputResults(argsInfo.file_output_arg, argsInfo.output_count_arg,
                cpuRandomValues, dfeRandomValues);
    }

    cmdline_parser_free(&argsInfo);
    if(NULL!=cpuRandomValues) {
        delete[] cpuRandomValues;
    }
    if(NULL!=dfeRandomValues) {
        delete[] dfeRandomValues;
    }

    return 0;
}


double* runCPU(unsigned int degree,
        unsigned int outputCount,
        double lambda,
        double& runtime) {

    double* randomValues = new double[outputCount];

    Clock clock;
    clock.start();

    NonCentralChiSquareCPU nccs(degree, lambda, 13);
    for(unsigned int i = 0; i < outputCount; ++i) {
        randomValues[i] = nccs.get();
    }

    clock.stop();
    runtime = clock.getNanoSeconds()/1.0e9;

    return randomValues;
}


double* runDFE(unsigned int degree,
        unsigned int outputCount,
        double lambda,
        double& runtime) {

    outputCount = getPciePaddedCount<double>(outputCount);
    double* randomValues = new double[outputCount];

    Clock clock;
    clock.start();

    const double sqrtLambda = sqrt(lambda);
    NCChiSquare(degree, outputCount, sqrtLambda, randomValues);

    clock.stop();
    runtime = clock.getNanoSeconds()/1.0e9;

    return randomValues;
}


template<typename T>
static unsigned int getPciePaddedCount(unsigned int count) {
    // byte-alignment for PCIE communication
    const unsigned int pcieAlignment = 16;
    // number of bytes in given count
    long unsigned int byteCount = count * sizeof(T);

    // could be done with least common multiple starting from count, but fact
    // is that the iterative method is simpler and very probably more than
    // adequate for the domain of values for which this function is expected
    // to be used
    for(; byteCount%pcieAlignment!=0; byteCount += sizeof(T)) {
    }

    return byteCount/sizeof(T);
}


void printResultComparisonReport(const size_t resultCount,
        const double* const cpuResults,
        const double* const dfeResults) {

    double accumDiff=0.0;
    double maxDiff=-1.0;
    double accumCpu=0.0;
    double accumDfe=0.0;

    for(size_t i=0; i<resultCount; ++i) {
        double absDiff = fabs(cpuResults[i] - dfeResults[i]);
        maxDiff = std::max(maxDiff, absDiff);
        accumCpu += cpuResults[i];
        accumDfe += dfeResults[i];
    }

    std::cout.precision(3);
    std::cout << std::scientific << "\nOutput values comparison:\n"
            "\tCPU average:\t\t\t" << accumCpu/resultCount << "\n" <<
            "\tDFE average:\t\t\t" << accumDfe/resultCount << "\n" <<
            "\tMaximum absolute difference:\t" << maxDiff << "\n" <<
            "\tAverage absolute difference:\t" << accumDiff/resultCount <<
            std::endl;
}


void outputResults(const char* const fileName,
        size_t resultCount,
        const double* const cpuResults,
        const double* const dfeResults) {

    std::ofstream fout(fileName);
    if(!fout) {
        std::cerr << strerror(errno) << std::endl;
        exit(EXIT_FAILURE);
    }
    fout.precision(std::numeric_limits<double>::digits10);
    fout << std::scientific;

    std::string separator = "";
    if(NULL!=cpuResults && NULL!=dfeResults) {
        separator = ",";
    }

    if(NULL!=cpuResults) {
        fout << "CPU";
    }
    fout << separator;
    if(NULL!=dfeResults) {
        fout << "DFE";
    }
    fout << "\n";

    for(size_t i=0; i<resultCount; ++i) {
        if(NULL!=cpuResults) {
            fout << cpuResults[i];
        }
        fout << separator;
        if(NULL!=dfeResults) {
            fout << dfeResults[i];
        }
        fout << "\n";
    }
    fout.close();
}


MersenneTwister::MersenneTwister(uint32_t seed) : index(0) {
    initSeeds(seed);
}


void MersenneTwister::initSeeds(const uint32_t seed) {
    const uint32_t init_key[] = { 0x123 + seed, 0x234 + seed,
            0x345 + seed, 0x456 + seed };
    const size_t key_length = sizeof(init_key)/sizeof(*init_key);
    static const size_t N = DEGREES;

    MT[0] = 19650218L;
    for (size_t mti = 1; mti < N; mti++) {
        MT[mti] = (1812433253UL * (MT[mti - 1] ^ (MT[mti - 1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect */
        /* only MSBs of the array mt[]. */
        /* 2002/01/09 modified by Makoto Matsumoto */
        MT[mti] &= 0xffffffffL;
        /* for >32 bit machines */
    }

    size_t i = 1;
    size_t j = 0;
    for (ssize_t k = std::max(N, key_length); k > 0; k--) {
        MT[i] = (MT[i] ^ ((MT[i - 1] ^ (MT[i - 1] >> 30)) * 1664525UL))
                + init_key[j] + j; /* non linear */
        MT[i] &= 0xffffffffL; /* for WORDSIZE > 32 machines */
        i++;
        j++;
        if (i >= N) {
            MT[0] = MT[N - 1];
            i = 1;
        }
        if (j >= key_length)
            j = 0;
    }
    for (ssize_t k = N - 1; k > 0; k--) {
        MT[i] = (MT[i] ^ ((MT[i - 1] ^ (MT[i - 1] >> 30)) * 1566083941UL))
                - i; /* non linear */
        MT[i] &= 0xffffffffL; /* for WORDSIZE > 32 machines */
        i++;
        if (i >= N) {
            MT[0] = MT[N - 1];
            i = 1;
        }
    }

    MT[0] = 0x80000000L; /* MSB is 1; assuring non-zero initial array */
}


inline void MersenneTwister::updateSeeds() {
    // Co-efficient of rational normal form matrix
    const static uint32_t COEFFICIENT = 0x9908b0dfU;

    for(size_t i = 0; i < DEGREES; ++i) {
        uint32_t y = (MT[i] & 0x80000000U) +
                (MT[(i+1) % DEGREES] & 0x7FFFFFFFU);
        MT[i] = MT[(i + 397) % DEGREES] ^ (y >> 1);
        if(y % 2 != 0) { // y is odd
            MT[i] = MT[i] ^ COEFFICIENT;
        }
    }
}


inline uint32_t MersenneTwister::get() {
    // TGFSR(R) tempering bitmasks
    const static uint32_t MASK1 = 0x9d2c5680U;
    const static uint32_t MASK2 = 0xefc60000U;
    // TGFSR(R) tempering bit shifts
    const static uint8_t TGFSR_SHIFT1 = 7;
    const static uint8_t TGFSR_SHIFT2 = 15;
    // Additional Mersenne Twister tempering bit shifts
    const static uint8_t MTADD_SHIFT1 = 11;
    const static uint8_t MTADD_SHIFT2 = 18;

    uint32_t y = MT[index];
    y = y ^  (y >> MTADD_SHIFT1);
    y = y ^ ((y << TGFSR_SHIFT1) & MASK1);
    y = y ^ ((y << TGFSR_SHIFT2) & MASK2);
    y = y ^  (y >> MTADD_SHIFT2);

    index = (index + 1) % DEGREES;
    if(index == 0) {
        updateSeeds();
    }
    return y;
}


Gaussian::Gaussian(double mean, double stdDeviation, uint32_t seed) :
        mu(mean),
        sigma(stdDeviation),
        z(0.0),
        hasCache(false),
        mt(seed) {
}


inline double Gaussian::get() {
    if(hasCache) {
        hasCache = false;
        return z;
    }
    const double a = mt.get()/double(1L << 32);
    const double b = mt.get()/double(1L << 32);
    const double R = sqrt(-2.0*log(a));
    const double c = cos(2.0*M_PI*b);
    const double s = sin(2.0*M_PI*b);
    z = (R * s * sigma + mu);
    hasCache = true;

    return (R * c * sigma + mu);
}


CentralChiSquare::CentralChiSquare(unsigned int k, uint32_t seed) :
        isOdd(k%2!=0),
        g(0.0, 1.0, seed) {

    static const uint32_t SEED_TEMPERING = 397432319U;
    const size_t twisterCount = k/2;

    if(twisterCount > 0) {
        mt.reserve(twisterCount);
    }

    for(unsigned int i = 0; i < twisterCount; ++i) {
        mt.push_back(MersenneTwister(seed + (i+1) * SEED_TEMPERING));
    }
}


inline double CentralChiSquare::get() {
    double z = 0.0;
    if(isOdd) {
        const double x = g.get();
        z += x*x;
    }
    for(std::vector<MersenneTwister>::iterator it = mt.begin();
            it!=mt.end(); ++it) {

        const double x = it->get()/double(1L << 32);
        z += -2.0*log(x);
    }

    return z;
}


NonCentralChiSquareCPU::NonCentralChiSquareCPU(unsigned int k,
        double lambda,
        uint32_t seed) :
        k(k),
        sqrtLambda(sqrt(lambda)),
        ccs(k-1, seed),
        g(0.0, 1.0, seed) {
}


inline double NonCentralChiSquareCPU::get() {
    const double a = g.get() + sqrtLambda;
    const double b = ccs.get();

    return a*a + b;
}


Clock::Clock() : total(0L) {
}


inline void Clock::start() {
    clock_gettime(CLOCK_REALTIME, &ts);
}


inline void Clock::stop() {
    struct timespec te;
    clock_gettime(CLOCK_REALTIME, &te);
    long long diff = ((long long)te.tv_sec * (long long)1e9 + te.tv_nsec)
         - ((long long)ts.tv_sec * (long long)1e9 + ts.tv_nsec);
    total += diff;
}


long long Clock::getNanoSeconds() {
    return total;
}


void warmup() {
    double output[2];
    NCChiSquare(0, sizeof(output)/sizeof(*output), 1.0, output);
}
