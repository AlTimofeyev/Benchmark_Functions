/**
 * @file BenchmarkFunctions.h
 * @author  Al Timofeyev
 * @date    March 28, 2019
 * @brief   A library of benchmark functions.
 */

#ifndef BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H
#define BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H



//#define _USE_MATH_DEFINES // Uncomment if cmath constants are desirable, like M_PI.

#include <vector>
#include <math.h>   // Sine and Cosine, square root
#include <cmath>    // Absolute value of doubles

/** The pi constant used for calculations. */
const double pi = 3.14159265358979323846;

//using namespace std;

double schefelsFunc(std::vector<double> &vect, int size);
double deJongsFunc(std::vector<double> &vect, int size);
double rosenbrockFunc(std::vector<double> &vect, int size);
double rastriginFunc(std::vector<double> &vect, int size);
double griewangkFunc(std::vector<double> &vect, int size);
double sineEnvelopeSineWaveFunc(std::vector<double> &vect, int size);
double stretchedVSineWaveFunc(std::vector<double> &vect, int size);
double ackleysOneFunc(std::vector<double> &vect, int size);
double ackleysTwoFunc(std::vector<double> &vect, int size);
double eggHolderFunc(std::vector<double> &vect, int size);
double ranaFunc(std::vector<double> &vect, int size);
double pathologicalFunc(std::vector<double> &vect, int size);
double michalewiczFunc(std::vector<double> &vect, int size);
double mastersCosWaveFunc(std::vector<double> &vect, int size);
double quarticFunc(std::vector<double> &vect, int size);
double levyFunc(std::vector<double> &vect, int size);
double stepFunc(std::vector<double> &vect, int size);
double alpineFunc(std::vector<double> &vect, int size);

#endif //BENCHMARKFUNCTIONS_BENCHMARKFUNCTIONS_H