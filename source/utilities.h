/**
 * @file utilities.h
 * @author  Al Timofeyev
 * @date    March 28, 2019
 * @brief   This utilities file is used to create matricies using the
 *          Mersenne Twister and store them in Excel files.
 */

#ifndef BENCHMARKFUNCTIONS_UTILITIES_H
#define BENCHMARKFUNCTIONS_UTILITIES_H

#include <string>
#include <vector>

std::vector<double> parseString(std::string list, char delimiter); // Parses a string of numbers into a vector of doubles.

#endif //BENCHMARKFUNCTIONS_UTILITIES_H
