/**
 * @file utilities.cpp
 * @author  Al Timofeyev
 * @date    March 28, 2019
 * @brief   This utilities file is used to create matricies using the
 *          Mersenne Twister and store them in Excel files.
 */

#include "utilities.h"

// Utility Functions.
// --------------------------------------------------------------------------------------------------------------------
/**
 * @brief Parses a string of numbers into a vector of doubles.
 *
 * Constructs and returns a vector of doubles, given a string list of
 * numbers and a delimiter.
 *
 * @note The input string list MUST be a list of doubles!
 *
 * @param list  A string list of numbers.
 * @param delimiter A character used to separate the numbers in the string list.
 *
 * @return Returns a vector filled with doubles that were extracted from the string list.
 */
std::vector<double> parseString(std::string list, char delimiter)
{
    std::vector<double> numList;

    int startIndex = 0;
    int endIndex = list.find_first_of(delimiter, startIndex);

    while(endIndex < list.size() && endIndex != std::string::npos)
    {
        // Only convert to int and store in numList if start/end index aren't equal.
        if(startIndex != endIndex)
            numList.push_back(stod(list.substr(startIndex, endIndex - startIndex)));

        // Update the start and end indexes.
        startIndex = endIndex + 1;
        endIndex = list.find_first_of(delimiter, startIndex);
    }

    // Add the last element to the numList.
    if(startIndex < list.size())
        numList.push_back(stod(list.substr(startIndex, endIndex - startIndex)));

    // Return the vector of doubles.
    return numList;
}