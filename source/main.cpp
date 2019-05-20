/*
 * Author:  Al Timofeyev
 * Date:    3/27/2019
 * Desc:    Given 18 benchmark functions, create matrices for each one
 *          and test them. Output the results and choose the best fitness.
 *          Each function performs calculations on a matrix of 30 pseudo-random
 *          10, 20, or 30 multidimensional vectors. So, a matrix of 30 vectors
 *          that have either 10, 20, or 30 dimensions (elements) per vector.
 */

#include <iostream>
#include "ProcessFunctions.h"

using namespace std;

int main()
{
    ProcessFunctions procFuncs;
    procFuncs.setNumOfDimensions(30);

    procFuncs.constructMatrix(1, -512, 512);
    procFuncs.constructMatrix(2, -100, 100);
    procFuncs.constructMatrix(3, -100, 100);
    procFuncs.constructMatrix(4, -30, 30);
    procFuncs.constructMatrix(5, -500, 500);
    procFuncs.constructMatrix(6, -30, 30);
    procFuncs.constructMatrix(7, -30, 30);
    procFuncs.constructMatrix(8, -32, 32);
    procFuncs.constructMatrix(9, -32, 32);
    procFuncs.constructMatrix(10, -500, 500);
    procFuncs.constructMatrix(11, -500, 500);
    procFuncs.constructMatrix(12, -100, 100);
    procFuncs.constructMatrix(13, 0, pi);
    procFuncs.constructMatrix(14, -30, 30);
    procFuncs.constructMatrix(15, -100, 100);
    procFuncs.constructMatrix(16, -10, 10);
    procFuncs.constructMatrix(17, -100, 100);
    procFuncs.constructMatrix(18, -100, 100);


    procFuncs.calculateFitnessOfAllMatrices();
    procFuncs.saveAllProcessedDataToFile();

    procFuncs.analyzeAllFunctionResults();
    procFuncs.saveAllAnalyzedDataToFile();
}
