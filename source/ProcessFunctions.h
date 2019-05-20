/**
 * @file ProcessFunctions.h
 * @author  Al Timofeyev
 * @date    April 4, 2019
 * @brief   A class used to process matrices against Benchmark Functions
 *          and analyze the results.
 */

#ifndef BENCHMARKFUNCTIONS_PROCESSFUNCTIONS_H
#define BENCHMARKFUNCTIONS_PROCESSFUNCTIONS_H

#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include "utilities.h"
#include "BenchmarkFunctions.h"

// -------------- CONSTANTS --------------
/** The default minimum number of dimensions. */
#define DEFAULT_NUM_OF_DIMENSIONS 30
/** The default number of vectors per matrix. */
#define DEFAULT_NUM_OF_VECTORS 30
/** The default minimum boundary for the elements generated. */
#define BOUNDARY_MIN -500.0
/** The default maximum boundary for the elements generated. */
#define BOUNDARY_MAX 500.0

class ProcessFunctions{
private:
    // --------------------------- Structures ---------------------------
    /**
     * @brief  Function Data
     * Function Data Structure, to keep track of all the data
     * used for the Benchmark Functions.
     */
    struct FunctionData
    {
        int functionID;                                     /**< The ID used to determine which of the 18 Benchmark Functions to use.*/
        std::vector<double> fitness;                        /**< The list of fitness for each vector in the matrix.*/
        std::vector<std::vector<double>> functionMatrix;    /**< The matrix of double vectors.*/
        double timeToExecute;                               /**< This is time in ms to process all 30 rows.*/
    };

    /**
     * @brief  Function Analysis
     * Function Analysis Structure, to keep track of the analysis
     * performed on each FunctionData structure. Basically, it compiles
     * and holds the averages of the calculations performed for each function.
     */
    struct FunctionAnalysis
    {
        std::string header = "Function ID,Avg Fitness,Range(min),Range(max),Median,Time(ms)\n"; /**< Header used when saving the data.*/
        std::vector<int> functionIDs;               /**< List of function IDs.*/
        std::vector<double> avgFuntionFitness;      /**< List of the average fitness per FunctionData structure.*/
        std::vector<std::vector<double>> ranges;    /**< List of ranges for each fitness result in resultsOfFunctions.*/
        std::vector<double> medianFunctionFitness;  /**< List of the Median fitness from each FunctionData structure.*/
        std::vector<double> processTimes;           /**< List of process times in ms for all functions.*/
    };

    // --------------------------- Variables ----------------------------
    int numOfDimensions;
    std::vector<FunctionData> resultsOfFunctions;
    FunctionAnalysis analysis;

    // --------------------- Functions Declarations ---------------------
    FunctionData generateMatrix(double minBoundary, double maxBoundary);

    double calculateFitnessOfVector(std::vector<double> &vect, int functionID);    /**< Calculates the fitness of a single vector.*/
    void calculateFitnessOfMatrix(FunctionData &data);                             /**< Calculates the fitness of all vectors in matrix.*/
    void analyzeFunctionResults(FunctionData &data);                               /**< Analyzes the results of the functions.*/
    double calculateAvgFitness(FunctionData &data);
    double getMinFitness(FunctionData &data);
    double getMaxFitness(FunctionData &data);

    void saveFunctionMatrixToFile(std::string filename, FunctionData &data);    /**< Saves the matrix to file.*/
    void saveAllFunctionDataToFile(std::string filename, FunctionData &data);   /**< Saves the results of the function and it's data to file.*/

    void quicksort(FunctionData &data, int L, int R);
    void swap(FunctionData &data, int x, int y);

public:
    // --------------------- Constructor Declarations ---------------------
    ProcessFunctions(); // Sets the number of dimensions to 0;

    // --------------------- Functions Declarations ---------------------
    void setNumOfDimensions(int dimensions);    /**< Sets the number of dimensions.*/
    int getNumOfDimensions();                   /**< Returns the number of dimensions.*/

    void constructMatrix();                                                      /**< Uses all default constants, or previously user-set dimensions.*/
    void constructMatrix(int funcID, double minBoundary, double maxBoundary);    /**< Uses default number of dimensions.*/

    void calculateFitnessOfAllMatrices();   /**< Calculates Fitness for all matrices in resultsOfFunctions.*/
    void analyzeAllFunctionResults();       /**< Analyzes all the results from resultsOfFunctions.*/

    void saveAllMatricesToFile();       /**< Saves all the matrices in resultsOfFunctions to files.*/
    void saveAllProcessedDataToFile();  /**< Saves all the data in resultsOfFunctions to files.*/
    void saveAllAnalyzedDataToFile();   /**< Saves all analyzed data in analysis to file.*/

    void printAllFunctionIDs();             /**< Prints all the possible Function IDs to the screen.*/
    void printFunctionResults();            /**< Prints all the FunctionData structures in resultsOfFunctions.*/
    void printFunctionResultsAnalysis();    /**< Prints all the Analysis Results in analysis.*/
};


#endif //BENCHMARKFUNCTIONS_PROCESSFUNCTIONS_H
