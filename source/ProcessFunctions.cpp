/**
 * @file ProcessFunctions.cpp
 * @class ProcessFunctions ProcessFunctions.h "ProcessFunctions.h"
 * @author  Al Timofeyev
 * @date    April 4, 2019
 * @brief   A class used to process matrices against Benchmark Functions
 *          and analyze the results.
 */

#include "FilenameConstants.h"
#include "ProcessFunctions.h"

// ----------------------------------------------
// ---------------- CONSTRUCTORS ----------------
// ----------------------------------------------
/**
 * @brief The default constructor for the ProcessFunctions class.
 * The default constructor only initializes the numOfDimensions
 * variable to 0;
 */
ProcessFunctions::ProcessFunctions()
{
    numOfDimensions = 0;
}

// -------------------------------------------------------------------------------------------
// --------------------------------- PUBLIC FUNCTIONS BELOW ----------------------------------
// -------------------------------------------------------------------------------------------
/**
 * @brief Sets the number of dimensions for the ProcessFunctions object.
 *
 * After setting the new number of dimensions, the resultsOfFunctions vector
 * that held all the previous data, for the previous number of dimensions, is
 * also reset to 0.
 *
 * @param dimensions    The number of dimensions in the matrix data
 *                      (dimensions = size of each vector in the matrix).
 */
void ProcessFunctions::setNumOfDimensions(int dimensions)
{
    numOfDimensions = dimensions;
    resultsOfFunctions.resize(0);
}

/**
 * @brief Returns the number of dimensions used for the matrix.
 * @return The value stored in the numOfDimensions variable.
 */
int ProcessFunctions::getNumOfDimensions()
{
    return numOfDimensions;
}

/**
 * @brief Generates a matrix using Mersenne Twister.
 *
 * A matrix is constructed using the default number of dimensions, or a
 * previously user-set number of dimensions, and the default minimum and
 * maximum bound. Saves the constructed matrix to variable resultsOfFunctions.
 */
void ProcessFunctions::constructMatrix()
{
    // If the number of dimensions is 0, set it to the default value.
    if(numOfDimensions == 0)
        setNumOfDimensions(DEFAULT_NUM_OF_DIMENSIONS);

    // Create a Mersenne Twister pseudo-random number generator.
    // Generate a random function ID.
    std::mt19937 randGenerator(time(NULL));
    std::uniform_int_distribution<int> dis(1, 18);
    int funcID = dis(randGenerator);

    // Construct a matrix with user-provided boundaries.
    FunctionData funcData = generateMatrix(BOUNDARY_MIN, BOUNDARY_MAX);
    funcData.functionID = funcID;

    // Save the constructed matrix to resultsOfFunctions vector.
    resultsOfFunctions.push_back(funcData);
}

/**
 * @brief Generates a matrix using Mersenne Twister.
 *
 * A matrix is constructed using the default value of 30 dimensions, or a
 * previously user-set number of dimensions, and a user-provided minimum
 * and maximum bound. Saves the constructed matrix to variable resultsOfFunctions.
 *
 * @param funcID The function ID for which Benchmark Function the matrix is generated for.
 * @param minBoundary, maxBoundary  The minimum and maximum boundaries for the values
 *                                  in the matrix.
 */
void ProcessFunctions::constructMatrix(int funcID, double minBoundary, double maxBoundary)
{
    // If the function ID is out of range, notify user of Function IDs and exit.
    if(funcID < 1 || funcID > 18)
    {
        std::cout <<"\n******* ";
        std::cout << "Cannot generate matrix for Function ID " << funcID;
        std::cout <<" *******";
        printAllFunctionIDs();
        return;
    }

    // If the number of dimensions is 0, set it to the default value.
    if(numOfDimensions == 0)
        setNumOfDimensions(DEFAULT_NUM_OF_DIMENSIONS);

    // Construct a matrix with user-provided boundaries.
    FunctionData funcData = generateMatrix(minBoundary, maxBoundary);
    funcData.functionID = funcID;

    // Save the constructed matrix to resultsOfFunctions vector.
    resultsOfFunctions.push_back(funcData);
}

/**
 * @brief Calculates the fitness of all Matrices in resultsOfFunctions vector.
 */
void ProcessFunctions::calculateFitnessOfAllMatrices()
{
    for(int numOfData = 0; numOfData < resultsOfFunctions.size(); numOfData++)
    {
        // Record the start and end time of executing the benchmark function on matrix.
        auto startTime = std::chrono::high_resolution_clock::now();
        calculateFitnessOfMatrix(resultsOfFunctions[numOfData]);
        auto endTime = std::chrono::high_resolution_clock::now();

        // Calculate elapsed time in milliseconds it took to execute the benchmark function.
        auto elapsedTime = endTime - startTime;
        double elapsedTimeMS = std::chrono::duration_cast<std::chrono::milliseconds>(elapsedTime).count();

        // Save elapsed time to the FunctionData struct.
        resultsOfFunctions[numOfData].timeToExecute = elapsedTimeMS;

        // Sort the results by fitness.
        quicksort(resultsOfFunctions[numOfData], 0, resultsOfFunctions[numOfData].fitness.size()-1);
    }
}

void ProcessFunctions::analyzeAllFunctionResults()
{
    for(int numOfData = 0; numOfData < resultsOfFunctions.size(); numOfData++)
        analyzeFunctionResults(resultsOfFunctions[numOfData]);
}

/**
 * @brief Calculates the average fitness.
 * @param data The FunctionData structure that contains the list of fitness values.
 * @return The average fitness from the list of fitness values.
 */
double ProcessFunctions::calculateAvgFitness(FunctionData &data)
{
    double average;
    double summedUp = 0;

    // Sum up all the fitness values.
    for(int row = 0; row < data.fitness.size(); row++)
        summedUp += data.fitness[row];

    // Calculate average.
    average = summedUp / data.fitness.size();

    // Return the average.
    return average;
}

/**
 * @brief Returns the minimum fitness of the data in FunctionData struct.
 * @param data The FunctionData structure that contains a list of fitness values.
 * @return The Minimum fitness in FunctionaData data structure.
 */
double ProcessFunctions::getMinFitness(FunctionData &data)
{
    return data.fitness[0];
}

/**
 * @brief Returns the maximum fitness of the data in FunctionData struct.
 * @param data The FunctionData structure that contains a list of fitness values.
 * @return The Maximum fitness in FunctionaData data structure.
 */
double ProcessFunctions::getMaxFitness(FunctionData &data)
{
    return data.fitness[data.fitness.size() - 1];
}

/**
 * @brief Saves all the matrices in resultsOfFunctions vector to files.
 */
void ProcessFunctions::saveAllMatricesToFile()
{
    for(int numOfResults = 0; numOfResults < resultsOfFunctions.size(); numOfResults++)
    {
        FunctionData data = resultsOfFunctions[numOfResults];
        int funcID = data.functionID;

        // Save with filename referenced by function IDs.
        switch(funcID)
        {
            case 1:
                saveFunctionMatrixToFile(in_schefelsFilename, data);
                break;
            case 2:
                saveFunctionMatrixToFile(in_deJongsFilename, data);
                break;
            case 3:
                saveFunctionMatrixToFile(in_rosenbrockFilename, data);
                break;
            case 4:
                saveFunctionMatrixToFile(in_rastriginFilename, data);
                break;
            case 5:
                saveFunctionMatrixToFile(in_griewangkFilename, data);
                break;
            case 6:
                saveFunctionMatrixToFile(in_sEnvSWaveFilename, data);
                break;
            case 7:
                saveFunctionMatrixToFile(in_strchVSinWaveFilename, data);
                break;
            case 8:
                saveFunctionMatrixToFile(in_ackleys1Filename, data);
                break;
            case 9:
                saveFunctionMatrixToFile(in_ackleys2Filename, data);
                break;
            case 10:
                saveFunctionMatrixToFile(in_eggHolderFilename, data);
                break;
            case 11:
                saveFunctionMatrixToFile(in_ranaFilename, data);
                break;
            case 12:
                saveFunctionMatrixToFile(in_pathologicalFilename, data);
                break;
            case 13:
                saveFunctionMatrixToFile(in_michalewiczFilename, data);
                break;
            case 14:
                saveFunctionMatrixToFile(in_mastersCosWaveFilename, data);
                break;
            case 15:
                saveFunctionMatrixToFile(in_quarticFilename, data);
                break;
            case 16:
                saveFunctionMatrixToFile(in_levyFilename, data);
                break;
            case 17:
                saveFunctionMatrixToFile(in_stepFilename, data);
                break;
            case 18:
                saveFunctionMatrixToFile(in_alpineFilename, data);
                break;
            default:
                std::cout << "Cannot Save Matrix for FunctionData->Function ID: " << funcID << std::endl;
                break;
        }
    }
}

/**
 * @brief Saves all the data in resultsOfFunctions to files.
 */
void ProcessFunctions::saveAllProcessedDataToFile()
{
    for(int numOfResults = 0; numOfResults < resultsOfFunctions.size(); numOfResults++)
    {
        FunctionData data = resultsOfFunctions[numOfResults];
        int funcID = data.functionID;

        // Save with filename referenced by function IDs.
        switch(funcID)
        {
            case 1:
                saveAllFunctionDataToFile(out_schefelsFilename, data);
                break;
            case 2:
                saveAllFunctionDataToFile(out_deJongsFilename, data);
                break;
            case 3:
                saveAllFunctionDataToFile(out_rosenbrockFilename, data);
                break;
            case 4:
                saveAllFunctionDataToFile(out_rastriginFilename, data);
                break;
            case 5:
                saveAllFunctionDataToFile(out_griewangkFilename, data);
                break;
            case 6:
                saveAllFunctionDataToFile(out_sEnvSWaveFilename, data);
                break;
            case 7:
                saveAllFunctionDataToFile(out_strchVSinWaveFilename, data);
                break;
            case 8:
                saveAllFunctionDataToFile(out_ackleys1Filename, data);
                break;
            case 9:
                saveAllFunctionDataToFile(out_ackleys2Filename, data);
                break;
            case 10:
                saveAllFunctionDataToFile(out_eggHolderFilename, data);
                break;
            case 11:
                saveAllFunctionDataToFile(out_ranaFilename, data);
                break;
            case 12:
                saveAllFunctionDataToFile(out_pathologicalFilename, data);
                break;
            case 13:
                saveAllFunctionDataToFile(out_michalewiczFilename, data);
                break;
            case 14:
                saveAllFunctionDataToFile(out_mastersCosWaveFilename, data);
                break;
            case 15:
                saveAllFunctionDataToFile(out_quarticFilename, data);
                break;
            case 16:
                saveAllFunctionDataToFile(out_levyFilename, data);
                break;
            case 17:
                saveAllFunctionDataToFile(out_stepFilename, data);
                break;
            case 18:
                saveAllFunctionDataToFile(out_alpineFilename, data);
                break;
            default:
                std::cout << "Cannot Save Function Results for FunctionData->Function ID: " << funcID << std::endl;
                break;
        }
    }
}

/**
 * @brief Saves all analyzed data in analysis to file.
 */
void ProcessFunctions::saveAllAnalyzedDataToFile()
{
    std::string filename = "zz_AnalyzedData.csv";

    // Rows and Columns.
    int rows = analysis.functionIDs.size(); // Fitness IDs dictates the number of rows.

    // Create the file to where the matrix is saved.
    std::ofstream outputFile;
    outputFile.open (filename);

    // If there are more than 0 fitness IDs, save the header line first.
    if(analysis.functionIDs.size() > 0)
        outputFile << analysis.header;

    // Save data to file.
    std::string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness ID.
        line += std::to_string(analysis.functionIDs[row]) + ",";

        // Save the average fitness.
        line += std::to_string(analysis.avgFuntionFitness[row]) + ",";

        // Save the range.
        line += std::to_string(analysis.ranges[row][0]) + ",";
        line += std::to_string(analysis.ranges[row][1]) + ",";

        // Save the median.
        line += std::to_string(analysis.medianFunctionFitness[row]) + ",";

        // Save the execution time.
        line += std::to_string(analysis.processTimes[row]) + "\n";

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Prints all the possible Function IDs to the screen.
 *
 * Prints all possible Function ID, as well as the funtions they
 * reference, to the screen.
 */
void ProcessFunctions::printAllFunctionIDs()
{
    std::cout << "\n********************************************************\n";
    std::cout << "All Possible Function IDs and Their Respective Functions";
    std::cout << "\n--------------------------------------------------------\n";
    std::cout << "Function ID: 1\tFunction Name: Schwefel’s function\n";
    std::cout << "Function ID: 2\tFunction Name: 1st De Jong’s function\n";
    std::cout << "Function ID: 3\tFunction Name: Rosenbrock function\n";
    std::cout << "Function ID: 4\tFunction Name: Rastrigin function\n";
    std::cout << "Function ID: 5\tFunction Name: Griewangk function\n";
    std::cout << "Function ID: 6\tFunction Name: Sine Envelope Sine Wave function\n";
    std::cout << "Function ID: 7\tFunction Name: Stretched V Sine Wave function\n";
    std::cout << "Function ID: 8\tFunction Name: Ackley’s One function\n";
    std::cout << "Function ID: 9\tFunction Name: Ackley’s Two function\n";
    std::cout << "Function ID: 10\tFunction Name: Egg Holder function\n";
    std::cout << "Function ID: 11\tFunction Name: Rana function\n";
    std::cout << "Function ID: 12\tFunction Name: Pathological function\n";
    std::cout << "Function ID: 13\tFunction Name: Michalewicz function\n";
    std::cout << "Function ID: 14\tFunction Name: Masters Cosine Wave function\n";
    std::cout << "Function ID: 15\tFunction Name: Quartic function\n";
    std::cout << "Function ID: 16\tFunction Name: Levy function\n";
    std::cout << "Function ID: 17\tFunction Name: Step function\n";
    std::cout << "Function ID: 18\tFunction Name: Alpine function\n";
    std::cout << "********************************************************\n\n";
}


// -------------------------------------------------------------------------------------------
// --------------------------------- PRIVATE FUNCTIONS BELOW ---------------------------------
// -------------------------------------------------------------------------------------------
/**
 * @brief Generates a DEFAULT_NUM_OF_VECTORS by DEFAULT_NUM_OF_DIMENSIONS matrix using Mersenne Twister.
 *
 * A matrix is constructed using the default value of 30 dimensions
 * and a user-provided minimum and maximum bound.
 *
 * @note DEFAULT_NUM_OF_VECTORS is currently set to 30 (as of April 4, 2019).
 * @note DEFAULT_NUM_OF_DIMENSIONS is currently set to 30 (as of April 4, 2019).
 *
 * @param minBoundary, maxBoundary The max/min boundaries are the range
 *                                 range in which to generate numbers.
 * @return  The struct that contains the constructed matrix and an empty
 *          list of function fitness results.
 */
ProcessFunctions::FunctionData ProcessFunctions::generateMatrix(double minBoundary, double maxBoundary)
{
    FunctionData generatedData;

    // Create a Mersenne Twister pseudo-random number generator.
    std::mt19937 randGenerator(time(NULL));
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Create the Matrix (30 rows, 30 columns).
    double num, temp;
    for(int row = 0; row < DEFAULT_NUM_OF_VECTORS; row++)
    {
        std::vector<double> matrixRow;
        for (int col = 0; col < numOfDimensions; col++)
        {
            // Generate a random number using Mersenne Twister
            temp = dis(randGenerator);

            // Normalize the random number to the bounds.
            num = (maxBoundary-minBoundary) / (1 - 0) * (temp-0) + minBoundary;

            // Add value to line.
            matrixRow.push_back(num);
        }
        generatedData.functionMatrix.push_back(matrixRow);
    }

    // Return the generated data.
    return generatedData;
}

/**
 * @brief Calculates the fitness of a vector.
 *
 * The fitness of a vector is calculated by the Benchmark Function
 * referenced by the functionID.
 *
 * @param vect The vector of elements on which the Benchmark Functions operate.
 * @param functionID The ID that references which Benchmark Function to use.
 * @return The fitness of the vector.
 */
double ProcessFunctions::calculateFitnessOfVector(std::vector<double> &vect, int functionID)
{
    switch(functionID)
    {
        case 1:
            return schefelsFunc(vect, vect.size());
        case 2:
            return deJongsFunc(vect, vect.size());
        case 3:
            return rosenbrockFunc(vect, vect.size());
        case 4:
            return rastriginFunc(vect, vect.size());
        case 5:
            return griewangkFunc(vect, vect.size());
        case 6:
            return sineEnvelopeSineWaveFunc(vect, vect.size());
        case 7:
            return stretchedVSineWaveFunc(vect, vect.size());
        case 8:
            return ackleysOneFunc(vect, vect.size());
        case 9:
            return ackleysTwoFunc(vect, vect.size());
        case 10:
            return eggHolderFunc(vect, vect.size());
        case 11:
            return ranaFunc(vect, vect.size());
        case 12:
            return pathologicalFunc(vect, vect.size());
        case 13:
            return michalewiczFunc(vect, vect.size());
        case 14:
            return mastersCosWaveFunc(vect, vect.size());
        case 15:
            return quarticFunc(vect, vect.size());
        case 16:
            return levyFunc(vect, vect.size());
        case 17:
            return stepFunc(vect, vect.size());
        case 18:
            return alpineFunc(vect, vect.size());

        default:
            std::cout << "Fitness Process Failed for Function ID: " << functionID << std::endl;
            std::cout << "Possible Function IDs: 1 - 18\n\n";
            return 0.0;
    }
}

/**
 * @brief Calculates the fitness of all vectors of a matrix.
 *
 * Calculates the fitness of all the vectors of the matrix stored
 * in a FunctionData structure. All the fitness results are stored
 * in the fitness vector variable of the same FunctionData structure.
 *
 * @param data The FunctionData structure that contains the matrix.
 */
void ProcessFunctions::calculateFitnessOfMatrix(FunctionData &data)
{
    // Variables to hold function ID and fitness of each vector.
    int funcID = data.functionID;
    double fit;

    // Calculate the fitness of all rows in matrix.
    for(int row = 0; row < data.functionMatrix.size(); row++)
    {
        fit = calculateFitnessOfVector(data.functionMatrix[row], funcID);
        data.fitness.push_back(fit);
    }
}

void ProcessFunctions::analyzeFunctionResults(FunctionData &data) /**< Analyzes the results of the functions.*/
{
    int fitnessSize = data.fitness.size();

    // Save the function ID.
    analysis.functionIDs.push_back(data.functionID);

    // Save the average fitness of tempData.
    double averageFitness = calculateAvgFitness(data);
    analysis.avgFuntionFitness.push_back(averageFitness);

    // Save the fitness ranges.
    std::vector<double> range;
    range.push_back(getMinFitness(data));
    range.push_back(getMaxFitness(data));
    analysis.ranges.push_back(range);

    // Save the median fitness of tempData.
    analysis.medianFunctionFitness.push_back(data.fitness[fitnessSize/2]);

    // Save the execution time of tempData.
    analysis.processTimes.push_back(data.timeToExecute);
}

/**
 * @brief Saves the matrix of the FunctionData to file.
 *
 * @param filename  The filename where to store the matrix. Should be
 *                  a Excel file (.csv).
 * @param data  A FunctionData struct that contains all the data of the function,
 *              including the matrix that was used as well as the fitness
 *              result of that function.
 */
void ProcessFunctions::saveFunctionMatrixToFile(std::string filename, FunctionData &data)
{
    // Rows and Columns of matrix.
    int rows = data.functionMatrix.size();
    int columns = data.functionMatrix[0].size();

    // Create the file to where the matrix is saved.
    std::ofstream outputFile;
    outputFile.open (filename);

    // Save data to file.
    std::string line = "";
    for(int row = 0; row < rows; row++)
    {
        for (int col = 0; col < columns; col++)
        {
            // Add value to line.
            line += std::to_string(data.functionMatrix[row][col]);

            // If the end of the row has been reached.
            if(col == columns-1)
                line += "\n"; // Add a newline.
            else
                line += ","; // Else, add a comma.
        }

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}

/**
 * @brief Saves all the data of the function to file.
 *
 * @param filename  The filename where to store the matrix. Should be
 *                  a Excel file (.csv).
 * @param data  A FunctionData struct that contains all the data of the function,
 *              including the matrix that was used as well as the fitness
 *              result of that function.
 */
void ProcessFunctions::saveAllFunctionDataToFile(std::string filename, FunctionData &data)
{
    // Rows and Columns.
    int rows = data.fitness.size();                 // Fitness dictates the number of rows.
    int columns = data.functionMatrix[0].size();    // Matrix dictates the number of dimensions.

    // Create the file to where the matrix is saved.
    std::ofstream outputFile;
    outputFile.open (filename);

    // Save data to file.
    std::string line = "";
    for(int row = 0; row < rows; row++)
    {
        // Save the fitness.
        line += std::to_string(data.fitness[row]) + ",";

        for (int col = 0; col < columns; col++)
        {
            // Add value to line.
            line += std::to_string(data.functionMatrix[row][col]);

            // If the end of the row has been reached.
            if(col == columns-1)
                line += "\n"; // Add a newline.
            else
                line += ","; // Else, add a comma.
        }

        // Save the row to file and clear the line string.
        outputFile << line;
        line = "";
    }

    // Close the file.
    outputFile.close();
}


// *****************************************************************
void ProcessFunctions::quicksort(FunctionData &data, int L, int R) {
    int i, j, mid;
    double piv;
    i = L;
    j = R;
    mid = L + (R - L) / 2;
    piv = data.fitness[mid];

    while (i<R || j>L)
    {
        while (data.fitness[i] < piv)
            i++;
        while (data.fitness[j] > piv)
            j--;

        if (i <= j)
        {
            swap(data, i, j); //error=swap function doesnt take 3 arguments
            i++;
            j--;
        }
        else
        {
            if (i < R)
                quicksort(data, i, R);
            if (j > L)
                quicksort(data, L, j);
            return;
        }
    }
}

void ProcessFunctions::swap(FunctionData &data, int x, int y)
{
    // Swap fitness values.
    double fitTemp = data.fitness[x];
    data.fitness[x] = data.fitness[y];
    data.fitness[y] = fitTemp;

    // Swap vector values.
    std::vector<double> vectTemp = data.functionMatrix[x];
    data.functionMatrix[x] = data.functionMatrix[y];
    data.functionMatrix[y] = vectTemp;
}