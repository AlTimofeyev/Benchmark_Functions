/**
 * @file FilenameConstants.h
 * @author  Al Timofeyev
 * @date    March 28, 2019
 * @brief   A list of input and output filenames. The input files
 *          are where the matrices are stored. The output files
 *          are where the results from the benchmark functions
 *          are stored.
 */

#ifndef BENCHMARKFUNCTIONS_FILENAMECONSTANTS_H
#define BENCHMARKFUNCTIONS_FILENAMECONSTANTS_H


#include <string>

// ---------- INPUT FILES ----------
std::string in_schefelsFilename = "z_schefels-Data.csv";
std::string in_deJongsFilename = "z_deJongs-Data.csv";
std::string in_rosenbrockFilename = "z_rosenbrock-Data.csv";
std::string in_rastriginFilename = "z_rastrigin-Data.csv";
std::string in_griewangkFilename = "z_griewangk-Data.csv";
std::string in_sEnvSWaveFilename = "z_sEnvSWave-Data.csv";
std::string in_strchVSinWaveFilename = "z_strchVSinWave-Data.csv";
std::string in_ackleys1Filename = "z_ackleys1-Data.csv";
std::string in_ackleys2Filename = "z_ackleys2-Data.csv";
std::string in_eggHolderFilename = "z_eggHolder-Data.csv";
std::string in_ranaFilename = "z_rana-Data.csv";
std::string in_pathologicalFilename = "z_pathological-Data.csv";
std::string in_michalewiczFilename = "z_michalewicz-Data.csv";
std::string in_mastersCosWaveFilename = "z_mastersCosWave-Data.csv";
std::string in_quarticFilename = "z_quartic-Data.csv";
std::string in_levyFilename = "z_levy-Data.csv";
std::string in_stepFilename = "z_step-Data.csv";
std::string in_alpineFilename = "z_alpine-Data.csv";

// ---------- OUTPUT FILES ----------
std::string out_schefelsFilename = "z_schefels-Output.csv";
std::string out_deJongsFilename = "z_deJongs-Output.csv";
std::string out_rosenbrockFilename = "z_rosenbrock-Output.csv";
std::string out_rastriginFilename = "z_rastrigin-Output.csv";
std::string out_griewangkFilename = "z_griewangk-Output.csv";
std::string out_sEnvSWaveFilename = "z_sEnvSWave-Output.csv";
std::string out_strchVSinWaveFilename = "z_strchVSinWave-Output.csv";
std::string out_ackleys1Filename = "z_ackleys1-Output.csv";
std::string out_ackleys2Filename = "z_ackleys2-Output.csv";
std::string out_eggHolderFilename = "z_eggHolder-Output.csv";
std::string out_ranaFilename = "z_rana-Output.csv";
std::string out_pathologicalFilename = "z_pathological-Output.csv";
std::string out_michalewiczFilename = "z_michalewicz-Output.csv";
std::string out_mastersCosWaveFilename = "z_mastersCosWave-Output.csv";
std::string out_quarticFilename = "z_quartic-Output.csv";
std::string out_levyFilename = "z_levy-Output.csv";
std::string out_stepFilename = "z_step-Output.csv";
std::string out_alpineFilename = "z_alpine-Output.csv";

// ---------- FITNESS RESULTS FILES ----------
std::string out_Fitness10Dimensions = "10DimensionFitness-ResultsAnalysis";
std::string out_Fitness20Dimensions = "20DimensionFitness-ResultsAnalysis";
std::string out_Fitness30Dimensions = "30DimensionFitness-ResultsAnalysis";


#endif //BENCHMARKFUNCTIONS_FILENAMECONSTANTS_H