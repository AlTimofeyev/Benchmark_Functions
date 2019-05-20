******************************************************************************************************************
Author:	Al Timofeyev
Date:	April 5, 2019
Desc:	This is how to compile and execute the code.
******************************************************************************************************************

Environment Used to Code:
Windows 10
CLion version 2019.1
cygwin version 3.0.4
cygwin GDB version 8.1.1
gcc version 7.4.0 
g++ version 7.4.0

****** NOTE:
CLion generated a CMakeLists.txt file included with the source code.
cmake_minimum_required(VERSION 3.13)


To Execute:
run main.cpp

To Compile:
You could use CMake to compile CMakeLists.txt OR
Compile main.cpp using g++

****** NOTE:
1) main.cpp is just a test driver for the ProcessFunctions class.
   Feel free to change and add to it to generate resutls.
2) Every time main.cpp is executed, the previous files that held
   results will be overwritten.
3) BenchmarkFunctions.h/cpp are stand-alone benchmark testing functions
   that can be imported and used in other projects/files.
4) ProcessFunctions.h/cpp is used to test the Benchmark Functions.
5) ALL source files need to be included for this project (main.cpp) to work.