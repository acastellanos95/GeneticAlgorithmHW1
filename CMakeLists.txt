cmake_minimum_required(VERSION 3.22)
project(GeneticAlgorithm)

set(CMAKE_CXX_STANDARD 17)

add_executable(GeneticAlgorithm main.cpp Individual.h lib/VariadicTable.h lib/RandomGA.h Utils.h)
