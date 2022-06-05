#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include "Individual.h"
#include "lib/RandomGA.h"
#include "Utils.h"

int main(int argc, char *argv[]) {
  double pm, pc;
  size_t Gmax, populationSize;
  std::string debug;
  if(argc == 1){
    std::cout << "Introduzca la probabilidad de mutación (0.01 <= pm <= 0.1): ";
    std::cin >> pm;
    std::cout << "Introduzca la probabilidad de cruza (0.6 <= pm <= 1.0): ";
    std::cin >> pc;
    std::cout << "Introduzca el número máximo de generaciones (>=100): ";
    std::cin >> Gmax;
    std::cout << "Introduzca el tamaño máximo de población (>=50): ";
    std::cin >> populationSize;
    std::cout << "Semilla inicial (0.0 <= seed <= 1.0): ";
    std::cin >> Rseed;
    std::cout << "Debug(Y/N): ";
    std::cin >> debug;
  } else if(argc == 7) {
    pm = std::strtod(argv[1], nullptr);
    pc = std::strtod(argv[2], nullptr);
    Gmax = std::strtoul(argv[3], nullptr,0);
    populationSize = std::strtoul(argv[4], nullptr,0);
    Rseed = std::strtof(argv[5], nullptr);
    debug = std::string(argv[6]);
  } else {
    throw std::runtime_error("Formato esperado: probabilidad de mutación (0.01 <= pm <= 0.1), probabilidad de cruza (0.6 <= pm <= 1.0), número máximo de generaciones (>=100), tamaño máximo de población (>=50), Semilla inicial (0.0 <= seed <= 1.0), Modo debug (Y/N)");
  }

  // Verificación de condiciones (petición de clase)
  if(pm > 0.1)
    std::cout << "Advertencia!!!! probabilidad de mutación mayor que la sugerida (0.01 <= pm <= 0.1)\n";
  if(pm > 1.0 || pm < 0.0)
    throw std::runtime_error("Probabilidad de mutación no puede ser < 0.0 ó > 1.0");
  if(pc < 0.6)
    std::cout << "Advertencia!!!! probabilidad de cruza menor que la sugerida (0.6 <= pm <= 1.0)\n";
  if(pm > 1.0 || pm < 0.0)
    throw std::runtime_error("Probabilidad de cruza no puede ser < 0.0 ó > 1.0");
  if(Gmax < 100)
    std::cout << "Advertencia!!!! número de generaciones menor que la sugerida (>=100)\n";
  if(populationSize < 50)
    std::cout << "Advertencia!!!! tamaño de población menor que la sugerida (>=50)\n";
  if(Rseed < 0.0 || Rseed > 1.0)
    throw std::runtime_error("Semilla no puede ser < 0.0 ó > 1.0");

  // Inicializa población y generador de números aleatorios
  randomize();
  std::vector<individual> oldPopulation(populationSize);
  for(auto &individuo: oldPopulation){
    // Inicializar cromosoma
    std::string chrom(numVariables*numBitsVariable, '0');
    for(size_t alleleIndex = 0; alleleIndex < numVariables*numBitsVariable; ++alleleIndex){
      // fitness 0.1*maxfx
      // 0.05 0.95 0.2154 5000 it, 500 pob (1.09497e+00,1.18746e+00,9.76564e-01,9.39749e-01)
      if(flip(0.1)){
        chrom[alleleIndex]= '1';
      } else {
        chrom[alleleIndex]= '0';
      }
    }
    individuo.chromosome = chrom;

    // Decodificar
    individuo.x = decode(individuo.chromosome);
  }
  std::vector<individual> newPopulation(populationSize);

  // Reporte Inicial

  std::cout.precision(5);
  std::cout <<std::scientific;
  std::cout << "probabilidad de mutación: " << pm << ", probabilidad de cruza: " << pc << ", número de generaciones: " << Gmax << ", tamaño de población: " << populationSize << '\n';

  // Generaciones y variables globales para generación
  std::pair<size_t, individual> bestIndividual;
  // Solo para punto 5
  std::vector<std::pair<size_t, double>> convergence;
  double maxFx = -std::numeric_limits<double>::max();
  for(size_t generationIndex = 1; generationIndex <= Gmax; ++generationIndex){
    std::cout << "--------------------- Generación " << generationIndex << " ---------------------\n";
    // Selección proporcional de ruleta (maximización) TODO: minimización
    std::vector<double> fitness(populationSize, 0.0);
    std::vector<double> expected(populationSize, 0.0);
    individual generationBestIndividual;

    // Valores de f(x) a fitness y buscar el máximo valor
    for(size_t popIndex = 0; popIndex < populationSize; ++popIndex){
      fitness[popIndex] = function(oldPopulation[popIndex].x);
    }

    auto maxIndividualFx = std::max_element(fitness.begin(), fitness.end());
    maxFx = *maxIndividualFx > maxFx ? *maxIndividualFx: maxFx;

    for(size_t popIndex = 0; popIndex < populationSize; ++popIndex){
      fitness[popIndex] = maxFx - fitness[popIndex] + 0.1*maxFx;
      oldPopulation[popIndex].fitness = fitness[popIndex];
    }

    double avgfitness = avg(fitness);
    size_t maxElementIndex = std::max_element(fitness.begin(), fitness.end()) - fitness.begin();
    double maxfitness = *std::max_element(fitness.begin(), fitness.end());
    double minfitness = *std::min_element(fitness.begin(), fitness.end());

    // Mejor individuo de la generación
    generationBestIndividual = oldPopulation[maxElementIndex];
    if(bestIndividual.second.chromosome.empty()){
      bestIndividual.second = generationBestIndividual;
      bestIndividual.first = generationIndex;
    } else {
      // Actualizar mejor individuo global para comparar en caso de que maxFx haya cambiado
      bestIndividual.second.fitness = maxFx - function(bestIndividual.second.x) + 0.1*maxFx;
      if(bestIndividual.second.fitness < generationBestIndividual.fitness)
        bestIndividual.second = generationBestIndividual;
      bestIndividual.first = generationIndex;
    }

    std::cout << "Media de aptitud de población: " << avgfitness << '\n';
    std::cout << "Aptitud máxima de población: " << maxfitness << '\n';
    std::cout << "Aptitud mínima de población: " << minfitness << '\n';
    std::cout << "Mejor individuo de la generación: " << generationBestIndividual.chromosome << ", (" << generationBestIndividual.x[0] << "," << generationBestIndividual.x[1] << "," << generationBestIndividual.x[2] << "," << generationBestIndividual.x[3] << "), fitness: " << generationBestIndividual.fitness << '\n';
    std::cout << "Mejor individuo global: " << bestIndividual.second.chromosome << ", (" << bestIndividual.second.x[0] << "," << bestIndividual.second.x[1] << "," << bestIndividual.second.x[2] << "," << bestIndividual.second.x[3] << "), fitness: " << bestIndividual.second.fitness << ", generación: " << bestIndividual.first << '\n';
    // Punto 5 gráfica de convergencia
    convergence.emplace_back(generationIndex, generationBestIndividual.fitness);

    double sumExpected = 0.0;
    for(size_t popIndex = 0; popIndex < populationSize; ++popIndex){
      expected[popIndex] = fitness[popIndex]/avgfitness;
      sumExpected += expected[popIndex];
    }

    // Imprimir tabla de aptitudes y valores esperados de la generación (punto 4)
    if(debug == "Y"){
      // Format: popIndex, fitness, expected value
      VariadicTable<size_t, std::string, std::string> vtExpectedTable({"# individuo", "Aptitud", "Valor esperado"});
      double sumfitness = 0.0;
      for(size_t popIndex = 0; popIndex < populationSize; ++popIndex){
        std::stringstream fitnessSS;
        fitnessSS << std::setprecision(5) << std::scientific;
        fitnessSS << fitness[popIndex];
        std::stringstream expectedSS;
        expectedSS << std::setprecision(5) << std::scientific;
        expectedSS << expected[popIndex];
        vtExpectedTable.addRow(popIndex, fitnessSS.str(), expectedSS.str());
        sumfitness += fitness[popIndex];
      }
      vtExpectedTable.print(std::cout);
      std::cout << "Suma de todas las aptitudes: " << sumfitness << ", promedio de las aptitudes: " << avgfitness << '\n';
      std::cout << "Suma de los valores esperados: " << sumExpected << "\n";
    }

    // Seleccionamos pares para cruza
    std::vector<std::pair<size_t,size_t>> selected;
    std::vector<size_t> freqSelected(populationSize, 0);
    for(size_t popIndex = 0; popIndex < populationSize; popIndex += 2){

      // Ilustración de selección (punto 4)
      // Format: popIndex, fitness, expected value
      VariadicTable<size_t, std::string, std::string, std::string> vtSelect1({"# individuo", "Aptitud", "Valor esperado", "Suma"});
      VariadicTable<size_t, std::string, std::string, std::string> vtSelect2({"# individuo", "Aptitud", "Valor esperado", "Suma"});

      double r1 = rndreal(0.0, sumExpected);
      size_t ind1 = 0;
      double sum1 = expected[ind1];
      if(debug == "Y"){
        std::cout << "r1: " << r1 << '\n';
        std::stringstream fitnessSS;
        fitnessSS << std::setprecision(5) << std::scientific;
        fitnessSS << fitness[ind1];
        std::stringstream expectedSS;
        expectedSS << std::setprecision(5) << std::scientific;
        expectedSS << expected[ind1];
        std::stringstream sumSS;
        sumSS << std::setprecision(5) << std::scientific;
        sumSS << sum1;

        vtSelect1.addRow(ind1, fitnessSS.str(), expectedSS.str(), sumSS.str());
      }

      while (sum1 < r1 || ind1 == populationSize){
        sum1 += expected[++ind1];

        if(debug == "Y"){
          std::stringstream fitnessSS;
          fitnessSS << std::setprecision(5) << std::scientific;
          fitnessSS << fitness[ind1];
          std::stringstream expectedSS;
          expectedSS << std::setprecision(5) << std::scientific;
          expectedSS << expected[ind1];
          std::stringstream sumSS;
          sumSS << std::setprecision(5) << std::scientific;
          sumSS << sum1;

          vtSelect1.addRow(ind1, fitnessSS.str(), expectedSS.str(), sumSS.str());
        }
      }

      double r2 = rndreal(0.0, sumExpected);
      size_t ind2 = 0;
      double sum2 = expected[ind2];
      if(debug == "Y"){
        std::cout << "r2: " << r2 << '\n';
        std::stringstream fitnessSS;
        fitnessSS << std::setprecision(5) << std::scientific;
        fitnessSS << fitness[ind2];
        std::stringstream expectedSS;
        expectedSS << std::setprecision(5) << std::scientific;
        expectedSS << expected[ind2];
        std::stringstream sumSS;
        sumSS << std::setprecision(5) << std::scientific;
        sumSS << sum2;

        vtSelect2.addRow(ind2, fitnessSS.str(), expectedSS.str(), sumSS.str());
      }

      while (sum2 < r2  || ind2 == populationSize){
        sum2 += expected[++ind2];

        if(debug == "Y"){
          std::stringstream fitnessSS;
          fitnessSS << std::setprecision(5) << std::scientific;
          fitnessSS << fitness[ind2];
          std::stringstream expectedSS;
          expectedSS << std::setprecision(5) << std::scientific;
          expectedSS << expected[ind2];
          std::stringstream sumSS;
          sumSS << std::setprecision(5) << std::scientific;
          sumSS << sum2;

          vtSelect2.addRow(ind2, fitnessSS.str(), expectedSS.str(), sumSS.str());
        }
      }

      selected.emplace_back(ind1,ind2);
      ++freqSelected[ind1];
      ++freqSelected[ind2];

      // Imprimir selección de parejas (punto 4)
      if(debug == "Y"){
        vtSelect1.print(std::cout);
        vtSelect2.print(std::cout);
        std::cout << "Seleccionamos " << ind1 << " y " << ind2 << "\n";
      }
    }
    // Imprimir frecuencias de selección
    if(debug == "Y"){
      // Format: popIndex, fitness, expected value
      VariadicTable<size_t, size_t> vtSelectionFreq({"# individuo", "número de copias"});
      for(size_t popIndex = 0; popIndex < populationSize; ++popIndex){
        vtSelectionFreq.addRow(popIndex, freqSelected[popIndex]);
      }
      vtSelectionFreq.print(std::cout);
    }

    // Cruza de un punto
    onePointCrossover(selected, oldPopulation, newPopulation, pc);

    std::vector<individual> newPopulationBeforeMutation = newPopulation;

    // Mutación
    mutations(newPopulation, pm);

    // Imprimir reporte de padres e hijos, puntos de cruza, mutaciones (punto 4)
    if(debug == "Y")
      reportGeneration(oldPopulation, newPopulation, fitness, newPopulationBeforeMutation);

    // Elitismo
    elitism(newPopulation, generationBestIndividual);

    // Actualizar x en la nueva población
    for(auto &individuo: newPopulation){
      // Decodificar
      individuo.x = decode(individuo.chromosome);
    }

    auto tmpPopulation = oldPopulation;
    oldPopulation = newPopulation;
    newPopulation = tmpPopulation;
  }

  // Punto 5 estadísticas
//  std::cout << "probabilidad de mutación: " << pm << ", probabilidad de cruza: " << pc << ", número de generaciones: " << Gmax << ", tamaño de población: " << populationSize << '\n';
//  std::cout << "Final f(x) " << function(bestIndividual.second.x) << '\n';
//  std::cout << "$(" << bestIndividual.second.x[0] << "," << bestIndividual.second.x[1] << "," << bestIndividual.second.x[2] << "," << bestIndividual.second.x[3] << ")$ & $" << function(bestIndividual.second.x) << "$ & $" << bestIndividual.second.fitness << "$\\\\\n";
  // Punto 5 convergencia de la mediana
//  std::ofstream file("convergence.dat");
//  for(auto &pair: convergence){
//    file << pair.first << ", " << std::to_string(pair.second) << '\n';
//  }
//  file.close();
  return 0;
}
