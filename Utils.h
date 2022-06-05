//
// Created by andre on 5/28/22.
//

#ifndef GENETICALGORITHM__UTILS_H_
#define GENETICALGORITHM__UTILS_H_

#include <vector>
#include <bitset>
#include <cmath>
#include "lib/VariadicTable.h"

/* Definiciones de problema */
const size_t numVariables = 4;
const size_t numBitsVariable = 39; /* log2([20.0-(-20.0)x10^10]) = 38.54 => floor(38.54) + 1 = 39 */
const double l_sup = 20.0;
const double l_inf = -20.0;

/// Función a minimizar f(x1,x2,x3,x4)=(10(x2-x1^2))^2 + (1-x1)^2 + 90(x4-x3^2)^2 + (1-x3)^2 + 10(x2 + x4-2)^2 + 0.1(x2-x4)^2
/// \param x
/// \return f(x)
double function(const std::vector<double> &x){
  double result = (10.0*(x[1]-x[0]*x[0]))*(10.0*(x[1]-x[0]*x[0])) + (1.0-x[0])*(1.0-x[0]) + 90.0*(x[3]-x[2]*x[2])*(x[3]-x[2]*x[2]) + (1.0-x[2])*(1.0-x[2]) + 10.0*(x[1]+x[3]-2.0)*(x[1]+x[3]-2.0) + 0.1*(x[1]-x[3])*(x[1]-x[3]);
  return result;
}

/// Decodificar cromosoma a x O(numVariables)
/// \param cromosoma: Cadena binaria que representa
/// \return Vector x
std::vector<double> decode(const std::string &cromosoma){
  // Resultado
  std::vector<double> result(numVariables);
  // Vector de bitset de cada variable x_i
  std::vector<std::bitset<numBitsVariable>> variableGenes(numVariables);

  // Copiar de cadena al vector de bitset
  for(size_t i = 0; i < numVariables; ++i){
    std::string geneString = cromosoma.substr(i*numBitsVariable,numBitsVariable);
    variableGenes[i] = std::bitset<numBitsVariable>(geneString);
  }

  // Decodificación
  for(size_t i = 0; i < numVariables; ++i){
    double x_i = l_inf + ( (double ) variableGenes[i].to_ullong()) * ((l_sup - l_inf) / (pow(2.0, numBitsVariable) - 1.0));
    result[i] = x_i;
  }

  return result;
}

/// Calcular el promedio del vector de valores x
/// \param x
/// \return avg(x)
double avg(const std::vector<double> &x){
  double sum = 0.0;
  for(auto &x_i: x){
    sum += x_i;
  }
  return sum/(double ) x.size();
}

/// TODO: reporte en modo debug
/// \param selected: arreglo de par de índices individuos en la población seleccionados para cruza
/// \param oldPopulation: población de padres
/// \param newPopulation: población de hijos
/// \param pc: probabilidad de cruza
void onePointCrossover(const std::vector<std::pair<size_t,size_t>> &selected, const std::vector<individual> &oldPopulation, std::vector<individual> &newPopulation, const double pc){
  size_t numberCrossover = 0;
  size_t newPopIndex = 0;
  for(auto &selectPair: selected){
    if(flip(pc)){
      ++numberCrossover;
      /* Punto de cruza aleatorio */
      size_t crosspoint = rnd(0, numVariables*numBitsVariable);

      newPopulation[newPopIndex].chromosome = oldPopulation[selectPair.first].chromosome.substr(0,crosspoint);
      newPopulation[newPopIndex + 1].chromosome = oldPopulation[selectPair.second].chromosome.substr(0,crosspoint);

      newPopulation[newPopIndex].chromosome += oldPopulation[selectPair.second].chromosome.substr(crosspoint,numVariables*numBitsVariable);
      newPopulation[newPopIndex + 1].chromosome += oldPopulation[selectPair.first].chromosome.substr(crosspoint,numVariables*numBitsVariable);

      newPopulation[newPopIndex].parents = {selectPair.first, selectPair.second};
      newPopulation[newPopIndex + 1].parents = {selectPair.first, selectPair.second};

      newPopulation[newPopIndex].crossoverPoint = crosspoint;
      newPopulation[newPopIndex + 1].crossoverPoint = crosspoint;

      newPopulation[newPopIndex].hasCrossover = true;
      newPopulation[newPopIndex + 1].hasCrossover = true;
    } else {
      newPopulation[newPopIndex].chromosome = oldPopulation[selectPair.first].chromosome;
      newPopulation[newPopIndex + 1].chromosome = oldPopulation[selectPair.second].chromosome;

      newPopulation[newPopIndex].parents = {selectPair.first, selectPair.second};
      newPopulation[newPopIndex + 1].parents = {selectPair.first, selectPair.second};

      newPopulation[newPopIndex].hasCrossover = false;
      newPopulation[newPopIndex + 1].hasCrossover = false;
    }
    newPopIndex += 2;
  }
  std::cout << "Número de cruzas: " << numberCrossover << '\n';
}

/// Muta los alelos de cada cromosoma de un individuo
/// \param newPopulation: población de hijos a mutar
/// \param pm: probabilidad de mutación
/// \return lista de mutaciones
void mutations(std::vector<individual> &newPopulation, const double pm){
  size_t numberMutations = 0;
  for(auto &individuo: newPopulation){
    size_t individualNumberMutations = 0;
    for(auto &allele: individuo.chromosome){
      if(flip(pm)){
        ++numberMutations;
        ++individualNumberMutations;
        allele = allele == '0'? '1' : '0';
      }
    }
    individuo.numberMutations = individualNumberMutations;
  }

  std::cout << "Número de mutaciones totales: " << numberMutations << '\n';
}

/// Tomar el más apto y ponerlo en el primer hijo
/// \param newPopulation: población de hijos
/// \param bestIndividual: individuo más apto de la población de padres
void elitism(std::vector<individual> &newPopulation, individual &bestIndividual){
  newPopulation[0] = bestIndividual;
}


void reportGeneration(const std::vector<individual> &oldPopulation, const std::vector<individual> &newPopulation, const std::vector<double> &fitness, const std::vector<individual> &newPopulationBeforeMutation){
  // Format: popIndex, chromosome, x value, fitness, puntos de cruza, padres, hijos producidos, numero de mutaciones, efectos de mutación si hay
  VariadicTable<size_t, std::string, std::string, std::string, std::string, std::string, std::string, size_t , std::string> vt({"# individuo", "Cromosoma", "x", "Aptitud", "Punto de cruza", "Padres", "Hijos antes de mutar", "# mutaciones", "Hijos después de mutar"});

  for(size_t popIndex = 0; popIndex < oldPopulation.size(); ++popIndex){
    // x value
    auto x = decode(oldPopulation[popIndex].chromosome);

    std::stringstream xss;
    xss << std::setprecision(5) << std::scientific;
    xss << "(" << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << ")";

    // fitness
    std::stringstream fitnessSS;
    fitnessSS << std::setprecision(5) << std::scientific;
    fitnessSS << oldPopulation[popIndex].fitness;

    // padres
    std::stringstream parentsSS;
    parentsSS << "(" << newPopulation[popIndex].parents[0] << "," << newPopulation[popIndex].parents[1] << ")";

    vt.addRow(popIndex, oldPopulation[popIndex].chromosome, xss.str(), fitnessSS.str(), newPopulation[popIndex].hasCrossover ? std::to_string(newPopulation[popIndex].crossoverPoint): "---", parentsSS.str(), newPopulationBeforeMutation[popIndex].chromosome, newPopulation[popIndex].numberMutations, newPopulation[popIndex].chromosome);
  }

  vt.print(std::cout);
}

#endif //GENETICALGORITHM__UTILS_H_
