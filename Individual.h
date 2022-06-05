//
// Created by andre on 5/28/22.
//

#ifndef GENETICALGORITHM__INDIVIDUAL_H_
#define GENETICALGORITHM__INDIVIDUAL_H_

#include <vector>

/**
 * Estructura para representar un individuo
 */
struct individual {
  std::string chromosome; /* Cromosoma */
  std::vector<double> x; /* Cromosoma decodificado */
  double fitness; /* Valor de aptitud */
  size_t crossoverPoint; /* Punto de cruza */
  bool hasCrossover; /* Punto de cruza */
  size_t numberMutations; /* Número de mutaciones en el cromosoma */
  std::vector<unsigned long> parents = std::vector<unsigned long>(2, 0); /* Padres de individuo */
};

#endif //GENETICALGORITHM__INDIVIDUAL_H_
