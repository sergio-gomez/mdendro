#include "Centroid.h"

#include <cmath>  // std::sqrt

#include "LanceWilliams.h"
#include "Matrix.h"

mdendro::Centroid::Centroid() :
  LanceWilliams() {}

mdendro::Centroid::Centroid(bool isWeighted, const Matrix& proximity,
    bool isDistance, int precision, bool isVariable) :
  LanceWilliams(isWeighted, proximity, isDistance, precision, isVariable) {
}

mdendro::Centroid::~Centroid() {}

double mdendro::Centroid::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  std::pair<int, int> smi = sumMembers(nni);
  std::pair<int, int> smj = sumMembers(nnj);
  double alphaij = alphaTerm(nni, nnj, smi, smj);
  double betai = betaTerm(nni, smi, smj);
  double betaj = betaTerm(nnj, smj, smi);
  // Return square root
  return std::sqrt(alphaij + betai + betaj);
}

double mdendro::Centroid::getAlphaProximity(int mi, int mj,
    std::pair<int, int> smi, std::pair<int, int> smj, double pij) {
  int smi1 = smi.first;
  int smj1 = smj.first;
  // Use squared distance
  return ((double)(mi * mj) / (double)(smi1 * smj1)) * pij * pij;
}

double mdendro::Centroid::getBetaProximity(int mi1, int mi2,
    std::pair<int, int> smi, std::pair<int, int> smj, double pi) {
  int smi1 = smi.first;
  // Use squared distance
  return -((double)(mi1 * mi2) / (double)(smi1 * smi1)) * pi * pi;
}
