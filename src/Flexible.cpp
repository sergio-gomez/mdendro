#include "Flexible.h"

#include <algorithm>  // std::max, std::min

#include "LanceWilliams.h"
#include "Matrix.h"

mdendro::Flexible::Flexible() :
  LanceWilliams() {
  this->beta = 0.0;
}

mdendro::Flexible::Flexible(double beta, bool isWeighted,
    const Matrix& proximity, bool isDistance, int precision, bool isVariable) :
  LanceWilliams(isWeighted, proximity, isDistance, precision, isVariable) {
  // -1.0 <= beta <= +1.0
  this->beta = beta;
  this->beta = std::max(this->beta, -1.0);
  this->beta = std::min(this->beta, +1.0);
}

mdendro::Flexible::~Flexible() {}

double mdendro::Flexible::getAlphaProximity(int mi, int mj,
    std::pair<int, int> smi, std::pair<int, int> smj, double pij) {
  int smi1 = smi.first;
  int smj1 = smj.first;
  return (1.0 - this->beta) * ((double)(mi * mj) / (double)(smi1 * smj1)) * pij;
}

double mdendro::Flexible::getBetaProximity(int mi1, int mi2,
    std::pair<int, int> smi, std::pair<int, int> smj, double pi) {
  int smi1 = smi.first;
  int smj1 = smj.first;
  int smi2 = smi.second;
  int smj2 = smj.second;
  int sigmai = (smi1 * smi1 - smi2) / 2;
  int sigmaj = (smj1 * smj1 - smj2) / 2;
  return this->beta * ((double)(mi1 * mi2) / (double)(sigmai + sigmaj)) * pi;
}
