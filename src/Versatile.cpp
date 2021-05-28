#include "Versatile.h"

#include <algorithm>  // std::max, std::min
#include <cmath>  // std::pow

#include "Matrix.h"
#include "Sahn.h"

mdendro::Versatile::Versatile() :
  Sahn() {
  this->power = +1.0;
}

mdendro::Versatile::Versatile(double power, bool isWeighted,
    const Matrix& proximity, bool isDistance, int precision, bool isVariable) :
  Sahn(isWeighted, proximity, isDistance, precision, isVariable) {
  this->power = power;
}

mdendro::Versatile::~Versatile() {}

double mdendro::Versatile::minimumProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  double prox = +INF;
  std::list<int>::const_iterator iti = nni.begin();
  while (iti != nni.end()) {
    int i = *iti;
    std::list<int>::const_iterator itj = nnj.begin();
    while (itj != nnj.end()) {
      int j = *itj;
      double proxij = this->proximity.getValue(i, j);
      prox = std::min(prox, proxij);
      itj ++;
    }
    iti ++;
  }
  return prox;
}

double mdendro::Versatile::maximumProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  double prox = -INF;
  std::list<int>::const_iterator iti = nni.begin();
  while (iti != nni.end()) {
    int i = *iti;
    std::list<int>::const_iterator itj = nnj.begin();
    while (itj != nnj.end()) {
      int j = *itj;
      double proxij = this->proximity.getValue(i, j);
      prox = std::max(prox, proxij);
      itj ++;
    }
    iti ++;
  }
  return prox;
}

double mdendro::Versatile::geometricMean(const std::list<int>& nni,
    const std::list<int>& nnj) {
  std::pair<int, int> smi = sumMembers(nni);
  std::pair<int, int> smj = sumMembers(nnj);
  int smi1 = smi.first;
  int smj1 = smj.first;
  double prox = 1.0;
  std::list<int>::const_iterator iti = nni.begin();
  while (iti != nni.end()) {
    int i = *iti;
    int mi = (this->isWeighted)? 1 : this->clusters[i].nMembers;
    std::list<int>::const_iterator itj = nnj.begin();
    while (itj != nnj.end()) {
      int j = *itj;
      int mj = (this->isWeighted)? 1 : this->clusters[j].nMembers;
      double proxij = this->proximity.getValue(i, j);
      prox *= std::pow(proxij, (double)(mi * mj) / (double)(smi1 * smj1));
      itj ++;
    }
    iti ++;
  }
  return prox;
}

double mdendro::Versatile::generalizedMean(const std::list<int>& nni,
    const std::list<int>& nnj) {
  std::pair<int, int> smi = sumMembers(nni);
  std::pair<int, int> smj = sumMembers(nnj);
  int smi1 = smi.first;
  int smj1 = smj.first;
  double prox = 0.0;
  std::list<int>::const_iterator iti = nni.begin();
  while (iti != nni.end()) {
    int i = *iti;
    int mi = (this->isWeighted)? 1 : this->clusters[i].nMembers;
    std::list<int>::const_iterator itj = nnj.begin();
    while (itj != nnj.end()) {
      int j = *itj;
      int mj = (this->isWeighted)? 1 : this->clusters[j].nMembers;
      double proxij = this->proximity.getValue(i, j);
      prox += ((double)(mi * mj) / (double)(smi1 * smj1))
          * std::pow(proxij, this->power);
      itj ++;
    }
    iti ++;
  }
  prox = std::pow(prox, 1.0 / this->power);
  return prox;
}

double mdendro::Versatile::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  double prox;
  if (this->power == -INF) {
    prox = minimumProximity(nni, nnj);
  } else if (this->power == +INF) {
    prox = maximumProximity(nni, nnj);
  } else if (this->power == 0.0) {
    prox = geometricMean(nni, nnj);
  } else {
    prox = generalizedMean(nni, nnj);
  }
  return prox;
}
