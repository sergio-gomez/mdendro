#include "LanceWilliams.h"

#include <iterator>  // std::next

#include "Matrix.h"
#include "Sahn.h"

mdendro::LanceWilliams::LanceWilliams() :
  Sahn() {}

mdendro::LanceWilliams::LanceWilliams(bool isWeighted, const Matrix& proximity,
    bool isDistance, int precision, bool isVariable) :
  Sahn(isWeighted, proximity, isDistance, precision, isVariable) {
}

mdendro::LanceWilliams::~LanceWilliams() {}

double mdendro::LanceWilliams::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  std::pair<int, int> smi = sumMembers(nni);
  std::pair<int, int> smj = sumMembers(nnj);
  double alphaij = alphaTerm(nni, nnj, smi, smj);
  double betai = betaTerm(nni, smi, smj);
  double betaj = betaTerm(nnj, smj, smi);
  return alphaij + betai + betaj;
}

double mdendro::LanceWilliams::alphaTerm(const std::list<int>& nni,
    const std::list<int>& nnj, std::pair<int, int> smi,
    std::pair<int, int> smj) {
  double alphaij = 0.0;
  std::list<int>::const_iterator iti = nni.begin();
  while (iti != nni.end()) {
    int i = *iti;
    int mi = (this->isWeighted) ? 1 : this->clusters[i].nMembers;
    std::list<int>::const_iterator itj = nnj.begin();
    while (itj != nnj.end()) {
      int j = *itj;
      int mj = (this->isWeighted) ? 1 : this->clusters[j].nMembers;
      double pij = this->proximity.getValue(i, j);
      double adij = getAlphaProximity(mi, mj, smi, smj, pij);
      alphaij = alphaij + adij;
      itj ++;
    }
    iti ++;
  }
  return alphaij;
}

double mdendro::LanceWilliams::betaTerm(const std::list<int>& nni,
    std::pair<int, int> smi, std::pair<int, int> smj) {
  double betai = 0.0;
  std::list<int>::const_iterator it1 = nni.begin();
  while (it1 != nni.end()) {
    int i1 = *it1;
    int mi1 = (this->isWeighted) ? 1 : this->clusters[i1].nMembers;
    std::list<int>::const_iterator it2 = it1;
    it2 ++;
    while (it2 != nni.end()) {
      int i2 = *it2;
      int mi2 = (this->isWeighted) ? 1 : this->clusters[i2].nMembers;
      double pi = this->proximity.getValue(i1, i2);
      double bdi = getBetaProximity(mi1, mi2, smi, smj, pi);
      betai = betai + bdi;
      it2 ++;
    }
    it1 ++;
  }
  return betai;
}
