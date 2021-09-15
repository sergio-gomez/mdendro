#include "Ultrametricity.h"

#include <algorithm>  // std::max, std::min
#include <cmath>  // std::log, std::sqrt
#include <list>  // std::list
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"

mdendro::Ultrametricity::Ultrametricity() {
  this->nObjects = 0;
  this->cophCorr = NOT_A_NUMBER;
  this->distortion = NOT_A_NUMBER;
  this->agglomerative = NOT_A_NUMBER;
  this->chaining = NOT_A_NUMBER;
  this->balance = NOT_A_NUMBER;
}

mdendro::Ultrametricity::Ultrametricity(const Matrix& iniProx,
	const std::vector<Merger>& mergers, double bottomHgt) {
  this->nObjects = iniProx.rows();
  calcCopheneticProximity(mergers);
  calcCopheneticMeasures(iniProx);
  calcAgglomerativeMeasures(mergers, bottomHgt);
}

mdendro::Matrix mdendro::Ultrametricity::getCopheneticProximity() const {
  return this->cophProx;
}

double mdendro::Ultrametricity::getCopheneticCorrCoeff() const {
  return this->cophCorr;
}

double mdendro::Ultrametricity::getSpaceDistortion() const {
  return this->distortion;
}

double mdendro::Ultrametricity::getAgglomerativeCoeff() const {
  return this->agglomerative;
}

double mdendro::Ultrametricity::getChainingCoeff() const {
  return this->chaining;
}

double mdendro::Ultrametricity::getTreeBalance() const {
  return this->balance;
}

void mdendro::Ultrametricity::calcCopheneticProximity(
    const std::vector<Merger>& mergers) {
  this->cophProx = Matrix(this->nObjects);
  std::vector< std::list<int> > members(this->nObjects);
  for (int i = 0; i < (int)members.size(); i ++) {
    members[i].push_back(i);
  }
  for (int k = 0; k < (int)mergers.size(); k ++) {
    double prox = mergers[k].getHeight();
    std::list<int> clusters = mergers[k].getClusters();
    std::list<int>::const_iterator iti = clusters.begin();
    while (iti != clusters.end()) {
      int ci = *iti;
      std::list<int>::const_iterator itj = iti;
      itj ++;
      while (itj != clusters.end()) {
        int cj = *itj;
        groupPair(members[ci], members[cj], prox);
        itj ++;
      }
      iti ++;
    }
    //  Add members
    iti = clusters.begin();
    int ci = *iti;
    iti ++;
    while (iti != clusters.end()) {
      int cj = *iti;
      members[ci].splice(members[ci].end(), members[cj]);
      iti ++;
    }
  }
}

void mdendro::Ultrametricity::groupPair(const std::list<int>& sci,
    const std::list<int>& scj, double prox) {
  std::list<int>::const_iterator iti = sci.begin();
  while (iti != sci.end()) {
    int i = *iti;
    std::list<int>::const_iterator itj = scj.begin();
    while (itj != scj.end()) {
      int j = *itj;
      this->cophProx.setValue(i, j, prox);
      itj ++;
    }
    iti ++;
  }
}

void mdendro::Ultrametricity::calcCopheneticMeasures(const Matrix& iniProx) {
  const double INF = std::numeric_limits<double>::infinity();
  double pMin = +INF;
  double pMax = -INF;
  double pSum = 0.0;
  double ppSum = 0.0;
  double cMin = +INF;
  double cMax = -INF;
  double cSum = 0.0;
  double ccSum = 0.0;
  double pcSum = 0.0;
  for (int i = 0; i < iniProx.rows(); i ++) {
    for (int j = i + 1; j < iniProx.rows(); j ++) {
      double pij = iniProx.getValue(i, j);
      pMin = std::min(pMin, pij);
      pMax = std::max(pMax, pij);
      pSum += pij;
      ppSum += pij * pij;
      double cij = this->cophProx.getValue(i, j);
      cMin = std::min(cMin, cij);
      cMax = std::max(cMax, cij);
      cSum += cij;
      ccSum += cij * cij;
      pcSum += pij * cij;
    }
  }
  int nValues = (iniProx.rows() - 1) * iniProx.rows() / 2;
  ppSum *= (double)nValues;
  ccSum *= (double)nValues;
  pcSum *= (double)nValues;
  double ppSum2 = pSum * pSum;
  double ccSum2 = cSum * cSum;
  double pcSum2 = pSum * cSum;
  double ppSigma2 = ppSum - ppSum2;
  double ccSigma2 = ccSum - ccSum2;
  double pcSigma2 = pcSum - pcSum2;
  this->cophCorr = pcSigma2 / std::sqrt(ppSigma2 * ccSigma2);
  this->distortion = (cMax - cMin) / (pMax - pMin);
}

void mdendro::Ultrametricity::calcAgglomerativeMeasures(
    const std::vector<Merger>& mergers, double bottomHgt) {
  double sumHgts = 0.0;
  int diffMembers = 0;
  this->balance = 0.0;
  std::vector<int> nMembers(this->nObjects, 1);
  for (int k = 0; k < (int)mergers.size(); k ++) {
    int maxMembers = 0;
    int minMembers = this->nObjects;
    int sumMembers = 0;
    std::list<int> clusters = mergers[k].getClusters();
    std::list<int>::const_iterator it = clusters.begin();
    while (it != clusters.end()) {
      int i = *it;
      if (nMembers[i] == 1) {
        sumHgts += mergers[k].getHeight() - bottomHgt;
      }
      maxMembers = std::max(maxMembers, nMembers[i]);
      minMembers = std::min(minMembers, nMembers[i]);
      sumMembers += nMembers[i];
      it ++;
    }
    diffMembers += maxMembers - minMembers;
    this->balance += entropy(clusters, sumMembers, nMembers);
    nMembers[clusters.front()] = sumMembers;
  }
  double topHgt = mergers.back().getHeight() - bottomHgt;
  this->agglomerative = 1.0 - sumHgts / ((double)this->nObjects * topHgt);
  if (this->nObjects < 3) {
    this->chaining = 0.0;
  } else {
    int maxDiff = (this->nObjects - 2) * (this->nObjects - 1) / 2;
    this->chaining = (double)diffMembers / (double)maxDiff;
  }
  this->balance /= (double)mergers.size();
}

double mdendro::Ultrametricity::entropy(const std::list<int>& clusters,
    int sMembers, const std::vector<int>& nMembers) const {
  double h = 0.0;
  std::list<int>::const_iterator it = clusters.begin();
  while (it != clusters.end()) {
    int i = *it;
    double p = (double)nMembers[i] / (double)sMembers;
    h += - p * std::log(p);
    it ++;
  }
  h /= std::log((double)clusters.size());
  return h;
}
