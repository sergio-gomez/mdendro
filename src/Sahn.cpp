#include "Sahn.h"

#include <algorithm>  // std::find, std::max, std::min
#include <cmath>  // std::abs, std::floor, std::log10, std::pow, std::round
#include <iterator>  // std::next
#include <list>  // std::list
#include <stack>  // std::stack
#include <utility>  // std::pair
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"

mdendro::Sahn::Cluster::Cluster() {
  this->prevAgglomerable = -1;
  this->nextAgglomerable = -1;
  this->nMembers = 1;
  this->proximity = 0.0;
}

bool mdendro::Sahn::Cluster::hasNearestNeighbor(int i) const {
  std::list<int>::const_iterator it = std::find(this->nearestNeighbors.begin(),
      this->nearestNeighbors.end(), i);
  bool found = (it != this->nearestNeighbors.end());
  return found;
}

mdendro::Sahn::Sahn() {
  this->isWeighted = false;
  this->isVariable = true;
  this->nObjects = 0;
  this->isDistance = true;
  this->epsilon = std::pow(1.0, -(double)MAX_DIGITS);
  this->pow10precision = 1e6;
  this->firstAgglomerable = -1;
}

mdendro::Sahn::Sahn(bool isWeighted, const Matrix& proximity, bool isDistance,
    int precision, bool isVariable) {
  // Initializations
  this->isWeighted = isWeighted;
  this->isVariable = isVariable;
  this->proximity = Matrix(proximity);
  this->nObjects = proximity.rows();
  this->isDistance = isDistance;
  double maxProx = std::max(std::abs(proximity.getMaximumValue()), 1.0);
  int intDigits = 1 + (int)std::floor(std::log10(maxProx));
  int maxPrecision = MAX_DIGITS - intDigits;
  this->epsilon = std::pow(10.0, -(double)maxPrecision);
  // 0 <= precision <= maxPrecision
  precision = std::max(precision, 0);
  precision = std::min(precision, maxPrecision);
  this->pow10precision = std::pow(10.0, (double)precision);
  this->mergers.reserve(this->nObjects - 1);
  // Initial partition
  this->clusters = std::vector<Cluster>(this->nObjects);
  for (int i = 0; i < this->nObjects; i ++) {
    this->clusters[i].prevAgglomerable = i - 1;
    this->clusters[i].nextAgglomerable = i + 1;
    this->clusters[i].proximity = isDistance? +INF : -INF;
  }
  this->firstAgglomerable = 0;
  for (int i = 0; i < this->nObjects - 1; i ++) {
    setNearestNeighbors(i);
  }
}

mdendro::Sahn::~Sahn() {}

void mdendro::Sahn::build() {
  int nAgglomerated = 0;
  // Repeat while there are clusters to agglomerate
  while (nAgglomerated < this->nObjects - 1) {
    double pnext;
    std::list<int> inext;
    getNextProximity(pnext, inext);
    std::vector<bool> connected = connectNeighbours(pnext, inext);
    int nNew = createAgglomerations(pnext, inext);
    nAgglomerated = nAgglomerated + nNew;
    updateProximity(inext, connected);
    updateNeighbors(connected);
  }
}

std::vector<mdendro::Merger> mdendro::Sahn::getMergers() const {
  return this->mergers;
}

std::pair<int, int> mdendro::Sahn::sumMembers(const std::list<int>& nn) const {
  int sm1;
  int sm2;
  if (this->isWeighted) {
    sm1 = nn.size();
    sm2 = nn.size();
  } else {
    sm1 = 0;
    sm2 = 0;
    std::list<int>::const_iterator it = nn.begin();
    while (it != nn.end()) {
      int i = *it;
      int mi = this->clusters[i].nMembers;
      sm1 = sm1 + mi;
      sm2 = sm2 + mi * mi;
      it ++;
    }
  }
  return std::pair<int, int>(sm1, sm2);
}

void mdendro::Sahn::getNextProximity(double& pnext, std::list<int>& inext) {
  pnext = this->isDistance? +INF : -INF;
  int i = this->firstAgglomerable;
  while (i < this->nObjects) {
    double pi = precisionRound(this->clusters[i].proximity);
    if (( this->isDistance && (pi < pnext)) ||
        (!this->isDistance && (pi > pnext))) {
      inext.clear();
      inext.push_back(i);
      pnext = pi;
    } else if ((pi == pnext) && this->isVariable) {
      inext.push_back(i);
    }
    i = this->clusters[i].nextAgglomerable;
  }
}

std::vector<bool> mdendro::Sahn::connectNeighbours(double pnext,
    std::list<int>& inext) {
  std::vector<bool> connected = std::vector<bool>(this->nObjects, false);
  std::list<int>::iterator itnext = inext.begin();
  while (itnext != inext.end()) {
    int ii = *itnext;
    if (connected[ii]) {
      itnext = inext.erase(itnext);
    } else {
      connected[ii] = true;
      std::list<int> nnNext;
      std::list<int>::const_iterator itnn =
          this->clusters[ii].nearestNeighbors.begin();
      while (itnn != this->clusters[ii].nearestNeighbors.end()) {
        int j = *itnn;
        connectComponent(j, pnext, connected, nnNext);
        itnn ++;
      }
      clearNearestNeighbors(ii);
      itnn = nnNext.begin();
      while (itnn != nnNext.end()) {
        int j = *itnn;
        this->clusters[ii].nearestNeighbors.push_back(j);
        this->clusters[j].nearestNeighborOf.push_back(ii);
        itnn ++;
      }
      itnext ++;
    }
  }
  return connected;
}

void mdendro::Sahn::connectComponent(int jj, double pnext,
    std::vector<bool>& connected, std::list<int>& nnNext) {
  std::stack<int> st;
  st.push(jj);
  while (!st.empty()) {
    int j = st.top();
    st.pop();
    if (!connected[j]) {
      nnNext.push_back(j);
      connected[j] = true;
      double pj = precisionRound(this->clusters[j].proximity);
      removeAgglomerable(j);
      if (this->isVariable) {
        if (pj == pnext) {
          std::list<int>::const_iterator itnn =
              this->clusters[j].nearestNeighbors.begin();
          while (itnn != this->clusters[j].nearestNeighbors.end()) {
            int k = *itnn;
            st.push(k);
            itnn ++;
          }
        }
        std::list<int>::iterator itnnof =
            this->clusters[j].nearestNeighborOf.begin();
        while (itnnof != this->clusters[j].nearestNeighborOf.end()) {
          int i = *itnnof;
          double pi = precisionRound(this->clusters[i].proximity);
          if (pi == pnext) {
            st.push(i);
          }
          itnnof ++;
        }
      }
    }
  }
}

void mdendro::Sahn::removeAgglomerable(int j) {
  int i = this->clusters[j].prevAgglomerable;
  int k = this->clusters[j].nextAgglomerable;
  if (i < 0) {
    this->firstAgglomerable = k;
  } else {
    this->clusters[i].nextAgglomerable = k;
  }
  if (k < this->nObjects) {
    this->clusters[k].prevAgglomerable = i;
  }
  this->clusters[j].prevAgglomerable = -1;
  this->clusters[j].nextAgglomerable = -1;
  this->clusters[j].proximity = 0.0;
}

int mdendro::Sahn::createAgglomerations(double pnext,
    const std::list<int>& inext) {
  int nNew = 0;
  std::list<int>::const_iterator itnext = inext.begin();
  while (itnext != inext.end()) {
    int ii = *itnext;
    Merger merger(pnext, ii);
    std::list<int>::const_iterator itnn =
        this->clusters[ii].nearestNeighbors.begin();
    while (itnn != this->clusters[ii].nearestNeighbors.end()) {
      int i = *itnn;
      merger.addCluster(i);
      itnn ++;
    }
    // Calculate aggloremation range
    double range = 0.0;
    std::list<int> clus = merger.getClusters();
    std::list<int>::const_iterator iti = clus.begin();
    while (iti != clus.end()) {
      int ci = *iti;
      std::list<int>::const_iterator itj = iti;
      itj ++;
      while (itj != clus.end()) {
        int cj = *itj;
        double pij = precisionRound(this->proximity.getValue(ci, cj));
        range = std::max(range, std::abs(pij - pnext));
        itj ++;
      }
      iti ++;
    }
    merger.setRange(range);
    this->mergers.push_back(merger);
    nNew += this->clusters[ii].nearestNeighbors.size();
    itnext ++;
  }
  return nNew;
}

void mdendro::Sahn::setNearestNeighbors(int i) {
  int jnext = -1;
  double pnext = this->isDistance? +INF : -INF;
  int j = this->clusters[i].nextAgglomerable;
  while (j < this->nObjects) {
    double pij = precisionRound(this->proximity.getValue(i, j));
    if (( this->isDistance && (pij < pnext)) ||
        (!this->isDistance && (pij > pnext))) {
      jnext = j;
      pnext = pij;
    }
    j = this->clusters[j].nextAgglomerable;
  }
  clearNearestNeighbors(i);
  this->clusters[i].proximity = pnext;
  if (jnext > -1) {
    if (this->isVariable) {
      int j = jnext;
      while (j < this->nObjects) {
        double pij = precisionRound(this->proximity.getValue(i, j));
        if (pij == pnext) {
          this->clusters[i].nearestNeighbors.push_back(j);
          this->clusters[j].nearestNeighborOf.push_back(i);
        }
        j = this->clusters[j].nextAgglomerable;
      }
    } else {
      this->clusters[i].nearestNeighbors.push_back(jnext);
      this->clusters[jnext].nearestNeighborOf.push_back(i);
    }
  }
}

void mdendro::Sahn::clearNearestNeighbors(int i) {
  std::list<int>::const_iterator it =
    this->clusters[i].nearestNeighbors.begin();
  while (it != this->clusters[i].nearestNeighbors.end()) {
    int j = *it;
    this->clusters[j].nearestNeighborOf.remove(i);
    it ++;
  }
  this->clusters[i].nearestNeighbors.clear();
}

void mdendro::Sahn::updateProximity(const std::list<int>& inext,
    const std::vector<bool>& connected) {
  std::list<int>::const_iterator iti = inext.begin();
  while (iti != inext.end()) {
    int ii = *iti;
    std::list<int> nni = nearestNeighbors(connected, ii);
    std::list<int>::const_iterator itj = iti;
    itj ++;
    while (itj != inext.end()) {
      int jj = *itj;
      std::list<int> nnj = nearestNeighbors(connected, jj);
      double pij = newProximity(nni, nnj);
      this->proximity.setTriangularValue(ii, jj, pij);
      itj ++;
    }
    int k = this->firstAgglomerable;
    while (k < this->nObjects) {
      if (!connected[k]) {
        std::list<int> nnk = nearestNeighbors(connected, k);
        double pik = newProximity(nni, nnk);
        this->proximity.setTriangularValue(ii, k, pik);
        if (k < ii) {
          // Ensure that correct nearest neighbours are found when
          // non-monotone (e.g. centroid) methods are used
          pik = precisionRound(pik);
          double pk = precisionRound(this->clusters[k].proximity);
          if (pik <= pk) {
            if (pik < pk) {
              clearNearestNeighbors(k);
              this->clusters[k].proximity = pik;
            }
            this->clusters[k].nearestNeighbors.push_back(ii);
            this->clusters[ii].nearestNeighborOf.push_back(k);
          }
        }
      }
      k = this->clusters[k].nextAgglomerable;
    }
    std::list<int>::const_iterator itnn =
        this->clusters[ii].nearestNeighbors.begin();
    while (itnn != this->clusters[ii].nearestNeighbors.end()) {
      int i = *itnn;
      this->clusters[ii].nMembers += this->clusters[i].nMembers;
      clearNearestNeighbors(i);
      itnn ++;
    }
    setNearestNeighbors(ii);
    iti ++;
  }
}

void mdendro::Sahn::updateNeighbors(const std::vector<bool>& connected) {
  int i = this->firstAgglomerable;
  while (i < this->nObjects) {
    bool hasnn = false;
    std::list<int>::const_iterator it =
        this->clusters[i].nearestNeighbors.begin();
    while ((!hasnn) && (it != this->clusters[i].nearestNeighbors.end())) {
      int j = *it;
      hasnn = hasnn || connected[j];
      it ++;
    }
    if (hasnn) {
      // Redetermine NN of i
      setNearestNeighbors(i);
    }
  i = this->clusters[i].nextAgglomerable;
  }
}

std::list<int> mdendro::Sahn::nearestNeighbors(
    const std::vector<bool>& connected, int i) {
  std::list<int> nn;
  nn.push_back(i);
  if (connected[i]) {
    std::list<int>::const_iterator it =
        this->clusters[i].nearestNeighbors.begin();
    while (it != this->clusters[i].nearestNeighbors.end()) {
      int j = *it;
      nn.push_back(j);
      it ++;
    }
  }
  return nn;
}

double mdendro::Sahn::precisionRound(double value) const {
  // Add epsilon to avoid 0.49999999999999... being rounded to 0
  value += (value >= 0.0)? +this->epsilon : -this->epsilon;
  return std::round(value * this->pow10precision) / this->pow10precision;
}
