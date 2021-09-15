#ifndef MDENDRO_UTILS_ULTRAMETRICITY_H_
#define MDENDRO_UTILS_ULTRAMETRICITY_H_

#include <list>  // std::list
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"

namespace mdendro {

  class Ultrametricity {
  public:
    Ultrametricity();
    Ultrametricity(const Matrix& iniProx, const std::vector<Merger>& mergers,
        double bottomHgt);
    Matrix getCopheneticProximity() const;
    double getCopheneticCorrCoeff() const;
    double getSpaceDistortion() const;
    double getAgglomerativeCoeff() const;
    double getChainingCoeff() const;
    double getTreeBalance() const;
  private:
    int nObjects;  // Number of objects being clustered
    Matrix cophProx;  // Cophenetic (dis)similarities
    double cophCorr;  // Cophenetic Correlation Coefficient
    double distortion;  // Space Distortion Ratio
    double agglomerative;  // Agglomerative Coefficient
    double chaining;  // Chaining Coefficient
    double balance;  // Tree Balance
    void calcCopheneticProximity(const std::vector<Merger>& mergers);
    void groupPair(const std::list<int>& sci, const std::list<int>& scj,
        double prox);
    void calcCopheneticMeasures(const Matrix& iniProx);
    void calcAgglomerativeMeasures(const std::vector<Merger>& mergers,
        double bottomHgt);
    double entropy(const std::list<int>& clusters, int sMembers,
        const std::vector<int>& nMembers) const;
  };

}

#endif /* MDENDRO_UTILS_ULTRAMETRICITY_H_ */
