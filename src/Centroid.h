#ifndef MDENDRO_SAHN_LANCEWILLIAMS_CENTROID_H_
#define MDENDRO_SAHN_LANCEWILLIAMS_CENTROID_H_

#include "LanceWilliams.h"
#include "Matrix.h"

namespace mdendro {

  class Centroid : public LanceWilliams {
  public:
    Centroid();
    Centroid(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~Centroid();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
    virtual double getAlphaProximity(int mi, int mj, std::pair<int, int> smi,
        std::pair<int, int> smj, double pij);
    virtual double getBetaProximity(int mi1, int mi2, std::pair<int, int> smi,
        std::pair<int, int> smj, double pi);
  };

}

#endif /* MDENDRO_SAHN_LANCEWILLIAMS_CENTROID_H_ */
