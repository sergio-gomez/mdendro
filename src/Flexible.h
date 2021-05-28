#ifndef MDENDRO_SAHN_LANCEWILLIAMS_FLEXIBLE_H_
#define MDENDRO_SAHN_LANCEWILLIAMS_FLEXIBLE_H_

#include "LanceWilliams.h"
#include "Matrix.h"

namespace mdendro {

  class Flexible : public LanceWilliams {
  public:
    Flexible();
    Flexible(double beta, bool isWeighted, const Matrix& proximity,
        bool isDistance, int precision, bool isVariable);
    virtual ~Flexible();
  private:
    double beta;  // Beta parameter
    virtual double getAlphaProximity(int mi, int mj, std::pair<int, int> smi,
        std::pair<int, int> smj, double pij);
    virtual double getBetaProximity(int mi1, int mi2, std::pair<int, int> smi,
        std::pair<int, int> smj, double pi);
  };

}

#endif /* MDENDRO_SAHN_LANCEWILLIAMS_FLEXIBLE_H_ */
