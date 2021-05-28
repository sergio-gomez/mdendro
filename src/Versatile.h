#ifndef MDENDRO_SAHN_VERSATILE_VERSATILE_H_
#define MDENDRO_SAHN_VERSATILE_VERSATILE_H_

#include "Matrix.h"
#include "Sahn.h"

namespace mdendro {

  class Versatile : public Sahn {
  public:
    Versatile();
    Versatile(double power, bool isWeighted, const Matrix& proximity,
        bool isDistance, int precision, bool isVariable);
    virtual ~Versatile();
  protected:
    double power;  // Power parameter for the generalized mean
    double minimumProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
    double maximumProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
    double geometricMean(const std::list<int>& nni,
        const std::list<int>& nnj);
    double generalizedMean(const std::list<int>& nni,
        const std::list<int>& nnj);
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_VERSATILE_H_ */
