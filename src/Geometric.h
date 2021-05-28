#ifndef MDENDRO_SAHN_VERSATILE_GEOMETRIC_H_
#define MDENDRO_SAHN_VERSATILE_GEOMETRIC_H_

#include "Matrix.h"
#include "Versatile.h"

namespace mdendro {

  class Geometric : public Versatile {
  public:
    Geometric();
    Geometric(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~Geometric();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_GEOMETRIC_H_ */
