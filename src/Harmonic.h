#ifndef MDENDRO_SAHN_VERSATILE_HARMONIC_H_
#define MDENDRO_SAHN_VERSATILE_HARMONIC_H_

#include "Matrix.h"
#include "Versatile.h"

namespace mdendro {

  class Harmonic : public Versatile {
  public:
    Harmonic();
    Harmonic(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~Harmonic();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_HARMONIC_H_ */
