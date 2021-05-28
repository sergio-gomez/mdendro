#ifndef MDENDRO_SAHN_VERSATILE_ARITHMETIC_H_
#define MDENDRO_SAHN_VERSATILE_ARITHMETIC_H_

#include "Matrix.h"
#include "Versatile.h"

namespace mdendro {

  class Arithmetic : public Versatile {
  public:
    Arithmetic();
    Arithmetic(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~Arithmetic();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_ARITHMETIC_H_ */
