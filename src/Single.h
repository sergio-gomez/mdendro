#ifndef MDENDRO_SAHN_VERSATILE_SINGLE_H_
#define MDENDRO_SAHN_VERSATILE_SINGLE_H_

#include "Matrix.h"
#include "Versatile.h"

namespace mdendro {

  class Single : public Versatile {
  public:
    Single();
    Single(const Matrix& proximity, bool isDistance, int precision,
        bool isVariable);
    virtual ~Single();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_SINGLE_H_ */
