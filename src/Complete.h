#ifndef MDENDRO_SAHN_VERSATILE_COMPLETE_H_
#define MDENDRO_SAHN_VERSATILE_COMPLETE_H_

#include "Matrix.h"
#include "Versatile.h"

namespace mdendro {

  class Complete : public Versatile {
  public:
    Complete();
    Complete(const Matrix& proximity, bool isDistance, int precision,
        bool isVariable);
    virtual ~Complete();
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
  };

}

#endif /* MDENDRO_SAHN_VERSATILE_COMPLETE_H_ */
