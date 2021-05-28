#include "Arithmetic.h"

#include "Matrix.h"
#include "Versatile.h"

mdendro::Arithmetic::Arithmetic() :
  Versatile() {
  this->power = +1.0;
}

mdendro::Arithmetic::Arithmetic(bool isWeighted, const Matrix& proximity,
    bool isDistance, int precision, bool isVariable) :
  Versatile(+1.0, isWeighted, proximity, isDistance, precision, isVariable) {}

mdendro::Arithmetic::~Arithmetic() {}

double mdendro::Arithmetic::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  return generalizedMean(nni, nnj);
}
