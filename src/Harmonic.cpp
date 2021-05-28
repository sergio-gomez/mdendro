#include "Harmonic.h"

#include "Matrix.h"
#include "Versatile.h"

mdendro::Harmonic::Harmonic() :
  Versatile() {
  this->power = -1.0;
}

mdendro::Harmonic::Harmonic(bool isWeighted, const Matrix& proximity,
    bool isDistance, int precision, bool isVariable) :
  Versatile(-1.0, isWeighted, proximity, isDistance, precision, isVariable) {}

mdendro::Harmonic::~Harmonic() {}

double mdendro::Harmonic::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  return generalizedMean(nni, nnj);
}
