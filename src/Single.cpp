#include "Single.h"

#include "Matrix.h"
#include "Versatile.h"

mdendro::Single::Single() :
  Versatile() {
  this->power = -INF;
}

mdendro::Single::Single(const Matrix& proximity, bool isDistance, int precision,
    bool isVariable) :
  Versatile(-INF, false, proximity, isDistance, precision, isVariable) {}

mdendro::Single::~Single() {}

double mdendro::Single::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  return minimumProximity(nni, nnj);
}
