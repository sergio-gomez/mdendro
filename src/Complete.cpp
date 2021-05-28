#include "Complete.h"

#include "Matrix.h"
#include "Versatile.h"

mdendro::Complete::Complete() :
  Versatile() {
  this->power = +INF;
}

mdendro::Complete::Complete(const Matrix& proximity, bool isDistance,
    int precision, bool isVariable) :
  Versatile(+INF, false, proximity, isDistance, precision, isVariable) {}

mdendro::Complete::~Complete() {}

double mdendro::Complete::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  return maximumProximity(nni, nnj);
}
