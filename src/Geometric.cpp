#include "Geometric.h"

#include "Matrix.h"
#include "Versatile.h"

mdendro::Geometric::Geometric() :
  Versatile() {
  this->power = 0.0;
}

mdendro::Geometric::Geometric(bool isWeighted, const Matrix& proximity,
    bool isDistance, int precision, bool isVariable) :
  Versatile(0.0, isWeighted, proximity, isDistance, precision, isVariable) {}

mdendro::Geometric::~Geometric() {}

double mdendro::Geometric::newProximity(const std::list<int>& nni,
    const std::list<int>& nnj) {
  return geometricMean(nni, nnj);
}
