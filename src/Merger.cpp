#include "Merger.h"

#include <limits>  // std::numeric_limits
#include <list>  // std::list

mdendro::Merger::Merger() {
  const double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();
  this->height = NOT_A_NUMBER;
  this->range = NOT_A_NUMBER;
}

mdendro::Merger::Merger(double height, int i) {
  this->height = height;
  this->range = 0.0;
  this->clusters.push_back(i);
}

double mdendro::Merger::getHeight() const {
  return this->height;
}

double mdendro::Merger::getRange() const {
  return this->range;
}

void mdendro::Merger::setRange(double range) {
  this->range = range;
}

std::list<int> mdendro::Merger::getClusters() const {
  return this->clusters;
}

void mdendro::Merger::addCluster(int i) {
  this->clusters.push_back(i);
}
