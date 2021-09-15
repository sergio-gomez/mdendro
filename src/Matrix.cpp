#include "Matrix.h"

#include <algorithm>  // std::max, std::min
#include <cmath>  // std::round, std::sqrt
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

mdendro::Matrix::Matrix() {
  this->nrows = 0;
  this->minvalue = +INF;
  this->maxvalue = -INF;
}

mdendro::Matrix::Matrix(const Matrix& other) {
  this->nrows = other.nrows;
  this->values = other.values;
  this->minvalue = +INF;
  this->maxvalue = -INF;
  for (int i = 0; i < (int)values.size(); i ++) {
    this->minvalue = std::min(this->minvalue, this->values[i]);
    this->maxvalue = std::max(this->maxvalue, this->values[i]);
  }
}

mdendro::Matrix::Matrix(const std::vector<double>& values) {
  int nvalues = (int)values.size();
  this->nrows = (1 + (int)std::round(std::sqrt((double)(1 + 8 * nvalues)))) / 2;
  this->values = values;
  this->minvalue = +INF;
  this->maxvalue = -INF;
  for (int i = 0; i < (int)values.size(); i ++) {
    this->minvalue = std::min(this->minvalue, this->values[i]);
    this->maxvalue = std::max(this->maxvalue, this->values[i]);
  }
}

mdendro::Matrix::Matrix(int nrows) {
  this->nrows = nrows;
  int nvalues = (nrows - 1) * nrows / 2;
  this->values = std::vector<double>(nvalues, NOT_A_NUMBER);
  this->minvalue = +INF;
  this->maxvalue = -INF;
}

int mdendro::Matrix::rows() const {
  return this->nrows;
}

std::vector<double> mdendro::Matrix::getValues() const {
  return this->values;
}

double mdendro::Matrix::getValue(int i, int j) const {
  double value;
  if (i == j) {
    value = NOT_A_NUMBER;
  } else {
    value = this->values[index(i, j)];
  }
  return value;
}

void mdendro::Matrix::setValue(int i, int j, double value) {
  if (i != j) {
    this->values[index(i, j)] = value;
    this->minvalue = std::min(this->minvalue, value);
    this->maxvalue = std::max(this->maxvalue, value);
  }
}

double mdendro::Matrix::getMinimumValue() const {
  return this->minvalue;
}

double mdendro::Matrix::getMaximumValue() const {
  return this->maxvalue;
}

int mdendro::Matrix::getPrecision() const {
  std::ostringstream oss;
  oss.precision(MAX_DIGITS);  // Modify the default precision
  int maxdecimals = 0;
  for (int i = 0; i < (int)values.size(); i ++) {
    oss.str("");  // Clear string stream
    oss << this->values[i];
    std::string s = oss.str();
    std::size_t found = s.find('.');
    int decimals = (found == std::string::npos)? 0 : (int)(s.size()-found)-1;
    maxdecimals = std::max(maxdecimals, decimals);
  }
  return maxdecimals;
}

int mdendro::Matrix::index(int i, int j) const {
  int k;
  if (i == j) {
    k = -1;
  } else if (i > j) {
    k = i + j * this->nrows - (j + 1) * (j + 2) / 2;
  } else { // (i < j)
    k = j + i * this->nrows - (i + 1) * (i + 2) / 2;
  }
  return k;
}
