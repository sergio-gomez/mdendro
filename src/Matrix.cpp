#include "Matrix.h"

#include <algorithm>  // std::max, std::min
#include <cmath>  // std::round, std::sqrt
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

mdendro::Matrix::Matrix() {
  this->nrows = 0;
  this->diagvalue = NOT_A_NUMBER;
  this->minvalue = +INF;
  this->maxvalue = -INF;
}

mdendro::Matrix::Matrix(const Matrix& other) {
  this->nrows = other.nrows;
  this->diagvalue = other.diagvalue;
  this->trivalues = other.trivalues;
  this->minvalue = +INF;
  this->maxvalue = -INF;
  for (int i = 0; i < (int)trivalues.size(); i ++) {
    this->minvalue = std::min(this->minvalue, this->trivalues[i]);
    this->maxvalue = std::max(this->maxvalue, this->trivalues[i]);
  }
}

mdendro::Matrix::Matrix(double diagvalue,
    const std::vector<double>& trivalues) {
  int nvalues = (int)trivalues.size();
  this->nrows = (1 + (int)std::round(std::sqrt((double)(1 + 8 * nvalues)))) / 2;
  this->diagvalue = diagvalue;
  this->trivalues = trivalues;
  this->minvalue = +INF;
  this->maxvalue = -INF;
  for (int i = 0; i < (int)trivalues.size(); i ++) {
    this->minvalue = std::min(this->minvalue, this->trivalues[i]);
    this->maxvalue = std::max(this->maxvalue, this->trivalues[i]);
  }
}

mdendro::Matrix::Matrix(int nrows, double diagvalue) {
  this->nrows = nrows;
  this->diagvalue = diagvalue;
  int nvalues = (nrows - 1) * nrows / 2;
  this->trivalues = std::vector<double>(nvalues, NOT_A_NUMBER);
  this->minvalue = +INF;
  this->maxvalue = -INF;
}

int mdendro::Matrix::rows() const {
  return this->nrows;
}

double mdendro::Matrix::getDiagonalValue() const {
  return this->diagvalue;
}

std::vector<double> mdendro::Matrix::getTriangularValues() const {
  return this->trivalues;
}

double mdendro::Matrix::getValue(int i, int j) const {
  double value;
  if (i == j) {
    value = this->diagvalue;
  } else {
    value = this->trivalues[index(i, j)];
  }
  return value;
}

void mdendro::Matrix::setTriangularValue(int i, int j, double trivalue) {
  if (i != j) {
    this->trivalues[index(i, j)] = trivalue;
    this->minvalue = std::min(this->minvalue, trivalue);
    this->maxvalue = std::max(this->maxvalue, trivalue);
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
  for (int i = 0; i < (int)trivalues.size(); i ++) {
    oss.str("");  // Clear string stream
    oss << this->trivalues[i];
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
