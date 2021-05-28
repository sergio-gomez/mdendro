#ifndef MDENDRO_UTILS_MATRIX_H_
#define MDENDRO_UTILS_MATRIX_H_

#include <limits>  // std::numeric_limits
#include <vector>  // std::vector

namespace mdendro {

  const double INF = std::numeric_limits<double>::infinity();
  const int MAX_DIGITS = std::numeric_limits<double>::digits10;
  const double NOT_A_NUMBER = std::numeric_limits<double>::quiet_NaN();

  // Symmetric matrix with null diagonal values
  class Matrix {
  public:
    Matrix();
    Matrix(const Matrix& other);
    Matrix(double diagvalue, const std::vector<double>& trivalues);
    Matrix(int nrows, double diagvalue);
    int rows() const;
    double getDiagonalValue() const;
	std::vector<double> getTriangularValues() const;
    double getValue(int i, int j) const;
    void setTriangularValue(int i, int j, double trivalue);
    double getMinimumValue() const;
    double getMaximumValue() const;
    int getPrecision() const;
  private:
    int nrows;  // Number of rows
    double diagvalue;  // Diagonal value
    std::vector<double> trivalues;  // Lower triangular values by columns
    double minvalue;  // Lower triangular minimum value
    double maxvalue;  // Lower triangular maximum value
    int index(int i, int j) const;
  };

}

#endif /* MDENDRO_UTILS_MATRIX_H_ */
