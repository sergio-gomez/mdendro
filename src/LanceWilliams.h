#ifndef MDENDRO_SAHN_LANCEWILLIAMS_LANCEWILLIAMS_H_
#define MDENDRO_SAHN_LANCEWILLIAMS_LANCEWILLIAMS_H_

#include "Matrix.h"
#include "Sahn.h"

namespace mdendro {

  class LanceWilliams : public Sahn {
  public:
    LanceWilliams();
    LanceWilliams(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~LanceWilliams();
  protected:
    double alphaTerm(const std::list<int>& nni, const std::list<int>& nnj,
        std::pair<int, int> smi, std::pair<int, int> smj);
    double betaTerm(const std::list<int>& nni, std::pair<int, int> smi,
        std::pair<int, int> smj);
  private:
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj);
    virtual double getAlphaProximity(int mi, int mj, std::pair<int, int> smi,
        std::pair<int, int> smj, double pij) = 0;
    virtual double getBetaProximity(int mi1, int mi2, std::pair<int, int> smi,
        std::pair<int, int> smj, double pi) = 0;
  };

}

#endif /* MDENDRO_SAHN_LANCEWILLIAMS_LANCEWILLIAMS_H_ */
