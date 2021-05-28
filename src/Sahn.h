#ifndef MDENDRO_SAHN_SAHN_H_
#define MDENDRO_SAHN_SAHN_H_

#include <list>  // std::list
#include <utility>  // std::pair
#include <vector>  // std::vector

#include "Matrix.h"
#include "Merger.h"

namespace mdendro {

  // Sequential agglomerative hierarchical non-overlapping (SAHN) clustering
  class Sahn {
  public:
    Sahn();
    Sahn(bool isWeighted, const Matrix& proximity, bool isDistance,
        int precision, bool isVariable);
    virtual ~Sahn();
    void build();
    std::vector<Merger> getMergers() const;
  protected:
    class Cluster {
    public:
      Cluster();
      int prevAgglomerable;  // Previous agglomerable cluster
      int nextAgglomerable;  // Next agglomerable cluster
      int nMembers;  // Cluster cardinality
      double proximity;  // (Dis)similarity with nearest neighbors TO THE RIGHT
      std::list<int> nearestNeighbors;  // Nearest neighbors TO THE RIGHT
      std::list<int> nearestNeighborOf;  // Clusters that current one is NN of
      bool hasNearestNeighbor(int i) const;
    };
    bool isWeighted;  // Weighted clustering method
    Matrix proximity;  // (Dis)similarities between objects
    bool isDistance;  // distance (or similarity) data type
    std::vector<Cluster> clusters;  // Clusters
    std::pair<int, int> sumMembers(const std::list<int>& nn) const;
  private:
    bool isVariable;  // Variable (or pair) grouping criterion
    int nObjects;  // Number of objects being clustered
    double epsilon;  // Very small number
    double pow10precision;  // 10 to the power of significant decimal digits
    int firstAgglomerable;  // First agglomerable cluster
    std::vector<Merger> mergers;  // History of mergers
    void getNextProximity(double& pnext, std::list<int>& inext);
    std::vector<bool> connectNeighbours(double pnext, std::list<int>& inext);
    void connectComponent(int jj, double pnext, std::vector<bool>& connected,
        std::list<int>& nnNext);
    void removeAgglomerable(int j);
    int createAgglomerations(double pnext, const std::list<int>& inext);
    void setNearestNeighbors(int i);
    void clearNearestNeighbors(int i);
    void updateProximity(const std::list<int>& inext,
        const std::vector<bool>& connected);
    void updateNeighbors(const std::vector<bool>& connected);
    std::list<int> nearestNeighbors(const std::vector<bool>& connected, int i);
    double precisionRound(double value) const;
    virtual double newProximity(const std::list<int>& nni,
        const std::list<int>& nnj) = 0;
  };

}

#endif /* MDENDRO_SAHN_SAHN_H_ */
