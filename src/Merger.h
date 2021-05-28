#ifndef MDENDRO_UTILS_MERGER_H_
#define MDENDRO_UTILS_MERGER_H_

#include <list>  // std::list

namespace mdendro {

  class Merger {
  public:
    Merger();
    Merger(double height, int i);
    double getHeight() const;
    double getRange() const;
    void setRange(double range);
    std::list<int> getClusters() const;
    void addCluster(int i);
  private:
    double height;  // Minimum distance or maximum similarity
    double range;  // Difference between maximum and minimum (dis)similarities
    std::list<int> clusters;  // Clusters merged
  };

}

#endif /* MDENDRO_UTILS_MERGER_H_ */
