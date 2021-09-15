#include <algorithm>  // std::max
#include <cmath>  // std::ceil, std::floor, std::log10, std::pow
#include <sstream>  // std::ostringstream
#include <string>  // std::string
#include <vector>  // std::vector

#include <Rcpp.h>

#include "Arithmetic.h"
#include "Centroid.h"
#include "Complete.h"
#include "Flexible.h"
#include "Geometric.h"
#include "Harmonic.h"
#include "Matrix.h"
#include "Merger.h"
#include "Sahn.h"
#include "Single.h"
#include "Ultrametricity.h"
#include "Versatile.h"
#include "Ward.h"

mdendro::Sahn* newLinkage(std::string method, double methodPar, bool isWeighted,
    const mdendro::Matrix& proxMatr, bool isDistance, int precision,
    bool isVariable);

Rcpp::List listMergers(int nObjects,
    const std::vector<mdendro::Merger>& mergers);

// [[Rcpp::export]]
Rcpp::List rcppLinkage(const Rcpp::NumericVector& prox, bool isDistance = true,
    int digits = -1, std::string method = "arithmetic",
    double methodPar = 0.0, bool isWeighted = false, bool isVariable = true) {
  mdendro::Matrix proxMatr(Rcpp::as< std::vector<double> >(prox));
  if (digits < 0) {
    digits = proxMatr.getPrecision();
  } else {
    // Check maximum precision
    double maxProx = std::max(proxMatr.getMaximumValue(), 1.0);
    int intDigits = 1 + (int)std::floor(std::log10(maxProx));
    int maxPrecision = mdendro::MAX_DIGITS - intDigits;
    if (digits > maxPrecision) {
      std::ostringstream oss;
      oss << "'digits' for these data must be less than or equal to "
          << maxPrecision;
      Rcpp::stop(oss.str());
    }
  }
  // Build agglomerative hierarchical clustering
  mdendro::Sahn* sahn = newLinkage(method, methodPar, isWeighted, proxMatr,
      isDistance, digits, isVariable);
  sahn->build();
  std::vector<mdendro::Merger> mergers = sahn->getMergers();
  delete sahn;
  // Save results
  int nObjects = proxMatr.rows();
  Rcpp::List lm = listMergers(nObjects, mergers);
  double bottomHgt = isDistance?
      0.0 : std::pow(10.0, std::ceil(std::log10(proxMatr.getMaximumValue())));
  mdendro::Ultrametricity ultra(proxMatr, mergers, bottomHgt);
  mdendro::Matrix cophMatr = ultra.getCopheneticProximity();
  Rcpp::NumericVector coph = Rcpp::wrap(cophMatr.getValues());
  Rcpp::List lnk = Rcpp::List::create(
      Rcpp::Named("digits") = digits,
      Rcpp::Named("merger") = lm["merger"],
      Rcpp::Named("height") = lm["height"],
      Rcpp::Named("range") = lm["range"],
      Rcpp::Named("coph") = coph,
      Rcpp::Named("cor") = ultra.getCopheneticCorrCoeff(),
      Rcpp::Named("sdr") = ultra.getSpaceDistortion(),
      Rcpp::Named("ac") = ultra.getAgglomerativeCoeff(),
      Rcpp::Named("cc") = ultra.getChainingCoeff(),
      Rcpp::Named("tb") = ultra.getTreeBalance());
  return lnk;
}

mdendro::Sahn* newLinkage(std::string method, double methodPar, bool isWeighted,
    const mdendro::Matrix& proxMatr, bool isDistance, int precision,
    bool isVariable) {
  mdendro::Sahn* sahn;
  if (method == "single") {
    sahn = new mdendro::Single(proxMatr, isDistance, precision, isVariable);
  } else if (method == "complete") {
    sahn = new mdendro::Complete(proxMatr, isDistance, precision, isVariable);
  } else if (method == "arithmetic") {
    sahn = new mdendro::Arithmetic(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  } else if (method == "geometric") {
    sahn = new mdendro::Geometric(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  } else if (method == "harmonic") {
    sahn = new mdendro::Harmonic(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  } else if (method == "versatile") {
    double power = methodPar;
    sahn = new mdendro::Versatile(power, isWeighted, proxMatr, isDistance,
        precision, isVariable);
  } else if (method == "flexible") {
    double beta = methodPar;
    sahn = new mdendro::Flexible(beta, isWeighted, proxMatr, isDistance,
        precision, isVariable);
  } else if (method == "ward") {
    sahn = new mdendro::Ward(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  } else if (method == "centroid") {
    sahn = new mdendro::Centroid(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  } else {  // default method
    sahn = new mdendro::Arithmetic(isWeighted, proxMatr, isDistance, precision,
        isVariable);
  }
  return sahn;
}

Rcpp::List listMergers(int nObjects,
    const std::vector<mdendro::Merger>& mergers) {
  Rcpp::NumericVector height = Rcpp::NumericVector::create();
  Rcpp::NumericVector range = Rcpp::NumericVector::create();
  Rcpp::List merger = Rcpp::List::create();
  std::vector<int> merged = std::vector<int>(nObjects);
  for (int i = 0; i < nObjects; i ++) {
    merged[i] = -(i+1);
  }
  for (int k = 0; k < (int)mergers.size(); k ++) {
    height.push_back(mergers[k].getHeight());
    range.push_back(mergers[k].getRange());
    Rcpp::IntegerVector kMerger = Rcpp::IntegerVector::create();
    std::list<int> clusters = mergers[k].getClusters();
    std::list<int>::const_iterator it = clusters.begin();
    while (it != clusters.end()) {
      int i = *it;
      kMerger.push_back(merged[i]);
      merged[i] = k+1;
      it ++;
    }
    merger.push_back(kMerger);
  }
  return Rcpp::List::create(
      Rcpp::Named("height") = height,
      Rcpp::Named("range") = range,
      Rcpp::Named("merger") = merger);
}
