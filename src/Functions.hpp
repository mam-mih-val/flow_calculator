//
// Created by mikhail on 8/1/21.
//

#ifndef FLOW_CALCULATOR_SRC_FUNCTIONS_HPP_
#define FLOW_CALCULATOR_SRC_FUNCTIONS_HPP_

#include <TObject.h>
#include <TFile.h>
#include "Correlation.hpp"

class Functions : public TObject {
public:
  Functions() = default;
  ~Functions() override = default;
  static Correlation Resolution3S(const Correlation&, const Correlation&, const Correlation&);
  static std::vector<Correlation> VectorResolutions3S(TFile* file,
                                                      const std::string& directory,
                                                      const std::string& ep_vector,
                                                      const std::vector<std::string>& res_vectors,
                                                      const std::vector<std::string>& comp_names);
private:
  ClassDefOverride(Functions, 1);
};

#endif // FLOW_CALCULATOR_SRC_FUNCTIONS_HPP_
