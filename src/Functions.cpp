//
// Created by mikhail on 8/1/21.
//

#include "Functions.hpp"
#include "Math/SpecFuncMathMore.h"

ClassImp(Functions);

Correlation Functions::Resolution3S(const Correlation& first, const Correlation& second, const Correlation& third) {
  auto result = Sqrt(first*second/third);
  return result;
}
std::vector<Correlation> Functions::VectorResolutions3S(TFile* file,
                                                    const std::string& directory,
                                                    const std::string& ep_vector,
                                                    const std::vector<std::string>& res_vectors,
                                                    const std::vector<std::string>& comp_names){
  std::vector<Correlation> res_vector;
  const auto& qa = ep_vector;
  auto res_v = res_vectors;
  auto it = std::find( res_v.begin(), res_v.end(), qa );
  if (it != res_v.end()){
    res_v.erase(it);
  }
  for( size_t i=0; i<res_v.size(); ++i ){
    auto qb = res_v.at(i);
    for( size_t j=i+1; j<res_v.size(); ++j ){
      auto qc = res_v.at(j);
      Correlation qa_qb;
      try {
        qa_qb = Correlation (file, directory, {qa, qb}, comp_names);
      } catch ( std::exception& ){
        try {
          qa_qb = Correlation(file, directory, {qb, qa}, comp_names);
        }catch ( std::exception& ){
          continue;
        }
      }
      Correlation qa_qc;
      try {
        qa_qc = Correlation(file, directory, {qa, qc}, comp_names);
      }catch ( std::exception& ){
        try {
          qa_qc = Correlation(file, directory, {qc, qa}, comp_names);
        }catch ( std::exception& ){
          continue;
        }
      }
      Correlation qb_qc;
      try {
        qb_qc = Correlation(file, directory, {qb, qc}, comp_names);
      }catch ( std::exception& ){
        try {
          qb_qc = Correlation(file, directory, {qc, qb}, comp_names);
        }catch ( std::exception& ){
          continue;
        }
      }

      auto res_qa = Functions::Resolution3S( qa_qb*2, qa_qc*2, qb_qc*2 );
      res_qa.SetTitle( qa+"("+qb+","+qc+")" );
      res_vector.emplace_back(res_qa);
    }
  }
  return res_vector;
}
std::vector<Correlation>
Functions::VectorResolutions4S(TFile *file, const std::string &directory,
                               const std::string &ep_vector,
                               const std::vector<std::string> &sub_vectors,
                               const std::vector<std::string> &res_vectors,
                               const std::vector<std::string> &comp_names) {
  std::vector<Correlation> result_vector;
  for( const auto& sub : sub_vectors ){
    Correlation ep_sub;
    try {
      ep_sub = Correlation(file, directory, {ep_vector, sub}, comp_names) * 2;
    } catch (std::exception&){
      ep_sub = Correlation(file, directory, {sub, ep_vector}, comp_names) * 2;
    }
    auto vec_res_sub = VectorResolutions3S( file, directory, sub, res_vectors, comp_names );
    for( const auto& R1_sub : vec_res_sub ){
      auto res_4sub = ep_sub / R1_sub;
      res_4sub.SetTitle( ep_vector+"."+R1_sub.Title() );
      result_vector.emplace_back( res_4sub );
    }
  }
  return result_vector;
}

Correlation
Functions::ExtrapolateToFullEvent(const Correlation &half_event_resolution,
                                  double order) {
  auto result = half_event_resolution;
  for( auto& component : result.GetComponents() ){
    for( auto& bin : component ){
      auto mean = bin.Mean();
      if( fabs( mean ) < std::numeric_limits<double>::min() )
        continue;
      auto chi = DichotomyResolutionSolver( mean, order, {-10.0, 10.0} );
      auto extrapolation = ResolutionFunction( sqrt(2)*chi, order, 0 );
      bin = bin * extrapolation / mean;
    }
  }
  return result;
}
double Functions::DichotomyResolutionSolver(double res,
                                            double order,
                                            std::vector<double> range) {
  auto a = range.at(0);
  auto b = range.at(1);
  while( fabs(a-b) > 1e-3 ){
    auto c = (a + b) / 2;
    auto fa = ResolutionFunction(a, order, res);
    auto fc = ResolutionFunction(c, order, res);
    if( fa*fc < 0 ){
      b = c;
      continue;
    }
    auto fb = ResolutionFunction(b, order, res);
    if( fb*fc < 0 ){
      a = c;
      continue;
    }
    throw std::runtime_error("The case is unsolvable in this range");
  }
  return (a+b)/2;
}
double Functions::ResolutionFunction(double chi, double k, double y) {
  auto chi2_over_2 = chi*chi / 2;
  auto f1 = sqrt( M_PI ) / 2 * chi * exp( -chi2_over_2 );
  auto f2 = ROOT::Math::cyl_bessel_i((k-1)/2, chi2_over_2);
  auto f3 = ROOT::Math::cyl_bessel_i((k+1)/2, chi2_over_2);
  auto f = f1*(f2+f3);
  return f - y;
};
