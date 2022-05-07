//
// Created by mikhail on 8/1/21.
//

#include "Functions.hpp"

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
  auto qa = ep_vector;
  auto res_v = res_vectors;
  auto it = std::find( res_v.begin(), res_v.end(), qa );
  if (it != res_v.end()){
    res_v.erase(it);
  }
  for( size_t i=0; i<res_v.size(); ++i ){
    auto qb = res_v.at(i);
    for( size_t j=i+1; j<res_v.size(); ++j ){
      auto qc = res_v.at(j);
      Correlation qa_qb( file, directory, {qa, qb}, comp_names);
      Correlation qa_qc( file, directory, {qa, qc}, comp_names);
      Correlation qb_qc( file, directory, {qb, qc}, comp_names);

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
  for( auto sub : sub_vectors ){
    auto ep_sub = Correlation( file, directory, {ep_vector, sub}, comp_names ) * 2;
    auto vec_res_sub = VectorResolutions3S( file, directory, sub, res_vectors, comp_names );
    for( auto R1_sub : vec_res_sub ){
      auto res_4sub = ep_sub / R1_sub;
      res_4sub.SetTitle( ep_vector+"."+R1_sub.Title() );
      result_vector.emplace_back( res_4sub );
    }
  }
  return result_vector;
}
