//
// Created by mikhail on 7/30/21.
//

#ifndef QNANALYSIS_SRC_QNANALYSISCALCULATE_CORRELATION_HPP_
#define QNANALYSIS_SRC_QNANALYSISCALCULATE_CORRELATION_HPP_

#include <QnTools/DataContainer.hpp>
#include <TFile.h>
#include <TObject.h>

class Correlation : public TObject {
 public:
  Correlation() = default;
  Correlation(TFile* file,
              const std::string& directory,
              const std::vector<std::string>& vector_names,
              const std::vector<std::string>& component_names);
  Correlation(const Correlation&) = default;
  ~Correlation() override = default;
  Qn::DataContainerStatCalculate operator[]( size_t idx ){return components_.at(idx);};
  [[nodiscard]] const std::string &Title() const { return title_; }
  [[nodiscard]] const std::vector<std::string> &GetComponentNames() const {
    return component_names_;
  }
  [[nodiscard]] const std::vector<std::string> &GetVectorNames() const {
    return vector_names_;
  }
  size_t Size(){ return components_.size(); };
  void SetTitle(const std::string &titles) { title_ = titles; }
  void Save( const std::string& name ){
    int idx=0;
    for( auto component : components_ ){
      component.Write( std::data(name+"."+component_names_.at(idx)) );
      idx++;
    }
  }
  friend Correlation operator+( const Correlation&, const Correlation& );
  friend Correlation operator-( const Correlation&, const Correlation& );
  friend Correlation operator*( const Correlation&, const Correlation& );
  friend Correlation operator*( const Correlation&, const Correlation& );
  friend Correlation operator*( const Correlation&, const std::vector<double>& );
  friend Correlation operator*( const Correlation&, double );
  friend Correlation operator/( const Correlation&, const Correlation& );
  friend Correlation operator/( const Correlation&, double );
  friend Correlation Sqrt( const Correlation& );
  friend Correlation MatrixMultiply( const Correlation&, const Correlation& );
protected:
  std::vector<Qn::DataContainerStatCalculate> components_;
  std::string title_;
  std::vector<std::string> component_names_;
  std::vector<std::string> vector_names_;
  ClassDefOverride(Correlation, 1)
};

#endif//QNANALYSIS_SRC_QNANALYSISCALCULATE_CORRELATION_HPP_
