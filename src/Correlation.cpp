//
// Created by mikhail on 7/30/21.
//

#include "Correlation.hpp"

ClassImp(Correlation);

Correlation::Correlation(TFile* file,
                         const std::string& directory,
                         const std::vector<std::string>& vector_names,
                         const std::vector<std::string>& component_names){
  component_names_ = component_names;
  vector_names_ = vector_names;
  Qn::DataContainerStatCalculate* container{nullptr};
  for( const auto& component : component_names ) {
    std::string name = directory + "/";
    for (const auto &vector : vector_names) {
      name += vector + ".";
    }
    name+=component;
    file->GetObject(name.c_str(), container);
    if( container ){
      components_.emplace_back(*container);
    } else{
      Qn::DataContainerStatCollect* container_collect{nullptr};
      file->GetObject(name.c_str(), container_collect);
      if( container_collect ) {
        components_.emplace_back(*container_collect);
      } else
        throw std::runtime_error( "Correlation::Correlation(): No such correlation "+name+" in file" );
    }
  }
}

Correlation operator+(const Correlation& first, const Correlation& second) {
  auto result = Correlation(first);
  if( first.components_.size() != second.components_.size() )
    throw std::runtime_error( "Correlation::operator+(): unequal size of vectors" );
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)+second.components_.at(i);
  return result;
}
Correlation operator-(const Correlation& first, const Correlation& second) {
  auto result = Correlation(first);
  if( first.components_.size() != second.components_.size() )
    throw std::runtime_error( "Correlation::operator-(): unequal size of vectors" );
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)-second.components_.at(i);
  return result;
}
Correlation operator*(const Correlation& first, const Correlation& second) {
  auto result = Correlation(first);
  if( first.components_.size() != second.components_.size() )
    throw std::runtime_error( "Correlation::operator*(): unequal size of vectors" );
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)*second.components_.at(i);
  return result;
}
Correlation operator*(const Correlation& first, const std::vector<double>& tensor) {
  auto result = Correlation(first);
  if( first.components_.size() != tensor.size() )
    throw std::runtime_error( "Correlation::operator*(): unequal size of vectors" );
  for( size_t i=0; i<first.components_.size(); i++ ) {
    result.components_.at(i) = first.components_.at(i) * tensor.at(i);
  }
  return result;
}
Correlation operator*(const Correlation& first, double num) {
  auto result = Correlation(first);
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)*num;
  return result;
}
Correlation operator/(const Correlation& first, const Correlation& second) {
  auto result = Correlation(first);
  if( first.components_.size() != second.components_.size() )
    throw std::runtime_error( "Correlation::operator/(): unequal size of vectors" );
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)/second.components_.at(i);
  return result;
}
Correlation operator/(const Correlation& first, double den) {
  auto result = Correlation(first);
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = first.components_.at(i)/den;
  return result;
}
Correlation Sqrt(const Correlation& first) {
  auto result = Correlation(first);
  for( size_t i=0; i<first.components_.size(); i++ )
    result.components_.at(i) = Sqrt(first.components_.at(i));
  return result;
}

Correlation MatrixMultiply( const Correlation& first, const Correlation& second ){
  Correlation result;
  int idx1=0;
  for( const auto& comp1 : first.components_ ){
    int idx2=0;
    for( const auto& comp2 : second.components_ ){
      result.components_.emplace_back( comp1*comp2 );
      result.component_names_.emplace_back( first.component_names_.at(idx1)+"."+second.component_names_.at(idx2) );
      idx2++;
    }
    idx1++;
  }
  return result;
}