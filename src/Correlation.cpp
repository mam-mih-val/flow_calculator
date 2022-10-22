//
// Created by mikhail on 7/30/21.
//

#include "Correlation.hpp"

Correlation::Correlation(TFile* file,
                         const std::string& directory,
                         const std::vector<std::string>& vector_names,
                         const std::vector<std::string>& component_names){
  component_names_ = component_names;
  vector_names_ = vector_names;
  Qn::DataContainerStatCalculate* container{nullptr};
  // Firstly attempting to read the correlation with straightforward vector name composition
  for (const auto &component : component_names) {
    std::string name = directory + "/";
    for( const auto& vec : vector_names ){
      name+="."+vec;
    }
    name+=component;
    file->GetObject(name.c_str(), container);
    if (container) {
      components_.emplace_back(*container);
      continue;
    } else {
      Qn::DataContainerStatCollect *container_collect{nullptr};
      file->GetObject(name.c_str(), container_collect);
      if (container_collect) {
        components_.emplace_back(*container_collect);
        continue;
      }
    }
  }
  if( component_names.size() == components_.size() )
    return;
  auto possible_correlation_names = GetNameCombinations( vector_names );
  for( const auto& possible_name : possible_correlation_names  ) {
    for (const auto &component : component_names) {
      std::string name = directory + "/" + possible_name + "." + component;
      file->GetObject(name.c_str(), container);
      if (container) {
        components_.emplace_back(*container);
        continue;
      } else {
        Qn::DataContainerStatCollect *container_collect{nullptr};
        file->GetObject(name.c_str(), container_collect);
        if (container_collect) {
          components_.emplace_back(*container_collect);
          continue;
        }
      }
    }
    if( component_names.size() == components_.size() )
      return;
  }
  throw std::runtime_error( std::string(__func__) + ": No such container with the name "
                           + possible_correlation_names.front()+": "+possible_correlation_names.size()
                           + " combinations were attempted" );
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
void Correlation::Rebin(std::vector<Qn::AxisD> axes) {
  for( const auto& axis : axes ){
    for( auto& correlation : components_ ){
      correlation = correlation.Rebin( axis );
    }
  }
}
void Correlation::Project(std::vector<std::string> axes) {
  for( auto& correlation : components_ ){
    correlation = correlation.Projection( axes );
  }
}
std::vector<std::string> Correlation::GetNameCombinations(std::vector<std::string> vector_names) {
  if( vector_names.size() == 1 ){
    return vector_names;
  }
  if( vector_names.size() == 2 ){
    std::vector<std::string> combinations;
    auto qa = vector_names.at(0);
    auto qb = vector_names.at(1);
    combinations.emplace_back( qa+"."+qb );
    combinations.emplace_back( qb+"."+qa );
    return combinations;
  }
  if( vector_names.size() > 2 ){
    std::vector<std::string> combinations;
    auto qa = vector_names.front();
    vector_names.erase(vector_names.begin());
    auto prev_combinations = GetNameCombinations( vector_names );
    for( const auto& comb : prev_combinations ){
      combinations.emplace_back( qa + "." +comb );
      combinations.emplace_back( comb + "." + qa );
    }
    return combinations;
  }
  return {};
}