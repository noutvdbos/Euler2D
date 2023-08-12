#pragma once

#include <string>
#include <iostream>
#include "utils/Vector.h"
#include "Boundary.h"


//-----------------------------------------------------------------------
//   class Config
//-----------------------------------------------------------------------


class Config 
{

public:

  // static variables
  
  // constructor
  Config(); 

  Config

    ( const std::string&      path,
      const std::string&      fileName );

  // destructor
  ~Config    ();

  // member functions

  std::string getMeshFile() const            { return meshFile_; }
  double      getStartTime() const           { return startTime_; }
  double      getEndTime() const             { return endTime_; }
  double      getTimeStep() const            { return dt_; }  
  bool        adjustTimeStep() const         { return adjustTimeStep_; }
  double      getMaxCourant() const          { return maxCourant_; }

  std::string getInitMethod() const          { return initMethod_; }
  std::string getInitFile() const            { return initFile_; }
  double      getInitDensity() const         { return initRho_; }
  double      getInitPressure() const        { return initP_; }
  double      getInitUVelocity() const       { return initU_; }
  double      getInitVVelocity() const       { return initV_; }
  double      getGamma() const               { return gamma_; }
  
  double      getRefDensity() const          { return refRho_; }
  double      getRefPressure() const         { return refP_; }
  double      getRefVelocity() const         { return refVmag_; }
  double      getRefLength() const           { return refLength_; }

  std::string getFluxScheme() const          { return fluxScheme_; }
  std::string getTimeScheme() const          { return timeScheme_; }
  std::string getResultsFolder() const       { return resultsFolder_; }
  std::string getResultsPrefix() const       { return resultsPrefix_; }

  std::string getWriteMethod() const         { return writeMethod_; }
  double      getWriteInterval() const       { return writeInterval_; }

  bool        getWriteCoefficients() const   { return writeCoefficients_; }
  int         getSelBoundarySize() const     { return selectedBoundaries_.rowSize(); }
  std::string getSelBoundary(int row) const  { return selectedBoundaries_[row]; }

  Boundary    getBoundary(int row) const     { return boundaries_[row]; }

  void getInitVariables   
  (
    Vector<double>& p,
    Vector<double>& rho,
    Vector<double>& u,
    Vector<double>& v,
    int         rowSize
  );

private: 
  // member functions

  void readConfigFile_    (std::ifstream& file);
  void readBoundaryFile_  ();
  void readSelectedBCs_   (std::string line);
  std::string toLower_    (std::string s); 
  std::string peekWord_   (std::ifstream& file);

  void ltrim_             (std::string &s); 
  void rtrim_             (std::string &s);
  void trim_              (std::string &s); 


  std::string path_;
  std::string meshFile_;
  std::string boundaryFile_;

  double startTime_;
  double endTime_;
  double dt_;
  bool   adjustTimeStep_;
  double maxCourant_;

  std::string initMethod_;
  std::string initFile_;

  double initRho_;
  double initP_;
  double initU_;
  double initV_;
  double initAlpha_;
  double gamma_;

  double refRho_;
  double refP_;
  double refVmag_;
  double refLength_;

  std::string fluxScheme_;
  std::string timeScheme_;

  std::string resultsFolder_; 		
  std::string resultsPrefix_; 		

  std::string writeMethod_;
  
  double writeInterval_; 	

  bool writeCoefficients_;
  Vector<std::string> selectedBoundaries_;


  Vector<Boundary> boundaries_;

};

