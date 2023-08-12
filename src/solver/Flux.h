#pragma once

#include <string>

#include "utils/Vector.h"
#include "Mesh.h"

//-----------------------------------------------------------------------
//   class Flux
//-----------------------------------------------------------------------


class Flux
{   
public:

  //Constructor
  Flux();
  Flux(double gamma, std::string fluxScheme );

  //Destructor
  ~Flux();

  // member functions
  void calcFlux
  (
      Vector<double>&         flux,    
      const Vector<double>&   cellL,
      const Vector<double>&   cellR,
      const Vector<double>&   normal
 );

 double getGamma() { return gamma_;}

private:

  void calcUpwindFlux_
  (
      Vector<double>&         flux,
      const Vector<double>&   cellL,
      const Vector<double>&   cellR
  );

  void calcAusmPlusFlux_
  (
      Vector<double>&         flux,
      const Vector<double>&   cellL,
      const Vector<double>&   cellR,
      const Vector<double>&   normal
  );


  double         getSOS_
  (
    double              enthalpy,
    double              uMag
  );  
  
  Vector<double> getSOS_
  (
    const Vector<double>&      enthalpy,
    const Vector<double>&      uMag
  );

  double        getFaceMach_
  (
    const Vector<double>& machs
  );

  double        getFaceMach_
  (
    double machL,
    double machR
  );

  double        getFacePressure_
  (
    const Vector<double>& machs,
    const Vector<double>& ps
  );

  double        getFacePressure_
  (
    double machL,
    double machR,
    double pL,
    double pR
  );

private:

  double gamma_;
  std::string fluxScheme_;

};