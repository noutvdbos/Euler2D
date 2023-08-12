#pragma once

#include <string>
#include <iostream>
#include "utils/Vector.h"
#include "Mesh.h"


//-----------------------------------------------------------------------
//   class Boundary
//-----------------------------------------------------------------------


class Boundary
{

public:

  // static variables
  
  // constructor
  Boundary(); 

  Boundary( int val );

  // destructor
  ~Boundary();

  // member variables
  std::string name;
  int key;
  std::string type;
  double u;
  double v;
  double p;
  double rho;

  Vector<double>& getOutletBoundary
  (
    Vector<double>& cell
  );

  Vector<double> getInletBoundary
  (
    double                gamma
  );

  static void getWallBoundary
  (
    Vector<double>&       bcCell,
    const Vector<double>& cell,
    const Vector<double>& normal
  );    

  void getSubInletBoundary
  (
    Vector<double>&       bcCell,
    const Vector<double>& cell,
    const Vector<double>& normal
  );  

  void getSubOutletBoundary
  (
    Vector<double>&       bcCell,
    const Vector<double>& cell,
    const Vector<double>& normal,
    double                gamma
  );  

  bool checkDirection
  (
    const Mesh&           mesh,
    int                   cellIdx,
    int                   node1,
    int                   node2
  );

  friend std::ostream& operator<<(std::ostream& os, const Boundary& bound);

};
