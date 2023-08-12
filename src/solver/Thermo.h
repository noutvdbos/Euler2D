#pragma once


//-----------------------------------------------------------------------------
//   class Thermo
//-----------------------------------------------------------------------------


class Thermo
{

public:

  // static functions

  static double getPressure
  (
      const Vector<double>& cellState,
      double                gamma
  );

  static Vector<double> getPressure
  (
      const Matrix<double>& cellStates,
      double                gamma
  );

  static double getRhoEnthalpy
  (
      const Vector<double>& cellState,
      double                gamma
  );

  static  Vector<double> getRhoEnthalpy
  (
      const Matrix<double>& cellStates,
      double                gamma
  );
  
  static double getInternalEnergy
  (
    double  p,
    double  rho,
    double  u,
    double  v,
    double  gamma
  );

  static void getInternalEnergy
  (
    Vector<double>&        E,
    const Vector<double>&  p,
    const Vector<double>&  rho,
    const Vector<double>&  u,
    const Vector<double>&  v,
    double                 gamma
  );

};