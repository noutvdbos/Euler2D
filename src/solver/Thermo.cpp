

#include "utils/Vector.h"
#include "utils/Matrix.h"
#include "utils/Constants.h"
#include "Thermo.h"


//=======================================================================
//   static class Thermo
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   getPressure()
//-----------------------------------------------------------------------

double Thermo::getPressure
(
    const Vector<double>& cellState,
    double         gamma
)
{
    return (gamma-1)*( cellState[3] - 0.5/cellState[0]
        *( std::pow(cellState[1],2)+std::pow(cellState[2],2) ));
}


//-----------------------------------------------------------------------
//   getPressure()
//-----------------------------------------------------------------------


Vector<double> Thermo::getPressure
(
    const Matrix<double>& cellStates,
    double         gamma
)
{
    Vector<double> vec ( cellStates.rowSize());
    for ( int i = 0; i < vec.rowSize(); i++)
    {
        vec[i] = getPressure(cellStates[i], gamma);
    }
    return vec;
}


//-----------------------------------------------------------------------
//   getRhoEnthalpy()
//-----------------------------------------------------------------------


double Thermo::getRhoEnthalpy
(
    const Vector<double>& cellState,
    double         gamma
)
{    
    double p = getPressure(cellState, gamma);
    
    return cellState[3] + p;
}


//-----------------------------------------------------------------------
//   getRhoEnthalpy()
//-----------------------------------------------------------------------


Vector<double> Thermo::getRhoEnthalpy
(
    const Matrix<double>& cellStates,
    double         gamma
)
{
    Vector<double> vec ( cellStates.rowSize());
    for ( int i = 0; i < vec.rowSize(); i++)
    {
        vec[i] = getRhoEnthalpy(cellStates[i], gamma);
    }
    return vec;
}


//-----------------------------------------------------------------------
//   getInternalEnergy()
//-----------------------------------------------------------------------


double  Thermo::getInternalEnergy
(
    double p,
    double rho,
    double u,
    double v,
    double gamma
)
{    
    return p/((gamma-1)*rho) + 0.5*( std::pow(u,2)+std::pow(v,2) );
}


//-----------------------------------------------------------------------
//   getInternalEnergy()
//-----------------------------------------------------------------------


void  Thermo::getInternalEnergy
(
    Vector<double>&        E,
    const Vector<double>&  p,
    const Vector<double>&  rho,
    const Vector<double>&  u,
    const Vector<double>&  v,
    double                 gamma
)
{    
    E.reinit( rho.rowSize() );

    for (int i = 0; i< E.rowSize(); i++ )
    {
        E[i] = p[i]/(std::max( (gamma-1)*rho[i], Constants::VERY_SMALL) )
             + 0.5*( std::pow(u[i],2)+std::pow(v[i],2) );
    }

}
