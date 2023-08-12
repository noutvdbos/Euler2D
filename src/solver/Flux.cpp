
#include "Flux.h"
#include "Mesh.h"
#include "Names.h"
#include "Thermo.h"
#include "utils/utils.h"
#include "utils/Constants.h"


//=======================================================================
//   class Flux
//=======================================================================


//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

Flux::Flux()
{
    gamma_ = 1.4;
    fluxScheme_ = FluxSchemeNames::AUSM_PLUS;
}

Flux::Flux
(
    double gamma,
    std::string fluxScheme
)
{
    gamma_ = gamma;
    fluxScheme_ = fluxScheme;
}

Flux::~Flux()
{}


//-----------------------------------------------------------------------
//   calcFlux()
//-----------------------------------------------------------------------

void Flux::calcFlux
(
    Vector<double>&         flux,
    const Vector<double>&   cell1,
    const Vector<double>&   cell2,
    const Vector<double>&   normal
)

{
    if (fluxScheme_ == FluxSchemeNames::UPWIND)
    {
        calcUpwindFlux_( flux, cell1, cell2 );
        return;
    }
    
    if (fluxScheme_ == FluxSchemeNames::AUSM_PLUS)
    {
        calcAusmPlusFlux_( flux, cell1, cell2, normal);
        return;

    }

    // if it gets here without returning anything, the function failed.
    std::exit(EXIT_FAILURE);
}


//-----------------------------------------------------------------------
//   calcUpwindFlux_()
//-----------------------------------------------------------------------

void Flux::calcUpwindFlux_
(
    Vector<double>&          flux,   
    const Vector<double>&    cellL,
    const Vector<double>&    cellR   
)

{
    // dummy return value for now
    flux = cellL;

    return;
}


//-----------------------------------------------------------------------
//   calcAusmPlusFlux_()
//-----------------------------------------------------------------------


void Flux::calcAusmPlusFlux_
(
    Vector<double>&          flux,   
    const Vector<double>&    cellL,
    const Vector<double>&    cellR,
    const Vector<double>&    normal        
)
{
    double pL = Thermo::getPressure(cellL, gamma_); 
    double pR = Thermo::getPressure(cellR, gamma_);

    double uNormL = (cellL[1]*normal[0]+ cellL[2]*normal[1])/cellL[0];
    double uNormR = (cellR[1]*normal[0]+ cellR[2]*normal[1])/cellR[0];

    double rhoEnthalpyL = Thermo::getRhoEnthalpy(cellL, gamma_);
    double rhoEnthalpyR = Thermo::getRhoEnthalpy(cellR, gamma_);
    
    double aL = getSOS_(rhoEnthalpyL/cellL[0], uNormL);
    double aR = getSOS_(rhoEnthalpyR/cellR[0], uNormR);

    double machL = uNormL/aL;
    double machR = uNormR/aR;    

    double mFace = getFaceMach_(machL, machR);
    double pFace = getFacePressure_(machL,machR,pL,pR);
    

    double m0 = 0.5*(mFace + std::abs(mFace));
    double m1 = 0.5*(mFace - std::abs(mFace));

    flux[0] = m0*aL*cellL[0] + m1*aR*cellR[0];
    flux[1] = m0*aL*cellL[1] + m1*aR*cellR[1] + pFace*normal[0];
    flux[2] = m0*aL*cellL[2] + m1*aR*cellR[2] + pFace*normal[1];
    flux[3] = m0*aL*rhoEnthalpyL + m1*aR*rhoEnthalpyR;

    return;
}


//-----------------------------------------------------------------------
//  getSOS_()
//-----------------------------------------------------------------------


double Flux::getSOS_
(
    double      enthalpy,
    double      uNorm
)

{
    double aCrit = std::sqrt(2*(gamma_-1)/(gamma_+1)*enthalpy);
    
    return std::pow(aCrit,2)/(std::max(aCrit,std::abs(uNorm)));
}


//-----------------------------------------------------------------------
//  getSOS_()
//-----------------------------------------------------------------------


Vector<double> Flux::getSOS_
(
    const Vector<double>&      enthalpy,
    const Vector<double>&      uNorm
)

{
    assert( enthalpy.rowSize() == uNorm.rowSize() );

    Vector<double> vec(enthalpy.rowSize());

    for (int i = 0; i < enthalpy.rowSize(); i++)
    {
        vec[i] = getSOS_(enthalpy[i],uNorm[i]);
    }
    
    return vec;
}


//-----------------------------------------------------------------------
//  getFaceMach_()
//-----------------------------------------------------------------------


double Flux::getFaceMach_
(
    const Vector<double>&     machs
)

{
    using std::abs;
    using std::pow;

    double mplus;
    double mmin;

    if (std::abs(machs[0]) <= 1 )
    {
        mplus = 0.25*pow(machs[0]+1, 2)  + 0.125*pow( pow(machs[0],2)-1, 2);
    }
    else
    {
        mplus = 0.5*(machs[0] + abs(machs[0])); 
    }

    if ( abs(machs[1])  <= 1 )
    {
        mmin  = -0.25*pow(machs[1]-1,2) - 0.125*pow( pow(machs[1],2)-1, 2);
    }
    else
    {
        mmin  = 0.5*(machs[1] - abs(machs[1]));
    }

    
    return mplus + mmin;
}


//-----------------------------------------------------------------------
//  getFaceMach_()
//-----------------------------------------------------------------------


double Flux::getFaceMach_
(
    double machL,
    double machR
)

{
    using std::abs;
    using std::pow;

    double mplus;
    double mmin;

    if (std::abs(machL) <= 1 )
    {
        mplus = 0.25*pow(machL+1, 2)  + 0.125*pow( pow(machL,2)-1, 2);
    }
    else
    {
        mplus = 0.5*(machL + abs(machL)); 
    }

    if ( abs(machR)  <= 1 )
    {
        mmin  = -0.25*pow(machR-1,2) - 0.125*pow( pow(machR,2)-1, 2);
    }
    else
    {
        mmin  = 0.5*(machR - abs(machR));
    }

    
    return mplus + mmin;
}


//-----------------------------------------------------------------------
//  getFacePressure_()
//-----------------------------------------------------------------------


double Flux::getFacePressure_
(
    const Vector<double>&     machs,
    const Vector<double>&     ps
)

{
    using std::abs;
    using std::pow;

    double pplus;
    double pmin;

    if ( abs(machs[0]) <= 1 )
    {
        pplus = 0.25*pow(machs[0]+1,2)*(2-machs[0]) +
                3/16*machs[0]*pow( pow(machs[0],2)-1,2);
    }
    else
    {
        pplus = 0.5*(machs[0] + abs(machs[0])) / machs[0];
    }

    if ( abs(machs[1]) <= 1 )
    {
        pmin  = 0.25*pow(machs[1]-1,2)*(2+machs[1]) -
                3/16*machs[1]*pow( pow(machs[1],2)-1,2);
    }
    else
    {
        pmin  = 0.5*(machs[1] - abs(machs[1])) / machs[1];
    }
    
    return ps[0]*pplus + ps[1]*pmin;
}


//-----------------------------------------------------------------------
//  getFacePressure_()
//-----------------------------------------------------------------------


double Flux::getFacePressure_
(
    double machL,
    double machR,
    double pL,
    double pR
)

{
    using std::abs;
    using std::pow;

    double pplus;
    double pmin;

    if ( abs(machL) <= 1 )
    {
        pplus = 0.25*pow(machL+1,2)*(2-machL) +
                3/16*machL*pow( pow(machL,2)-1,2);
    }
    else
    {
        pplus = 0.5*(machL + abs(machL)) / machL;
    }

    if ( abs(machR) <= 1 )
    {
        pmin  = 0.25*pow(machR-1,2)*(2+machR) -
                3/16*machR*pow( pow(machR,2)-1,2);
    }
    else
    {
        pmin  = 0.5*(machR - abs(machR)) / machR;
    }
    
    return pL*pplus + pR*pmin;
}



