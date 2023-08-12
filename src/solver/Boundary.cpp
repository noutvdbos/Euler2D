#include <iostream>
#include <algorithm>
#include <cmath>

#include "Boundary.h"
#include "Names.h"
#include "Thermo.h"
#include "utils/Constants.h"


//=======================================================================
//  class Boundary
//=======================================================================


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

Boundary::Boundary()
{}

Boundary::Boundary( int val )
{}

Boundary::~Boundary()
{}


//-----------------------------------------------------------------------
//   getOutletBoundary()
//-----------------------------------------------------------------------


Vector<double>& Boundary::getOutletBoundary
(
    Vector<double>& cell
)
{
    return cell;
}


//-----------------------------------------------------------------------
//   getInletBoundary()
//-----------------------------------------------------------------------


Vector<double> Boundary::getInletBoundary
(
    double gamma
)
{
    double E = Thermo::getInternalEnergy(p, rho, u, v, gamma );

    return { rho, rho*u, rho*v, rho*E};
}
 

//-----------------------------------------------------------------------
//   getWallBoundary()
//-----------------------------------------------------------------------


void Boundary::getWallBoundary
(
    Vector<double>&         bcCell,
    const Vector<double>&   cell,
    const Vector<double>&   normal
)

{    
    bcCell[0] = cell[0];
    bcCell[3] = cell[3];

    Vector<double> vDir = {cell[1], cell[2]};
    vDir = vDir / std::max(vDir.norm(),Constants::VERY_SMALL);

    Vector<double> rDir = vDir - 2*(normal*vDir).sum()*normal; 

    bcCell[1] = cell[1]*rDir[0];
    bcCell[2] = cell[2]*rDir[1];
}  


//-----------------------------------------------------------------------
//   getSubInletBoundary()
//-----------------------------------------------------------------------


void Boundary::getSubInletBoundary
(
    Vector<double>&         bcCell,
    const Vector<double>&   cell,
    const Vector<double>&   normal
)

{    
    bcCell[0] = rho;
    bcCell[1] = rho*u;
    bcCell[2] = rho*v;
    bcCell[3] = cell[3];
}  


//-----------------------------------------------------------------------
//   getSubOutletBoundary()
//-----------------------------------------------------------------------


void Boundary::getSubOutletBoundary
(
    Vector<double>&         bcCell,
    const Vector<double>&   cell,
    const Vector<double>&   normal,
    double                  gamma
)

{    
    bcCell[0] = cell[0];
    bcCell[1] = cell[1];
    bcCell[2] = cell[2];
    bcCell[3] = cell[0]* Thermo::getInternalEnergy( p, cell[0], 
                                                    cell[1]/cell[0], 
                                                    cell[2]/cell[0], gamma );
}  


//-----------------------------------------------------------------------
//   checkDirection()
//-----------------------------------------------------------------------


bool Boundary::checkDirection
(
    const Mesh&           mesh,
    int                   cellIdx,
    int                   node1,
    int                   node2
)
{
    // calculate the location of the ghost cell, xg, yg

    double a = mesh.points[node1][1] - mesh.points[node2][1];
    double b = mesh.points[node2][0] - mesh.points[node1][0];
    double c = mesh.points[node1][0] * mesh.points[node2][1] -
               mesh.points[node1][1] * mesh.points[node2][0];

    double x1 = mesh.cents[cellIdx][0];
    double y1 = mesh.cents[cellIdx][1];
    
    double xg = x1 - 2*(a*x1+b*y1+c)/(a*a+b*b)*a;
    double yg = y1 - 2*(a*x1+b*y1+c)/(a*a+b*b)*b;

    Vector<double> centGhost = {xg, yg};

    Vector<double> normal = mesh.calcNormal(node1,node2);
    std::tuple<int,int> slice = std::make_tuple( 0, normal.rowSize() );

    return ( (normal*mesh.cents(cellIdx,slice)).sum() >
            (normal*centGhost).sum()   );
}


//-----------------------------------------------------------------------
//   operator<<()
//-----------------------------------------------------------------------


std::ostream& operator<<
(
    std::ostream& os,
    const Boundary& bound
)
{
    os << bound.name << ":\n" << bound.key << "\t" << bound.type << "\t";
    
    return os;    
}

