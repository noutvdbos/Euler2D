
#include "Face.h"
#include "Mesh.h"
#include "Names.h"
#include "Boundary.h"
#include "Thermo.h"
#include "utils/utils.h"
#include "utils/Constants.h"

//=======================================================================
//   class Face
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


Face::Face()
{
  
  cellL_                = NULL;
  cellR_                = NULL;
  height_               = 0.0;

  hasNeighbour_         = false;
  bcIsRight_          = false;
}


Face::Face(int val)
{}


Face::Face
(
    const Mesh& mesh,
    const Matrix<double>& cellStates,
    const Config& config,
    int node1,
    int node2,
    int cellidx1,
    int cellidx2,
    int bcidx,
    bool hasNeighbour
)
{
    int idxL = cellidx1;
    int idxR = cellidx2;

    hasNeighbour_    = hasNeighbour;
    bcIsRight_       = false;

    cellL_ = NULL;
    cellR_ = NULL;

    cellBc_.reinit(cellStates.colSize()); 

    flux_ = Flux( config.getGamma(), config.getFluxScheme() );

    normal_ = mesh.calcNormal(node1, node2);
    height_ = (Vector<double> {mesh.points[node1] -
                                mesh.points[node2]}).norm();

    cent_ = (mesh.points[node1] + mesh.points[node2])*0.5;

    if ( hasNeighbour )
    {
        if( mesh.checkDirection(idxL,idxR,node1,node2) )                
        {
            int temp = idxR;
            idxR = idxL;
            idxL = temp;
        }       

        cellL_ = &cellStates[idxL];
        cellR_ = &cellStates[idxR];

    }
    else
    {
        //check which boundary it is
        
        int bcKey = mesh.bcKeys[ bcidx ];

        int idx = 0;
        bound_.key = bcKey*2;
        while (bound_.key != bcKey)
        {
            bound_ = config.getBoundary(idx);
            idx++;
        }

        //Freestream is simply an inlet or outlet, depending on the direction of the flow
        // and the direction of the boundary. So adjust the boundary type here.
        if ( bound_.type == BoundaryTypeNames::FREESTREAM )
        {
            double dotprod = bound_.u*normal_[0] + bound_.v*normal_[1];
            double mach = std::sqrt( (bound_.u*bound_.u+bound_.v*bound_.v) /
                                     (config.getGamma()*bound_.p/bound_.rho) );

            if ( ( bound_.checkDirection(mesh, cellidx1, node1, node2) ) == 
                 ( dotprod > 0 ) )
            {
                if (mach > 1)
                {
                    bound_.type = BoundaryTypeNames::SUPERSONIC_INLET;
                }
                else
                {
                    bound_.type = BoundaryTypeNames::SUBSONIC_INLET;
                }
            }
            else
            {
                if (mach > 1)
                {
                    bound_.type = BoundaryTypeNames::SUPERSONIC_OUTLET;
                }
                else
                {
                    bound_.type = BoundaryTypeNames::SUBSONIC_OUTLET;
                }
            }
        }
       

        if ( bound_.type == BoundaryTypeNames::SUPERSONIC_INLET )
        {
            if ( bound_.checkDirection(mesh, cellidx1, node1, node2) )
            {
                cellL_ =  new Vector<double>( bound_.getInletBoundary(config.getGamma()) );
                cellR_  = &cellStates[cellidx1];
            }
            else
            {
                cellL_  = &cellStates[cellidx1];
                cellR_ =  new Vector<double>(bound_.getInletBoundary(config.getGamma()) );
            }
        }

        else if ( bound_.type == BoundaryTypeNames::SUPERSONIC_OUTLET )
        {
            if ( bound_.checkDirection(mesh, cellidx1, node1, node2) )
            {
                cellL_ = &bound_.getOutletBoundary(cellStates[cellidx1]);
                cellR_ = &cellStates[cellidx1];
            }
            else
            {
                cellL_ = &cellStates[cellidx1];
                cellR_ = &bound_.getOutletBoundary(cellStates[cellidx1]);
            }
        }

        else //for subsonic bc's and walls, it is necessary to check where the bc is.
        {
            if ( bound_.checkDirection(mesh, cellidx1, node1, node2) )
            {
                bcIsRight_ = false;
                cellR_  = &cellStates[cellidx1];               
            }
            else
            { 
                bcIsRight_ = true;
                cellL_  = &cellStates[cellidx1];
            }
        }
    }
}

Face::~Face()
{
}


//-----------------------------------------------------------------------
//   operator=
//-----------------------------------------------------------------------

Face& Face::operator= ( const Face& face )
{
    cellL_ = face.getCellL();
    cellR_ = face.getCellR();

    cellBc_ = face.getCellBc();

    normal_ = face.getNormal();
    cent_   = face.getCentPos();
    height_ = face.getHeight();

    hasNeighbour_    = face.hasNeighbour();
    bcIsRight_     = face.bcIsRight();

    flux_  = face.getFluxObject();
    bound_ = face.getBoundary();

    return *this;
}


//-----------------------------------------------------------------------
//   getFlux();
//-----------------------------------------------------------------------


void Face::getFlux ( Vector<double>& flux, bool print )
{
    
    if ( bound_.type == BoundaryTypeNames::WALL )
    {
        if ( bcIsRight_ )
        {
            Boundary::getWallBoundary(cellBc_,*cellL_, normal_);

            flux_.calcFlux(flux, *cellL_, cellBc_, normal_);
        }
        else
        {
            Boundary::getWallBoundary(cellBc_,*cellR_, normal_);

            flux_.calcFlux(flux, cellBc_, *cellR_, normal_);
        }    
    }
    else if ( bound_.type == BoundaryTypeNames::SUBSONIC_INLET )
    {
        if ( bcIsRight_ )
        {
            bound_.getSubInletBoundary(cellBc_,*cellL_, normal_);

            flux_.calcFlux(flux, *cellL_, cellBc_, normal_);
        }
        else
        {
            bound_.getSubInletBoundary(cellBc_,*cellR_, normal_);

            flux_.calcFlux(flux, cellBc_, *cellR_, normal_);
        }
    }
    else if ( bound_.type == BoundaryTypeNames::SUBSONIC_OUTLET )
    {
        if ( bcIsRight_ )
        {
            bound_.getSubOutletBoundary(cellBc_,*cellL_, normal_, flux_.getGamma() );

            flux_.calcFlux(flux, *cellL_, cellBc_, normal_);
        }
        else
        {
            bound_.getSubOutletBoundary(cellBc_,*cellR_, normal_, flux_.getGamma());

            flux_.calcFlux(flux, cellBc_, *cellR_, normal_);
        }
    }
    else
    {
        flux_.calcFlux(flux, *cellL_, *cellR_, normal_);
    }
    
    flux*= height_;

    if (print)
        cellBc_.print();
}


//-----------------------------------------------------------------------
//   getWallPressure();
//-----------------------------------------------------------------------

void Face::getWallPressure ( Vector<double>& result )
{
    // returns the location of wall and pressure

    assert (bound_.type == BoundaryTypeNames::WALL);

    result[0] = cent_[0];
    result[1] = cent_[1];
    result[2] = cent_[2];

    if ( bcIsRight_ )
    {
        result[3] = Thermo::getPressure(*cellL_, flux_.getGamma());
    }
    else
    {
        result[3] = Thermo::getPressure(*cellR_, flux_.getGamma());
    }
}


//-----------------------------------------------------------------------
//   isWallBoundary();
//-----------------------------------------------------------------------

bool Face::isWallBoundary() const
{
    return ( bound_.type == BoundaryTypeNames::WALL );
}


//-----------------------------------------------------------------------
//   print();
//-----------------------------------------------------------------------


void Face::print (  )
{
    std::cout << "height: " << height_ << "\n";
    std::cout << "hasNeigbour: " << hasNeighbour_ << "\n";
}