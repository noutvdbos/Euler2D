#pragma once

#include <string>

#include "utils/Vector.h"
#include "Mesh.h"
#include "Config.h"
#include "Flux.h"
#include "Boundary.h"

//-----------------------------------------------------------------------
//   class Face
//-----------------------------------------------------------------------


class Face
{   
public:

  //Constructor
  Face ();
  Face (int val);

  Face
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
  );

  //Destructor
  ~Face();

public:

  Face&      operator=      ( const Face& face );

  void       getFlux        
  ( 
      Vector<double>& flux,
      bool print
  );

  void       getWallPressure( Vector<double>& result );
  bool       isWallBoundary() const;

  void       print();

  Vector<double>* getCellL()          const {return cellL_;}
  Vector<double>* getCellR()          const {return cellR_;}
  Vector<double>  getCellBc()         const {return cellBc_;}
  Vector<double>  getNormal()         const {return normal_;}
  Vector<double>  getCentPos()        const {return cent_;}
  double          getHeight()         const {return height_;}
  Flux            getFluxObject()     const {return flux_;}
  Boundary        getBoundary()       const {return bound_;}
  bool            hasNeighbour()      const {return hasNeighbour_;}
  bool            bcIsRight()         const {return bcIsRight_;}


private:

  //Do not delete these vectors, they point to outside entities
  Vector<double>*   cellL_;
  Vector<double>*   cellR_;
  Vector<double>    cellBc_;

  Vector<double>    normal_;
  Vector<double>    cent_;
  double            height_;

  bool              hasNeighbour_;
  bool              bcIsRight_;

  Flux              flux_;
  Boundary          bound_;


};