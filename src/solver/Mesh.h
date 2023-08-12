#pragma once

#include <fstream>
#include <string>
#include "utils/Matrix.h"


//-----------------------------------------------------------------------
//   class Mesh
//-----------------------------------------------------------------------


class Mesh
{   
  public:

    static const char* UNSTRUCTURED_GRID;
    static const int   MAX_NODES;
    //Constructor
    Mesh                    ();
    Mesh                    ( std::string& meshFile );

    //Destructor
    ~Mesh                   ();

    
    bool              hasNeighbour(int row, int col ) const;
    Vector<double>    calcNormal(int node1, int node2) const;
    bool              checkDirection(int idx1, int idx2, int node1, int node2) const;

    void              convertCellsToPoints
    ( 
      Vector<double>&       pointVals,
      const Vector<double>& cellVals 
    ) const;

  private:

    void readMesh_          ( std::string& meshFile );
    bool checkMesh_         ( std::ifstream& file );
    void readPoints_        ( std::ifstream& file );
    void readElems_         ( std::ifstream& file );
    void readElemTypes_     ( std::ifstream& file );
    void readPhysicTypes_   ( std::ifstream& file );

    void groupElems_       ();
    void calcCentroids_    ();
    void calcSurfaces_     ();
    void calcVertElems_    ();
    void calcVertBcs_      ();
    void calcFaceElems_    ();
    void calcFaceBcs_      ();

  public:

  Matrix<double> points;
  Matrix<double> cents;

  Matrix<int>    cells;
  Matrix<int>    elems;
  Matrix<int>    bcElems;
  Matrix<int>    elsuel;
  Matrix<int>    bcsuel;


  Vector<double> surfaces;
  
  Vector<int>    elsup1;
  Vector<int>    elsup2;
  Vector<int>    bcsup1;
  Vector<int>    bcsup2;
  Vector<int>    cellTypes;
  Vector<int>    physicTypes;  
  Vector<int>    bcKeys;

};

