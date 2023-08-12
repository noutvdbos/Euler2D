#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <algorithm> 
#include <tuple>

#include "Mesh.h"
#include "Names.h"
#include "utils/Matrix.h"
#include "utils/Vector.h"
#include "utils/utils.h"
#include "utils/Types.h"

//=======================================================================
//   class Mesh
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------

const char* Mesh::UNSTRUCTURED_GRID = "UNSTRUCTURED_GRID";
const int   Mesh::MAX_NODES         = 9;

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

Mesh::Mesh()
{

}

//TODO: Add function that checks for duplicate faces and gives warning for this.

Mesh::Mesh( std::string& meshFile )
{
    
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;  

    

    std::cout << 
    "\n---------------------------MESH-----------------------------\n";

    std::cout << "Using mesh file: " << meshFile << "\n\n";

    std::cout << "---Reading the mesh file--- \n";
    auto t0 = high_resolution_clock::now();
    auto t1 = t0;

    readMesh_( meshFile );

    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> time = t2 - t1;
    std::cout << "-This took: " << time.count()/1000 << "s.\n\n";  
    

    std::cout << "---Calculating the centroids--- \n";
    t1 = high_resolution_clock::now();

    calcCentroids_();

    t2 = high_resolution_clock::now();
    time = t2 - t1;
    std::cout << "-This took: " << time.count()/1000 << "s.\n\n";   

    std::cout << "---Calculating the element surfaces--- \n";
    t1 = high_resolution_clock::now();

    calcSurfaces_();

    t2 = high_resolution_clock::now();
    time = t2 - t1;
    std::cout << "-This took: " << time.count()/1000 << "s.\n\n";  


    std::cout << "---Calculating the vertex-elements connections--- \n";
    t1 = high_resolution_clock::now();

    calcVertElems_();
    calcVertBcs_();

    t2 = high_resolution_clock::now();
    time = t2 - t1;
    std::cout << "-This took: " << time.count()/1000 << "s.\n\n";  
    

    std::cout << "---Calculating the face-elements connections--- \n";
    t1 = high_resolution_clock::now();

    calcFaceElems_();
    calcFaceBcs_();

    t2 = high_resolution_clock::now();
    time = t2 - t1;
    std::cout << "-This took: " << time.count()/1000 << "s.\n\n";  

    time = t2 - t0;
    std::cout << "-Total mesh preprocessing time: " << time.count()/1000 << "s.\n\n";
}


Mesh::~Mesh()
{}


//-----------------------------------------------------------------------
//   hasNeighbour()
//-----------------------------------------------------------------------


bool Mesh::hasNeighbour( int row, int col ) const
{
    if ( row >= elsuel.rowSize() || col >= elsuel.colSize() )
    {
        return false;
    }

    return ( elsuel[row][col] >= 0 );
}


//-----------------------------------------------------------------------
//   calcNormal()
//-----------------------------------------------------------------------


Vector<double> Mesh::calcNormal( int node1, int node2 ) const
{
    double dx = points[node2][0]- points[node1][0];
    double dy = points[node2][1]- points[node1][1];

    Vector<double> normal = {-dy, dx};

    return normal/normal.norm();
}


//-----------------------------------------------------------------------
//   checkDirection()
//-----------------------------------------------------------------------


bool Mesh::checkDirection( int idx1, int idx2, int node1, int node2 ) const
{            
    Vector<double> normal = calcNormal(node1, node2);

    return ( ( normal[0]*cents[idx1][0] + normal[1]*cents[idx1][1] ) > 
             ( normal[0]*cents[idx2][0] + normal[1]*cents[idx2][1] ) );    
}


//-----------------------------------------------------------------------
//   convertCellsToPoints()
//-----------------------------------------------------------------------


void Mesh::convertCellsToPoints
( 
    Vector<double>&       pointVals,
    const Vector<double>& cellVals 
) const

{            
   pointVals.reinit( points.rowSize() );

   for ( int i = 0; i < pointVals.rowSize(); i++ )
   {
        auto slice = std::make_tuple(elsup2[i], elsup2[i+1]);

        Vector<int> cellidxs = elsup1(slice);

        pointVals[i] = cellVals(cellidxs).mean();       
   }
}


//-----------------------------------------------------------------------
//   readMesh_()
//-----------------------------------------------------------------------


void Mesh::readMesh_( std::string& meshFile )
{
    std::ifstream file( meshFile );
    std::string line;
    std::string word;

    if (file.fail())
    {
        std::cout << "Import of mesh file " << meshFile << " failed\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    // check if mesh is of type unstructured
    if ( !checkMesh_( file ) )
    {
        return;
    }

    //find the number of points
    readPoints_( file );

    // read elems
    readElems_( file );
    
    // read cell types
    //TODO: make fool proof if it isn't given
    readElemTypes_( file );

    //read the physic types (boundary conditions)
    readPhysicTypes_( file );

    // group elems together into elements and boundary lines
    groupElems_();

    //finished
    file.close();
    return;
}


//-----------------------------------------------------------------------
//   checkMesh_()
//-----------------------------------------------------------------------


bool Mesh::checkMesh_(std::ifstream& file)
{

    std::string word;

    // Check if it is an unstructured grid
    while( word != "DATASET")
    {
        file >> word;
    }
    file >> word;

    if (word != UNSTRUCTURED_GRID)
    {
      std::cout << "Grid type " << word << " is not supported.\n";
      return false;
    }

    return true;
}


//-----------------------------------------------------------------------
//   readPoints_()
//-----------------------------------------------------------------------


void Mesh::readPoints_( std::ifstream& file )
{
    std::string word;
    std::string line;

    int numPoints;

    while( word != "POINTS")
    {
        file >> word;
    }
    file >> numPoints;
    std::getline(file,line);

    // Save the points to the points matrix
    points.reinit(numPoints, 3, 0.0);

    for ( int i = 0; i< numPoints; i++ )
    {
        // coordinates are written in 3 dimensions, z is always zero
        for ( int j = 0; j < 3; j++ ) 
        {   
            file >> points(i,j);
        }
    }
}


//-----------------------------------------------------------------------
//   readElems_()
//-----------------------------------------------------------------------


void Mesh::readElems_( std::ifstream& file )
{
    std::string word;
    std::string line;

    int numElems;

    while( word != "CELLS" )
    {
        file >> word;
    }
    file >> numElems;
    std::getline(file,line);

    // assume that elements will not have more entries than 9 (polygonal)
    elems.reinit(numElems, MAX_NODES +1, -1 );


    for ( int i = 0; i< numElems; i++ )
    {
        int length;
        file >> length;

        elems(i,0) = length;
        // coordinates are written in 3 dimensions z is always zero
        for ( int j = 1; j < length +1; j++ ) 
        {
            file >> elems(i,j);
        }
    }
}


//-----------------------------------------------------------------------
//   readElemTypes_()
//-----------------------------------------------------------------------


void Mesh::readElemTypes_( std::ifstream& file )
{
    std::string word;

    int numElems;

    while( word != "CELL_TYPES" )
    {
        file >> word;
    }
    file >> numElems;

    cellTypes.reinit(numElems);

    for ( int i = 0; i< numElems; i++ )
    {
        file >> cellTypes[i];
    }
}


//-----------------------------------------------------------------------
//   readPhysicTypes_()
//-----------------------------------------------------------------------


void Mesh::readPhysicTypes_( std::ifstream& file )
{
    std::string word;
    std::string line;

    int numElems;
    while( word != "CELL_DATA" )
    {
        file >> word;

        if (file.peek() == EOF)
        {
            return;
        }
    }
    
    file >> numElems;
    std::getline(file,line);
    std::getline(file,line);
    std::getline(file,line);

    physicTypes.reinit(numElems);

    for ( int i = 0; i< numElems; i++ )
    {
        file >> physicTypes[i];
    }
}


//-----------------------------------------------------------------------
//   groupElems_()
//-----------------------------------------------------------------------


void Mesh::groupElems_()
{
    int numElems = elems.rowSize();
    int numBcElems = 0;
    int numCells;

    for (int i =0; i < numElems; i++)
    {
        // 2 means 2 points, thus it is a line, i.e it is a boundary line
        if ( (elems)(i,0) == 2) 
        {
            numBcElems++;
        }
    }
    numCells = numElems - numBcElems;

    // assume that elements will not have more entries than 9 (polygonal)
    cells.reinit(numCells, MAX_NODES +1, -1);
    bcElems.reinit(numBcElems, 2, -1);
    bcKeys.reinit(numBcElems, -1);

    int idxElems = 0;
    int idxBc    = 0;
    for (int i = 0; i < numElems; i++)
    {
        if ( elems(i,0) == 2) 
        {
            bcElems(idxBc,0) = (elems)(i,1);
            bcElems(idxBc,1) = (elems)(i,2);
            bcKeys[idxBc]    = physicTypes[i];
            idxBc++;
            continue;
        }
        else
        {
            int length = (elems)(i,0);

            for (int j = 0; j<length + 1; j++)
            {
                cells(idxElems,j) = (elems)(i,j);
            }
            idxElems++;
        }
        //TODO: resize cells to appropriate size
    }
}


//-----------------------------------------------------------------------
//   calcCentroids_()
//-----------------------------------------------------------------------


void Mesh::calcCentroids_()
{

    cents.reinit(cells.rowSize(),2);

    for(int i = 0; i< cells.rowSize(); i++)
    {

        int length = cells(i,0);
        for (int j = 1; j < length+1; j++)
        {
            cents(i,0) += points( cells(i,j), 0 )/length;
            cents(i,1) += points( cells(i,j), 1 )/length;
        } 
    }
}


//-----------------------------------------------------------------------
//   calcSurfaces_()
//-----------------------------------------------------------------------



void Mesh::calcSurfaces_()
{
    surfaces.reinit(cells.rowSize());

    for (int i = 0; i < cells.rowSize(); i++)
    {
        double surf = 0.0;
        for( int j = 1; j < cells(i,0) +1; j++)
        {
        
            int node1 = cells(i,j);
            int node2;

            if ( j < cells(i,0) )
            {
                node2 = cells(i, j+1);
            }
            else
            {
                node2 = cells(i,1);
            }

            surf+= ( points(node1,0)*points(node2,1) 
                 -   points(node1,1)*points(node2,0) ) /2.0;
        }
        surfaces[i] = std::abs(surf);
    }
}


//-----------------------------------------------------------------------
//   calcVertElems_()
//-----------------------------------------------------------------------


void Mesh::calcVertElems_()
{

/*

Create vertices-elements connections, average number
of elements using 1 node will not go higher than 10.

The vertices-elements connections stored as a linked list, so they  are 
divided into 2 arrays, elsup1(elements surrounding points 1) stores the 
elements, and elsup2 stores the start and end index of the elements that
surround a point. example: for point ipoint, the elements surrounding
this point are stored in elsup1( elsup2(ipoint) to elsup2(ipoint+1)-1 );

*/

elsup1.reinit(points.rowSize()*MAX_NODES);
elsup2.reinit(points.rowSize()+1);
// Pass 1: First count the elements connected to each point

for (int i = 0; i < cells.rowSize(); i++)
{
    for (int j = 1; j < cells(i,0) + 1; j++)
    {

        int ipoint = cells(i,j) + 1;
        elsup2[ipoint] = elsup2[ipoint] +1;
    }
}

// Reshuffle the first pass
for (int i = 1; i < points.rowSize() + 1; i++)
{
    elsup2[i] = elsup2[i] + elsup2[i-1];
}

// Pass 2: store the elements in esup1
int maxLen = 0;
for (int i = 0; i < cells.rowSize(); i++)
{
    for (int j = 1; j < cells(i,0) + 1; j++)
    {
        int ipoint = cells(i,j);
        int istore = elsup2[ipoint] +1;

        elsup2[ipoint] = istore;
        elsup1[istore -1] = i;

        maxLen = std::max(maxLen,istore);
    }
}

// Reshuffle the second pass
for (int i = points.rowSize(); i > 0; i--)
{
    elsup2[i] = elsup2[i-1];
}
elsup2[0] = 0;

elsup1.resize(maxLen);

}


//-----------------------------------------------------------------------
//   calcVertBcs_()
//-----------------------------------------------------------------------


void Mesh::calcVertBcs_()
{


bcsup1.reinit(points.rowSize()*MAX_NODES);
bcsup2.reinit(points.rowSize()+1);

// Pass 1: First count the elements connected to each point

for (int i = 0; i < bcElems.rowSize(); i++)
{
    for (int j = 0; j <  bcElems.colSize(); j++)
    {

        int ipoint = bcElems(i,j) + 1;
        bcsup2[ipoint] = bcsup2[ipoint] +1;
    }
}


// Reshuffle the first pass
for (int i = 1; i < points.rowSize() + 1; i++)
{
    bcsup2[i] = bcsup2[i] + bcsup2[i-1];
}

// Pass 2: store the elements in esup1
int maxLen = 0;
for (int i = 0; i < bcElems.rowSize(); i++)
{
    for (int j = 0; j <  bcElems.colSize(); j++)
    {
        int ipoint = bcElems(i,j);
        int istore = bcsup2[ipoint] +1;

        bcsup2[ipoint] = istore;
        bcsup1[istore -1] = i;

        maxLen = std::max(maxLen,istore);
    }
}

// Reshuffle the second pass
for (int i = points.rowSize(); i > 0; i--)
{
    bcsup2[i] = bcsup2[i-1];
}
bcsup2[0] = 0;

bcsup1.resize(maxLen);

}


//-----------------------------------------------------------------------
//   calcFaceElems_()
//-----------------------------------------------------------------------


void Mesh::calcFaceElems_()
{

elsuel.reinit( cells.rowSize(),cells.colSize(), -1 );

for (int i = 0; i < cells.rowSize(); i++)
{
    //number of faces per cell is equal to the number of points per cell
    elsuel(i,0) = cells(i,0);

    for (int j = 1; j < cells(i,0) + 1; j++)
    {
        int node1 = cells(i,j);
        int node2;

        if ( j < cells(i,0) )
        {
            node2 = cells(i, j+1);
        }
        else
        {
            node2 = cells(i,1);
        }

        for( int istore = elsup2[node1]; istore< elsup2[node1+1]; istore++ )
        {
            int jelem = elsup1[istore];

            if (jelem != i)
            {
                int count = 0;

                for (int jnodes = 1; jnodes < cells(jelem,0)+1; jnodes++ )
                {
                    if ( cells(jelem,jnodes) == node1 ||  
                         cells(jelem,jnodes) == node2 )
                    {
                        count++;
                    }
                }

                if (count == 2) //a face always has only 2 nodes
                {
                    elsuel(i,j) = jelem;
                }
            }
        }


    }
}

}


//-----------------------------------------------------------------------
//   calcFaceBcs_()
//-----------------------------------------------------------------------


void Mesh::calcFaceBcs_()
{

bcsuel.reinit( cells.rowSize(),cells.colSize(), -1 );

for (int i = 0; i < cells.rowSize(); i++)
{
    //number of faces per cell is equal to the number of points per cell
    bcsuel(i,0) = cells(i,0);

    for (int j = 1; j < cells(i,0) + 1; j++)
    {
        int node1 = cells(i,j);
        int node2;

        if ( j < cells(i,0) )
        {
            node2 = cells(i, j+1);
        }
        else
        {
            node2 = cells(i,1);
        }

        for( int istore = bcsup2[node1]; istore< bcsup2[node1+1]; istore++ )
        {
            int jelem = bcsup1[istore];

            int count = 0;

            for (int jnodes = 0; jnodes < bcElems.colSize(); jnodes++ )
            {
                if ( bcElems(jelem,jnodes) == node1 ||  
                        bcElems(jelem,jnodes) == node2 )
                {
                    count++;
                }
            }

            if (count == 2) //a face always has only 2 nodes
            {
                bcsuel(i,j) = jelem;
            }
            
        }
    }
}

}