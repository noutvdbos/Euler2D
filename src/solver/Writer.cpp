
#include <filesystem>
#include <fstream>

#include "Writer.h"
#include "Names.h"

//=======================================================================
//   class Writer
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


Writer::Writer
(
    Config config,
    Mesh*  mesh    
)
{
    meshFile_      = config.getMeshFile();
    resultsPrefix_ = config.getResultsPrefix();
    resultsFolder_ = config.getResultsFolder();

    mesh_ = mesh;

    writeBaseFile_();
}


Writer::~Writer()
{}


//-----------------------------------------------------------------------
//   write
//-----------------------------------------------------------------------

void Writer::write
(
    const Matrix<double>& results,
    std::string           suffix
) const

{
    using std::ofstream;

    std::string filePath = resultsFolder_ + "/" + resultsPrefix_ + "_" + suffix +
                           "." + InputNames::VTK_TYPE;
    
    if ( !std::filesystem::exists(resultsFolder_) )
    {
        std::filesystem::create_directory(resultsFolder_);
    }

    if ( std::filesystem::exists( filePath ) )
    {
        std::remove(filePath.c_str());
    }

    if (!std::filesystem::copy_file(baseFile_, filePath) )
    {
        std::cout << "Copying of base file " << baseFile_ << " failed\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    // The writer writes both point and cell data, so first the results must 
    // be converted into points.

    // The z-component of velocity is zero, but needs to be included for vtk.
    Vector<double> rho;
    Vector<double> u;
    Vector<double> v;
    Vector<double> w;  
    Vector<double> p;
    Vector<double> cp;
    Vector<double> mach;

    Vector<double> wCells;

    mesh_->convertCellsToPoints(rho, results[0]);
    mesh_->convertCellsToPoints(u, results[1]);
    mesh_->convertCellsToPoints(v, results[2]);
    mesh_->convertCellsToPoints(p, results[3]);
    mesh_->convertCellsToPoints(cp, results[4]);
    mesh_->convertCellsToPoints(mach, results[5]);

    // Set velocity in z-direction equal to zero.
    w =      Vector<double> ( u.rowSize(), 0.0 );
    wCells = Vector<double> (results[0].rowSize(), 0.0 );

    ofstream file;

    file.open( filePath, ofstream::app|ofstream::out );

    file << "\n" << "POINT_DATA " << rho.rowSize() << "\n";

    writeScalar_(rho, file, "rho");
    writeScalar_(p, file, "p");
    writeScalar_(cp, file, "Cp");
    writeScalar_(mach, file, "Mach");
    
    writeVector3D_(u,v,w,file, "U");

    
    file << "\n" << "CELL_DATA " <<results[0].rowSize() << "\n";
    
    writeScalar_(results[0], file, "rho");
    writeScalar_(results[3], file, "p");
    writeScalar_(results[4], file, "Cp");
    writeScalar_(results[5], file, "Mach");
    
    // Set velocity in z-direction equal to zero again.
    writeVector3D_(results[1],results[2],wCells,file, "U");

    file.close();
}


//-----------------------------------------------------------------------
//   writeWallCp()
//-----------------------------------------------------------------------


void Writer::writeWallCp
(
    const Matrix<double>& results,
    std::string           suffix
) const
{
    //Wall cp is written in csv for convenience

    using std::ofstream;

    std::string filePath = resultsFolder_ + "/" + resultsPrefix_ + "_" + suffix +
                           "." + InputNames::CSV_TYPE;
    
    if ( !std::filesystem::exists(resultsFolder_) )
    {
        std::filesystem::create_directory(resultsFolder_);
    }

    if ( std::filesystem::exists( filePath ) )
    {
        std::remove(filePath.c_str());
    }

    ofstream file;
    file.open( filePath, ofstream::out );

    file << "X" << " Y" << " Z" << " Cp\n";

    for (int i = 0; i < results.rowSize(); i++)
    {
        file << results[i][0] << " " << results[i][1]  << " " 
             << results[i][2] << " " << results[i][3]  << "\n" ;
    }

    file.close();
}


//-----------------------------------------------------------------------
//   writeWallCfs()
//-----------------------------------------------------------------------


void Writer::writeWallCfs
(
   double       cfx, 
   double       cfy, 
   double       time, 
   int          iter, 
   std::string  boundary
) const
{
    //Wall cfs is written in csv for convenience

    using std::ofstream;

    std::string filePath = resultsFolder_ + "/cf_" + boundary +
                           "." + InputNames::CSV_TYPE;
    
    if ( !std::filesystem::exists(resultsFolder_) )
    {
        std::filesystem::create_directory(resultsFolder_);
    }

    if ( !std::filesystem::exists( filePath ) )
    {
       ofstream file;
       file.open( filePath, ofstream::out );

       file << "time" << " iteration" << " Cfx" << " Cfy\n";

       file.close();
    }

    ofstream file;
    file.open( filePath, ofstream::app|ofstream::out );

    file << time << " " << iter << " " << cfx << " " << cfy <<"\n";

    file.close();
}


//-----------------------------------------------------------------------
//   writeBaseFile_
//-----------------------------------------------------------------------


void Writer::writeBaseFile_()
{
    using std::ofstream;

    baseFile_ = resultsFolder_ + "/" + "base_" + resultsPrefix_ 
                               + "." + InputNames::VTK_TYPE;
    
    if ( !std::filesystem::exists(resultsFolder_) )
    {
        std::filesystem::create_directory(resultsFolder_);
    }

    if ( std::filesystem::exists( baseFile_ ) )
    {
        std::remove(baseFile_.c_str());
    }

    ofstream file;
    file.open( baseFile_, ofstream::out );

    file << "# vtk DataFile Version 2.0\n"
         << "mesh, Created by Gmsh\n"
         << "ASCII\n"
         << "DATASET UNSTRUCTURED_GRID\n";

    //write point locations
    file << "POINTS " << mesh_->points.rowSize() << " double\n";

    for( int i = 0; i < mesh_->points.rowSize(); i++)
    {
        for (int j = 0; j < mesh_->points.colSize(); j++)
        {
            file << mesh_->points(i,j) << " ";
        }
        file << "\n";
    }
    

    int inputSum = 0;

    // calculate total sum of inputs for the cells.
    for (int i = 0; i < mesh_->cells.rowSize(); i++)
    {
        inputSum += mesh_->cells(i,0) + 1;
    }

    //write cell matrix
    file << "\n" << "CELLS " << mesh_->cells.rowSize() << " " << inputSum  << "\n";
    
    for( int i = 0; i < mesh_->cells.rowSize(); i++)
    {
        for (int j = 0; j < mesh_->cells(i,0) +1; j++)
        {
            file << mesh_->cells(i,j) << " ";
        }
        file << "\n";
    }
    
    // write cell types
    // only the inner cells are written, which are all of unstructured type
    // This type has the code 9 in vtk, so store 9 for each cell.

    file << "\n" << "CELL_TYPES " << mesh_->cells.rowSize() << "\n";
    
    for( int i = 0; i < mesh_->cells.rowSize(); i++)
    {
        file << 9 << "\n";
    }

    file.close();
}


//-----------------------------------------------------------------------
//   writeScalar_
//-----------------------------------------------------------------------


void Writer::writeScalar_
(
    const Vector<double>& vec,
    std::ofstream&        file,
    const char*           name
) const
{
    file << "SCALARS " << name << " double\n"
         << "LOOKUP_TABLE default\n"; 

    for (int i = 0; i < vec.rowSize(); i++)
    {
        file << vec[i] << "\n";
    }

    file << "\n";

}


//-----------------------------------------------------------------------
//   writeVector3D_
//-----------------------------------------------------------------------


void Writer::writeVector3D_
(
    const Vector<double>& vec1,
    const Vector<double>& vec2,
    const Vector<double>& vec3,
    std::ofstream&        file,
    const char*           name
) const
{
    assert(vec1.rowSize() == vec2.rowSize() && vec2.rowSize() == vec3.rowSize() );

    file << "VECTORS " << name << " double\n";

    for (int i = 0; i < vec1.rowSize(); i++)
    {
        file << vec1[i] << " " << vec2[i] << " " << vec3[i] << "\n" ;
    }

    file << "\n";

}
