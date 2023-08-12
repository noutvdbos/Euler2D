#pragma once

#include <string>
#include <filesystem>
#include <fstream>
#include "utils/Matrix.h"
#include "Config.h"
#include "Mesh.h"

//-----------------------------------------------------------------------
//   class Writer
//-----------------------------------------------------------------------


class Writer
{   

public:
  // Constructor & destructor
  Writer
  (
      Config config,
      Mesh*  mesh      
  );

  ~Writer();

  //Member functions:

  void write(const Matrix<double>& results, std::string suffix ) const;

  void writeWallCp(const Matrix<double>& results, std::string suffix ) const;

  void writeWallCfs( double cfx, double cfy, double time, int iter, std::string boundary ) const;

  std::string getTimeSuffix( double time, double endTime ) const;

private:
  
  //member functions:
  void writeBaseFile_();

  void writeScalar_
  ( 
      const Vector<double>& vec,
      std::ofstream& file,
      const char* name ) const;
    
  void writeVector3D_
  ( 
      const Vector<double>& vec1,
      const Vector<double>& vec2,
      const Vector<double>& vec3,
      std::ofstream& file,
      const char* name ) const;    


  //member variables:
  std::string meshFile_;
  std::string baseFile_;
  std::string resultsPrefix_;
  std::string resultsFolder_;

  Mesh* mesh_;

};