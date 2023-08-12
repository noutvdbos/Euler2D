#pragma once

#include <fstream>
#include "utils/Matrix.h"
#include "utils/Vector.h"
#include "Mesh.h"
#include "Flux.h"
#include "Face.h"
#include "Config.h"

//-----------------------------------------------------------------------
//   class Solver
//-----------------------------------------------------------------------


class Solver
{   
public:

  //Constructor
  Solver(Config config );

  //Destructor
  ~Solver                 ();


private:

  void      initCellStates_();
  void      initFaces_();
  void      calcTimeStep_();
  void      updateCellStates_();


public:

  Mesh*     getMesh() const { return mesh_; }

  double    calcMaxCourantNr() const;
  void      iterate();
  double    getTime() const;
  double    getTimeStep() const;
  int       getIteration() const;

  bool      hasFinished() const;
  void      getResults ( Matrix<double>& results ) const;

  void      getWallCoeffs 
            ( 
              Matrix<double>& results, 
              double& cfx,
              double& cfy,
              std::string boundary ) const;


private:
  Config            config_;
  double            time_;
  double            dt_;

  double            writeInt_;
  std::string       writeMethod_;

  int               iter_;
  int               wallCount_;

  Mesh*             mesh_;
  Matrix<double>    cellStates_;
  Matrix<Face>      faces_;
  double            gamma_;

  double            refP_;
  double            refRho_;
  double            refVmag_;
  double            refLength_;
  double            refMach_;
  double            refStagP_;
  double            refDynP_;



};

