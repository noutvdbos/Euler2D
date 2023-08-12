
#include "Names.h"

//=============================================================================
//   class InputNames
//=============================================================================


//mesh
  const char* InputNames::MESH_FILE              = "meshFile";
  const char* InputNames::BOUNDARY_FILE          = "boundaryFile";
  const char* InputNames::VTK_TYPE               = "vtk";
  const char* InputNames::CSV_TYPE               = "csv";

  //time
  const char* InputNames::START_TIME             = "startTime";
  const char* InputNames::END_TIME               = "endTime";
  const char* InputNames::TIME_STEP              = "timeStep";
  const char* InputNames::ADJUST_TIME_STEP       = "adjustTimeStep";
  const char* InputNames::MAX_COURANT            = "maxCourant";

  //initial conditions
  const char* InputNames::INIT_METHOD            = "initMethod";
  const char* InputNames::INIT_FILE              = "initFile";
  const char* InputNames::DENSITY                = "rho";
  const char* InputNames::PRESSURE               = "p";
  const char* InputNames::U_VELOCITY             = "u";
  const char* InputNames::V_VELOCITY             = "v";
  const char* InputNames::VELOCITY               = "U";
  const char* InputNames::SPECIFIC_HEAT_RATIO    = "specificHeatRatio";

  //reference conditions
  const char* InputNames::REF_DENSITY            = "rhoref";
  const char* InputNames::REF_PRESSURE           = "pref";
  const char* InputNames::REF_VELOCITY           = "Vref";
  const char* InputNames::REF_LENGTH             = "lref";

  //solver
  const char* InputNames::TIME_SCHEME            = "timeScheme";
  const char* InputNames::FLUX_SCHEME            = "fluxScheme";

  //results
  const char* InputNames::RESULTS_FOLDER         = "resultsFolder";
  const char* InputNames::RESULTS_PREFIX         = "resultsPrefix"; 
  const char* InputNames::WRITE_METHOD           = "writeMethod";
  const char* InputNames::WRITE_INTERVAL         = "writeInterval";
  
  const char* InputNames::WRITE_COEFFICIENTS     = "writeCoefficients";
  const char* InputNames::SELECTED_BOUNDARIES    = "selectedBoundaries";


//=============================================================================
//   class ProgramNames
//=============================================================================


const char* ProgramNames::PROGRAM                = "Euler2D";


//=============================================================================
//   class InitMethodNames
//=============================================================================


const char* InitMethodNames::UNIFORM             = "uniform";
const char* InitMethodNames::FILE                = "file";


//=============================================================================
//   class FluxSchemeNames
//=============================================================================


const char* FluxSchemeNames::AUSM_PLUS           = "AUSM+";
const char* FluxSchemeNames::UPWIND              = "upwind";


//=============================================================================
//   class TimeSchemeNames
//=============================================================================


const char* TimeSchemeNames::EXPLICIT_EULER      = "explicitEuler";


//=============================================================================
//   class WriteMethodNames
//=============================================================================


const char* WriteMethodNames::TIME               = "time";
const char* WriteMethodNames::ITERATION          = "iteration";


//=============================================================================
//   class BoundaryVarNames
//=============================================================================


const char* BoundaryVarNames::KEY                = "key";
const char* BoundaryVarNames::TYPE               = "type";
const char* BoundaryVarNames::U_VELOCITY         = "u";
const char* BoundaryVarNames::V_VELOCITY         = "v";
const char* BoundaryVarNames::PRESSURE           = "p";
const char* BoundaryVarNames::DENSITY            = "rho";


//=============================================================================
//   class BoundaryTypeNames
//=============================================================================


const char* BoundaryTypeNames::SUPERSONIC_INLET   = "supersonicInlet";
const char* BoundaryTypeNames::SUPERSONIC_OUTLET  = "supersonicOutlet";
const char* BoundaryTypeNames::SUBSONIC_INLET     = "subsonicInlet";
const char* BoundaryTypeNames::SUBSONIC_OUTLET    = "subsonicOutlet";
const char* BoundaryTypeNames::FREESTREAM         = "freestream";
const char* BoundaryTypeNames::WALL               = "wall";

