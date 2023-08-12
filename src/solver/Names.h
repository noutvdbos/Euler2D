#pragma once


//-----------------------------------------------------------------------------
//   class InputNames
//-----------------------------------------------------------------------------


class InputNames
{

public:

  // static variables

  //mesh
  static const char* MESH_FILE;
  static const char* BOUNDARY_FILE;
  static const char* VTK_TYPE;
  static const char* CSV_TYPE;

  //time
  static const char* START_TIME;
  static const char* END_TIME;
  static const char* TIME_STEP;
  static const char* ADJUST_TIME_STEP;
  static const char* MAX_COURANT;


  //initial conditions
  static const char* INIT_METHOD;
  static const char* INIT_METHOD_UNIFORM;
  static const char* INIT_METHOD_FILE;
  static const char* INIT_FILE;
  static const char* DENSITY;
  static const char* PRESSURE;
  static const char* U_VELOCITY;
  static const char* V_VELOCITY;
  static const char* VELOCITY;
  static const char* SPECIFIC_HEAT_RATIO;


  //reference conditions
  static const char* REF_DENSITY;
  static const char* REF_PRESSURE;
  static const char* REF_VELOCITY;
  static const char* REF_LENGTH;


  //solver
  static const char* TIME_SCHEME;
  static const char* FLUX_SCHEME;

  //results
  static const char* RESULTS_FOLDER;
  static const char* RESULTS_PREFIX;
  static const char* WRITE_METHOD;
  static const char* WRITE_INTERVAL;
  
  static const char* WRITE_COEFFICIENTS;
  static const char* SELECTED_BOUNDARIES;

};


//-----------------------------------------------------------------------------
//   class ProgramNames
//-----------------------------------------------------------------------------


class ProgramNames
{

public:

  // static variables

  static const char* PROGRAM;

};


//-----------------------------------------------------------------------------
//   class InitMethodNames
//-----------------------------------------------------------------------------


class InitMethodNames
{

public:

  // static variables

  static const char* UNIFORM;
  static const char* FILE;

};


//-----------------------------------------------------------------------------
//   class FluxSchemeNames
//-----------------------------------------------------------------------------


class FluxSchemeNames
{

public:

  // static variables

  static const char* AUSM_PLUS;
  static const char* UPWIND;

};


//-----------------------------------------------------------------------------
//   class TimeSchemeNames
//-----------------------------------------------------------------------------


class TimeSchemeNames
{

public:

  // static variables

  static const char* EXPLICIT_EULER;

};


//-----------------------------------------------------------------------------
//   class WriteMethodNames
//-----------------------------------------------------------------------------


class WriteMethodNames
{

public:

  // static variables

  static const char* TIME;
  static const char* ITERATION;

};


//-----------------------------------------------------------------------------
//   class BoundaryVarNames
//-----------------------------------------------------------------------------


class BoundaryVarNames
{

public:

  // static variables

  static const char* KEY;
  static const char* TYPE;
  static const char* U_VELOCITY;
  static const char* V_VELOCITY;
  static const char* PRESSURE;
  static const char* DENSITY;

};


//-----------------------------------------------------------------------------
//   class BoundaryTypeNames
//-----------------------------------------------------------------------------

class BoundaryTypeNames
{

public:

  // static variables

  static const char* SUPERSONIC_INLET;
  static const char* SUPERSONIC_OUTLET;
  static const char* SUBSONIC_INLET;
  static const char* SUBSONIC_OUTLET;
  static const char* FREESTREAM;
  static const char* WALL;
};