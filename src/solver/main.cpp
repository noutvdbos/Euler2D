#include <iostream>
#include <list>
#include <chrono>
#include <string>
#include <filesystem>

#include "Config.h"
#include "Solver.h"
#include "Writer.h"
#include "Mesh.h"
#include "Names.h"
#include "utils/Constants.h"
#include "utils/Matrix.h"
#include "utils/Vector.h"
#include "utils/utils.h"

//-----------------------------------------------------------------------
//   main()
//-----------------------------------------------------------------------


int main ( int argc, char** argv )
{
  using std::chrono::high_resolution_clock;
  using std::chrono::duration;  

  // declaration of helper function:
  void writeCoeffs
  ( 
    Matrix<double>& results,
    const Config& config,
    const Solver& solver,
    const Writer& writer,
    std::string suffix,
    double time,
    int iter 
  );

  // start the program by printing out some nice ASCII art
  std::cout <<  
  " .----------------.  .----------------.  .----------------.  .----------------.  .----------------.   .----------------.  .----------------. \n" << 
  "| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. | | .--------------. || .--------------. |\n" << 
  "| |  _________   | || | _____  _____ | || |   _____      | || |  _________   | || |  _______     | | | |    _____     | || |  ________    | |\n" << 
  "| | |_   ___  |  | || ||_   _||_   _|| || |  |_   _|     | || | |_   ___  |  | || | |_   __ \\    | | | |   / ___ `.   | || | |_   ___ `.  | |\n" << 
  "| |   | |_  \\_|  | || |  | |    | |  | || |    | |       | || |   | |_  \\_|  | || |   | |__) |   | | | |  |_/___) |   | || |   | |   `. \\ | |\n" << 
  "| |   |  _|  _   | || |  | '    ' |  | || |    | |   _   | || |   |  _|  _   | || |   |  __ /    | | | |   .'____.'   | || |   | |    | | | |\n" << 
  "| |  _| |___/ |  | || |   \\ `--' /   | || |   _| |__/ |  | || |  _| |___/ |  | || |  _| |  \\ \\_  | | | |  / /____     | || |  _| |___.' / | |\n" << 
  "| | |_________|  | || |    `.__.'    | || |  |________|  | || | |_________|  | || | |____| |___| | | | |  |_______|   | || | |________.'  | |\n" << 
  "| |              | || |              | || |              | || |              | || |              | | | |              | || |              | |\n" << 
  "| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' | | '--------------' || '--------------' |\n" << 
  " '----------------'  '----------------'  '----------------'  '----------------'  '----------------'   '----------------'  '----------------' \n\n\n";

  std::cout << "Welcome to Euler2D! Your simulation is starting now.\n\n";
  std::cout << std::endl;
  
  // Initialise timer variables
  auto t1 = high_resolution_clock::now();
  auto t2 = t1;
  auto tSolve = t2;
  duration<double, std::milli> time;
  
  std::string inputFile;
  std::string path;
  size_t loc;

  
  //First print absolute path
  auto absPath = std::filesystem::current_path();

  std::cout << "Your simulation is run from the location: " << absPath 
            << ". All shown file locations are relative to this path.\n\n";

  if (argc > 1)
  {
    inputFile = argv[1];
  }
  else
  {
    std::cout << "No input file is given, assuming the name "
              << "inputFile.euler\n\n";

    inputFile = "inputFile.euler";
  }


  // Get the file path of the input file:
  char ch = '/';
  loc = inputFile.find_last_of(ch);

  if (loc != std::string::npos)
  {
    path = inputFile.substr(0,loc+1);
  }

  Config config ( path, inputFile );
  Solver solver ( config );
  Writer writer ( config, solver.getMesh() );

  Matrix<double> results;
  Matrix<double> wallResults;
  double interval = config.getWriteInterval();

  //write the initial results;
  solver.getResults(results);
  writer.write(results, std::to_string(solver.getIteration()) );

  writeCoeffs( wallResults, config, solver, writer, 
               std::to_string(solver.getIteration()),
               solver.getTime(), solver.getIteration() );
  
  //Start solving
  
  std::cout << 
    "\n--------------------------SOLVER----------------------------\n";

  tSolve = high_resolution_clock::now();

  while ( !solver.hasFinished() )
  {
    //update the solution
    solver.iterate();

    // write the solution
    if 
    ( (config.getWriteMethod() == WriteMethodNames::TIME && 
        ( std::fmod( solver.getTime(), interval ) < Constants::VERY_SMALL || 
            std::abs(std::fmod( solver.getTime(), interval )-interval)
            < Constants::VERY_SMALL) )

        ||

      (config.getWriteMethod() == WriteMethodNames::ITERATION && 
        ( std::fmod( solver.getIteration(), interval ) < Constants::VERY_SMALL || 
            std::abs(std::fmod( solver.getIteration(), interval )-interval)
            < Constants::VERY_SMALL) )
    )
    {

      std::string suffix = std::to_string( solver.getIteration() );

      solver.getResults( results );
      writer.write( results, suffix );

      writeCoeffs( wallResults, config, solver, writer, suffix, 
                   solver.getTime(), solver.getIteration() );

      // give information about the time and cfl etc.
      
      t2 = high_resolution_clock::now();
      time = t2 - tSolve;

      std::cout << "---At time = " << solver.getTime() << "s of " << 
                   config.getEndTime() << "s--- \n";
      std::cout << "-Maximum Courant Number: " << solver.calcMaxCourantNr() << "\n";
      std::cout << "-Current Time Step: " << solver.getTimeStep() << "s. \n";
      std::cout << "-Wall time: " << time.count()/1000 << "s. \n\n";
    }
  }

  //write final results
  solver.getResults( results );
  writer.write( results, "_final" );
  
  writeCoeffs( wallResults, config, solver, writer, "_final", 
               solver.getTime(), solver.getIteration() );

  std::cout<<std::endl;

  t2 = high_resolution_clock::now();
  time = t2 - tSolve;

  std::cout << "\n---End of simulation---\n";
  std::cout << "-Total run time = " << time.count()/1000 << "s.\n";   

  time = t2-t1;

  std::cout << "-Total time = " << time.count()/1000 << "s.\n";  
}


// helper functions

void writeCoeffs
( 
  Matrix<double>& results,
  const Config& config,
  const Solver& solver,
  const Writer& writer,
  std::string suffix,
  double time,
  int iter 
)
{

  if ( config.getWriteCoefficients() )
  {
    for (int i = 0; i < config.getSelBoundarySize(); i++)
    {
      double cfx, cfy;
      std::string cpSuffix = "cp_" + config.getSelBoundary(i) + "_" + suffix;

      solver.getWallCoeffs(results, cfx, cfy, config.getSelBoundary(i));

      writer.writeWallCp(results, cpSuffix);
      writer.writeWallCfs( cfx, cfy, time, iter, config.getSelBoundary(i) );
    }
  }
}