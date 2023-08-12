#include <iostream>
#include <string>
#include <fstream>
#include <cctype>
#include <algorithm>

#include "Config.h"
#include "Names.h"
#include "utils/Constants.h"


//=======================================================================
//   class Config
//=======================================================================

//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------


Config::Config()
{   
    path_           = ""; 
    meshFile_       = "mesh.vtk";
    boundaryFile_   = "boundaries.dat";

    startTime_      = 0;
    endTime_        = 1;
    dt_             = 1;
    adjustTimeStep_ = false;
    maxCourant_     = 0.5;

    initMethod_     = InitMethodNames::UNIFORM;
    initFile_       = "init.vtk";

    initRho_        = 0;
    initP_          = 0;
    initU_          = 0;
    initV_          = 0;
    gamma_          = 0;

    refRho_         = 0;
    refP_           = 0;
    refVmag_        = 0;
    refLength_      = 0;

    fluxScheme_     = FluxSchemeNames::AUSM_PLUS;
    timeScheme_     = TimeSchemeNames::EXPLICIT_EULER;

    resultsFolder_  = "results"; 		
    resultsPrefix_  = "results"; 		

    writeMethod_    = WriteMethodNames::ITERATION;

    writeInterval_  = 1; 	

    writeCoefficients_ = false;

    boundaries_.reinit(1);
}

Config::Config 
(
    const std::string&      path,
    const std::string&      fileName
)
{
    Config();
    path_ = path;

    std::ifstream file(fileName);

    std::cout << 
    "\n-----------------------CONFIGURATION------------------------\n";

    //NOTE: It is important that the config file is read first.
    //This is because the boundary file reading relies on knowing 
    //the path to the boundary file

    

    std::cout << "Using input file: " << fileName << "\n\n";

    std::cout << "---Reading input file---\n";

    if (file.fail())
    {
        std::cout << "---Reading input file failed---\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    readConfigFile_(file);
    
    std::cout << "-Read input file succesfully. \n\n";

    std::cout << "---Reading boundary file---\n";
    readBoundaryFile_();
    
    std::cout << "-Read boundary file succesfully. \n\n";



    file.close();

}

Config::~Config()
{}


//-----------------------------------------------------------------------
//  getInitVariables()
//-----------------------------------------------------------------------


void Config::getInitVariables
(
    Vector<double>& p,
    Vector<double>& rho,
    Vector<double>& u,
    Vector<double>& v,
    int         rowSize
)

{
    std::ifstream file ( initFile_);
    
    if (file.fail())
    {
        std::cout << "ERROR: Could not open file: " << initFile_ << " \n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    std::string word;
    std::string line;

    int numCells;

    while( word != "CELL_DATA")
    {
        file >> word;
    }
    file >> numCells;


    //Check if the numCells is equal for the mesh and the initial conditions
    assert( numCells == rowSize );

    p.reinit(numCells);
    rho.reinit(numCells);
    u.reinit(numCells);
    v.reinit(numCells);

    int varCount = 0;
    while (file.peek() != EOF)
    {
        file >> word;

        //NOTE that capitals matter for vtk file info, as U is for the 
        if ( word == "SCALARS" )
        {
            file >> word;

            if ( word == InputNames::DENSITY )
            {
                varCount++;
                std::getline(file,line);
                std::getline(file,line);

                for ( int i = 0; i< numCells; i++ )
                {  
                    file >> rho[i];
                }
            }

            if ( word == InputNames::PRESSURE )
            {
                varCount++;
                std::getline(file,line);
                std::getline(file,line);

                for ( int i = 0; i< numCells; i++ )
                {  
                    file >> p[i];
                }
            }
        }

        if ( word == "VECTORS" )
        {
            varCount+=2;

            std::string temp;
            std::getline(file,line);

            // the third dimension is not used, so store this in temp.
            for ( int i = 0; i< numCells; i++ )
            {  
                file >> u[i] >> v[i] >> temp;
            }
        }
        
    }

    if ( varCount < 4 )
    {
        std::cout << "ERROR: Not all variables (u,v,p,rho) are in the initialisation file.\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

}


//-----------------------------------------------------------------------
//   readConfigFile_()
//-----------------------------------------------------------------------


void Config::readConfigFile_(std::ifstream& file)
{
    int  varCount    = 0;
    int  refVarCount = 0;
    bool declaredInitFile = false;

    while (file.peek() != EOF)
    {
        std::string word;
        std::string line;

        file >> word;


        // If it is a comment, skip the line or the word. Comments are the same
        // as in c++, using //comment or /* comment spanning multiple lines*/

        if ( word.find("//") == 0  )
        {
            std::getline(file,line);
        }

        if ( word.find("/*") == 0 )
        {
            while ( word.find("*/") == word.npos  && file.peek() != EOF )
            {
                file >> word;
            }

            file >> word;
        }

        // Read the input
        
        if ( toLower_(word) == toLower_(InputNames::MESH_FILE) )
        {
            file >> meshFile_;
            
            meshFile_ = path_ + meshFile_;

            std::string fileType =  meshFile_.substr(meshFile_.find_last_of('.') + 1);

            if (toLower_(fileType) != toLower_( InputNames::VTK_TYPE ) )
            {
                std::cout << "ERROR: Mesh format '." << fileType << 
                             "' is not supported. Exiting " << 
                             ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }

        }

        if ( toLower_(word) == toLower_(InputNames::BOUNDARY_FILE) )
        {
            file >> boundaryFile_;

            boundaryFile_ = path_ + boundaryFile_;
        }


        if ( toLower_(word) == toLower_(InputNames::START_TIME) )
        {
            file >> startTime_;

            //TODO: check if the input is indeed a float
        }

        if ( toLower_(word) == toLower_(InputNames::TIME_STEP) )
        {
            file >> dt_;

            //TODO: check if the input is indeed a float
        }

        if ( toLower_(word) == toLower_(InputNames::END_TIME) )
        {
            file >> endTime_;

            //TODO: check if the input is indeed a float
        }

        if ( toLower_(word) == toLower_(InputNames::ADJUST_TIME_STEP) )
        {
            std::string temp;

            file >> temp;

            if (toLower_(temp) == "true")
            {
                adjustTimeStep_ = true;
            }
            else if (toLower_(temp) == "false")
            {
                adjustTimeStep_ = false;
            }
            else
            {
                std::cout << "ERROR: " << temp << 
                             " is not an option for AdjustTimeStep. \n" <<
                              "Options are: (true/false). Exiting " << 
                             ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }

        }

        if ( toLower_(word) == toLower_(InputNames::MAX_COURANT) )
        {
            file >> maxCourant_;

            //TODO: check if the input is indeed a float
        }


        // INITIAL CONDITIONS

        if ( toLower_(word) == toLower_(InputNames::INIT_METHOD) )
        {
            std::string temp;

            file >> initMethod_;

            if ( ( toLower_(initMethod_) != InitMethodNames::UNIFORM ) && 
                 ( toLower_(initMethod_) != InitMethodNames::FILE ) )
            {
                 std::cout << "ERROR: " << initMethod_ << 
                             " is not an option for initMethod. \n" <<
                               "Exiting " << 
                             ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }
        }

        if ( toLower_(word) == toLower_(InputNames::INIT_FILE) )
        {
            declaredInitFile = true;

            file >> initFile_;

            initFile_ = path_ + initFile_;
        }

        if ( toLower_(word) == toLower_(InputNames::DENSITY) )
        {
            varCount++;
            if ( toLower_(initMethod_) == InitMethodNames::UNIFORM )
            {
                file >> initRho_;
            }
        }

        if ( toLower_(word) == toLower_(InputNames::PRESSURE) )
        {
            varCount++;
            if ( toLower_(initMethod_) == InitMethodNames::UNIFORM )
            {
                file >> initP_;
            }
        }

        if ( toLower_(word) == toLower_(InputNames::U_VELOCITY) )
        {
            varCount++;
            if ( toLower_(initMethod_) == InitMethodNames::UNIFORM )
            {
                file >> initU_;
            }
        }

        if ( toLower_(word) == toLower_(InputNames::V_VELOCITY) )
        {
            varCount++;
            if ( toLower_(initMethod_) == InitMethodNames::UNIFORM )
            {
                file >> initV_;
            }
        }

        if ( toLower_(word) == toLower_(InputNames::SPECIFIC_HEAT_RATIO) )
        {
            file >> gamma_;
        }
 

        // REFERENCE CONDITIONS

        if ( toLower_(word) == toLower_(InputNames::REF_DENSITY) )
        {
            refVarCount++;
            
            file >> refRho_;
            
        }

        if ( toLower_(word) == toLower_(InputNames::REF_PRESSURE) )
        {
            refVarCount++;
            
            file >> refP_;
            
        }

        if ( toLower_(word) == toLower_(InputNames::REF_VELOCITY) )
        {
            refVarCount++;
            
            file >> refVmag_;
            
        }        

        if ( toLower_(word) == toLower_(InputNames::REF_LENGTH) )
        {
            refVarCount++;
            
            file >> refLength_;
            
        } 


        // SOLVER

        if ( toLower_(word) == toLower_(InputNames::TIME_SCHEME) )
        {
            file >> timeScheme_;

            if ( toLower_(timeScheme_) != toLower_(TimeSchemeNames::EXPLICIT_EULER ) )
            {
                std::cout << "ERROR: " << timeScheme_ << 
                                " is not an option for timeScheme. \n" <<
                                "Exiting " << 
                                ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }
        }

        if ( toLower_(word) == toLower_(InputNames::FLUX_SCHEME) )
        {
            file >> fluxScheme_;

            if ( toLower_(fluxScheme_) != toLower_( FluxSchemeNames::AUSM_PLUS ) )
            {
                std::cout << "ERROR: " << fluxScheme_ << 
                                " is not an option for fluxScheme. \n" <<
                                "Exiting " << 
                                ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }
        }


        // RESULTS

        if ( toLower_(word) == toLower_(InputNames::RESULTS_FOLDER) )
        {
            file >> resultsFolder_;

            resultsFolder_ = path_ + resultsFolder_;
        }

        if ( toLower_(word) == toLower_(InputNames::RESULTS_PREFIX) )
        {
            file >> resultsPrefix_;
        }

        if( toLower_(word) == toLower_(InputNames::WRITE_METHOD) )
        {
            file >> writeMethod_;

             if ( toLower_(writeMethod_) != toLower_( WriteMethodNames::TIME ) && 
                  toLower_(writeMethod_) != toLower_( WriteMethodNames::ITERATION ) )
            {
                std::cout << "ERROR: " << writeMethod_ << 
                                " is not an option for writeMethod. \n" <<
                                "Exiting " << 
                                ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }
        }

        if ( toLower_(word) == toLower_(InputNames::WRITE_INTERVAL) )
        {
            file >> writeInterval_;
        }

        if ( toLower_(word) == toLower_(InputNames::WRITE_COEFFICIENTS) )
        {
            std::string temp;

            file >> temp;

            if (toLower_(temp) == "true")
            {
                writeCoefficients_ = true;
            }
            else if (toLower_(temp) == "false")
            {
                writeCoefficients_ = false;
            }
            else
            {
                std::cout << "ERROR: " << temp << 
                             " is not an option for writeCoefficients. \n" <<
                              "Options are: (true/false). Exiting " << 
                             ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }

        }

        if ( toLower_(word) == toLower_(InputNames::SELECTED_BOUNDARIES) )
        {
            std::getline(file,line);
            
            readSelectedBCs_(line);
        }


    }

    if ( varCount < 4 && 
         toLower_(initMethod_) == toLower_(InitMethodNames::UNIFORM) )
    {
        std::cout << "ERROR: Not all variables (u,v,p,rho) are in the input file.\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    if ( refVarCount < 4 && writeCoefficients_ )
    {
        std::cout << "ERROR: Not all reference variables are in the initialisation file. " 
                  << "They are needed for calculating force and pressure coefficients."
                  << "Needed are: pref, rhoref, Vref and lref.\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    if ( gamma_ <= Constants::VERY_SMALL )
    {
        std::cout << "ERROR: The value of gamma = " << gamma_ << 
                     " is smaller than zero, and thus invalid.\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    if ( !declaredInitFile && toLower_(initMethod_) == InitMethodNames::FILE )
    {
        std::cout << "ERROR: initMethod is set to file, but initFile is not specified.\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

}


//-----------------------------------------------------------------------
//   readBoundaryFile_()
//-----------------------------------------------------------------------


void Config::readBoundaryFile_()
{

    std::ifstream file(boundaryFile_);

    std::cout << "Using boundary file: " << boundaryFile_ << "\n\n";

    if (file.fail())
    {
        std::cout << "---Reading boundary file failed---\n";
        std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
        std::exit(EXIT_FAILURE);
    }

    std::string word;
    std::string oldWord;
    std::string line = "start";
    std::string oldLine;

    // assume that 20 boundaries are enough at the start;
    int maxBound = 20;
    boundaries_.reinit( maxBound ); 

    int count = 0;
    while (file.peek() != EOF)
    {

        oldLine = line;

        std::getline(file,line);
        
        if ( line.empty() )
        {
            continue;
        }

        if ( ( line.back() == '{' ) )
        {
            Boundary bound;

            if ( line.length() > 1 )
            {
                int pos = line.find("{");
                bound.name = line.substr(0,pos);
            }
            else
            {
                bound.name = oldLine;
            }

            file >> word;

            while ( word.back() != '}' )
            {
                //TODO: make error that } must be on a new line


                if ( toLower_(word) == toLower_(BoundaryVarNames::KEY ) )
                {
                    file >> bound.key;
                }

                if ( toLower_(word) == toLower_(BoundaryVarNames::TYPE ) )
                {
                    file >> bound.type;
                    //TODO: check types
                }

                if ( toLower_(word) == toLower_(BoundaryVarNames::U_VELOCITY ) )
                {
                    file >> bound.u;
                }

                if ( toLower_(word) == toLower_(BoundaryVarNames::V_VELOCITY ) )
                {
                    file >> bound.v;
                }

                if ( toLower_(word) == toLower_(BoundaryVarNames::PRESSURE ) )
                {
                    file >> bound.p;
                }

                if ( toLower_(word) == toLower_(BoundaryVarNames::DENSITY ) )
                {
                    file >> bound.rho;
                }

                file >> word;
            }

            if ( count < boundaries_.rowSize() )
            {
                boundaries_[count] = bound;
            }
            else
            {
                boundaries_.resize( 2*count );
                boundaries_[count] = bound;
            }

            count++;           
        }
    }

    boundaries_.resize( count );
}


//-----------------------------------------------------------------------
//   readSelectedBCs_()
//-----------------------------------------------------------------------


void Config::readSelectedBCs_( std::string line )
{
    char        delimiter = ',';
    std::string comment   = "//";

    // first remove everything that comes after a comment if there is one
    size_t pos = 0;
    if ( (pos = line.find(comment)) != std::string::npos ) 
    {
        line.erase(pos, std::string::npos);
    }

    // Initialise vector with correct number of bc's
    int count = 1;
    for (size_t i = 0; i < line.length(); i++)
    {
        if (line[i] == delimiter)
        {
            count++; 
        } 
    }

    // std::string can not be initialised with null, so use empty string literal
    selectedBoundaries_.reinit(count, "");

    // then save the selected boundaries in array
    int idx = 0;
    while ( (pos = line.find(delimiter)) != std::string::npos )
    {
        selectedBoundaries_[idx] = line.substr(0, pos);
        idx++;

        line.erase( 0, pos + 1 ); // +1 to account for the delimiter char
    } 

    // add last entry to the vector
    selectedBoundaries_[idx] = line.substr(0, pos);

    // from the entries, remove unnecessary whitespace
    for (int i = 0; i < selectedBoundaries_.rowSize(); i++)
    {
        trim_(selectedBoundaries_[i]);
    }
}


//-----------------------------------------------------------------------
//   toLower_()
//-----------------------------------------------------------------------


std::string Config::toLower_(std::string s) 
{
    std::transform(s.begin(), s.end(), s.begin(), 
                   [](unsigned char c){ return std::tolower(c); } 
                  );
    return s;
}


//-----------------------------------------------------------------------
//   peekWord_()
//-----------------------------------------------------------------------


std::string Config::peekWord_(std::ifstream& file) 
{
  std::string    result;
  std::streampos pOrig = file.tellg();

  file >> result;
  file.seekg(pOrig);

  return result;
}


//-----------------------------------------------------------------------
//   trim_() 
//-----------------------------------------------------------------------


// trim from start (in place)
void Config::ltrim_(std::string &s) 
{
    s.erase( s.begin(), std::find_if (s.begin(), s.end(), [](unsigned char ch) 
    {
        return !std::isspace(ch);
    } ) );
}

// trim from end (in place)
void Config::rtrim_(std::string &s) 
{
    s.erase(std::find_if (s.rbegin(), s.rend(), [](unsigned char ch) 
    {
        return !std::isspace(ch);
    } ).base(), s.end());
}

// trim from both ends (in place)
void Config::trim_(std::string &s) {
    rtrim_(s);
    ltrim_(s);
}
