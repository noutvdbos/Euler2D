#include <iostream>
#include <string>
#include <fstream>
#include <chrono>
#include <algorithm>
#include <cmath>

#include "Solver.h"    
#include "Mesh.h"
#include "Flux.h"
#include "Names.h"
#include "Thermo.h"
#include "utils/Constants.h"
#include "utils/Matrix.h"
#include "utils/Vector.h"
#include "utils/utils.h"
#include "utils/Types.h"

//=======================================================================
//   class Solver
//=======================================================================

//-----------------------------------------------------------------------
//   static data
//-----------------------------------------------------------------------


//-----------------------------------------------------------------------
//   constructor & destructor
//-----------------------------------------------------------------------

Solver::Solver( Config config)
{
    std::string meshFile = config.getMeshFile();

    dt_   = config.getTimeStep();
    time_ = config.getStartTime();

    writeInt_    = config.getWriteInterval();
    writeMethod_ = config.getWriteMethod();


    gamma_ = config.getGamma();

    mesh_ = new Mesh(meshFile);

    config_ = config;

    iter_ = 0;

    initCellStates_();

    initFaces_();

    refRho_    = config.getRefDensity();
    refP_      = config.getRefPressure();
    refVmag_   = config.getRefVelocity();
    refLength_ = config.getRefLength();


    refMach_  = refVmag_/ std::sqrt(gamma_*refP_/refRho_);
    refStagP_ = refP_* std::pow( ( 1 + (gamma_-1)/2*refMach_*refMach_ ), gamma_/(gamma_-1) );
    refDynP_  = 0.5*refRho_*refVmag_*refVmag_;

    std::cout << "reference Mach: \t" << refMach_ << "\n" <<
                 "reference Stag Pres: \t" << refStagP_ << "\n";
    
    //std::cout   << "reference dynamic pressure: " << 0.5*refRho_*refVmag_*refVmag_ << "\n"
    //            << "difference between stag ref and ref pressure: " << refStagP_ - refP_ << "\n";             


}

Solver::~Solver()
{
    delete mesh_;
}


//-----------------------------------------------------------------------
//   initCellStates_()
//-----------------------------------------------------------------------


void  Solver::initCellStates_( )
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration;  

    std::cout << "\n----------------------INITIALISATION------------------------\n";
    auto t0 = high_resolution_clock::now();
    auto t2 = t0;

    if ( config_.getInitMethod() == InitMethodNames::UNIFORM )
    {
        double u = config_.getInitUVelocity();
        double v = config_.getInitVVelocity();
        double p = config_.getInitPressure();
        double rho = config_.getInitDensity();

        double E = Thermo::getInternalEnergy(p,rho,u,v,gamma_);

        std::list<double> init = {rho, rho*u, rho*v, rho*E};

        cellStates_.reinit(mesh_->cells.rowSize(),init);

    }
    else if ( config_.getInitMethod() == InitMethodNames::FILE )
    {
        Vector<double> u;
        Vector<double> v;
        Vector<double> rho;
        Vector<double> p;
        Vector<double> e;


        std::cout << "Using initialisation file: " 
                  << config_.getInitFile() << "\n\n";

        std::cout << "---Reading initialisation file---\n";

        config_.getInitVariables( p,rho,u,v,mesh_->cells.rowSize() );
        
        t2 = high_resolution_clock::now();

        if (!(u.rowSize() == v.rowSize() && v.rowSize() == rho.rowSize() &&
              rho.rowSize() == p.rowSize() ) )
            {
                std::cout << "ERROR: input variables do not have the same length \n";
                std::cout << "lengths are: u: " << u.rowSize() << " v: " << v.rowSize()
                          << " rho: " << rho.rowSize() << " p: " << p.rowSize() << "\n";
                std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
                std::exit(EXIT_FAILURE);
            }

        Thermo::getInternalEnergy(e,p,rho,u,v,gamma_);

        // with a initialiser lists of vectors, the matrix has 4 rows, 
        // and colums as long as the rowSize of the vectors. For cellStates_,
        // this needs to be the other way around, so transpose after. 

        cellStates_ = { rho, rho*u, rho*v, rho*e};

        cellStates_ = cellStates_.transpose();

        if ( cellStates_.rowSize() != mesh_->cells.rowSize() )
        {
            std::cout << "ERROR: inital input files do not have the same "
            << " number of cells as mesh file. \n";
            std::cout << "Exiting " << ProgramNames::PROGRAM << " ... \n";
            std::exit(EXIT_FAILURE);
        } 
    }  
    
    std::cout << "Read initialisation file succesfully.\n";
    
    duration<double, std::milli> time = t2 - t0;
    std::cout << "\n-Total initialisation time: " << time.count()/1000 << "s.\n\n";
}


//-----------------------------------------------------------------------
//   initFaces_()
//-----------------------------------------------------------------------


void  Solver::initFaces_( )
{
    // Be careful of the off by one error in the faces. Cells have one more
    // column because the first column is the total column size. Faces don't
    // need this information. 
    wallCount_ = 0;

    faces_.reinit( mesh_->cells.rowSize(), mesh_->cells.colSize()-1 );

    for (int i = 0; i < mesh_->cells.rowSize(); i++)
    {
       for (int j = 1; j < mesh_->cells(i,0) +1; j++)
       {
            int cellidx1 = i;
            int cellidx2 = mesh_->elsuel(i,j);
            int bcidx    = mesh_->bcsuel(i,j);

            int node1 = mesh_->cells(i,j);
            int node2;

            if ( j < mesh_->cells(i,0) )
            {
                node2 = mesh_->cells(i, j+1);
            }
            else
            {
                node2 = mesh_->cells(i,1);
            }

            faces_[i][j-1] = Face( *mesh_,
                                 cellStates_,
                                 config_,
                                 node1,
                                 node2,
                                 cellidx1,
                                 cellidx2,
                                 bcidx,
                                 mesh_->hasNeighbour(i,j) );
            
            //store how many wall boundary faces there are, in case you want
            // to output them

            if (faces_[i][j-1].isWallBoundary() )
            {
                wallCount_ ++;
            }
       }
    }
        
}



//-----------------------------------------------------------------------
//   calcTimeStep_()
//-----------------------------------------------------------------------


void  Solver::calcTimeStep_( )
{
    if( config_.adjustTimeStep() )
    {
        dt_ *= config_.getMaxCourant()/calcMaxCourantNr();   

        if ( writeMethod_ == WriteMethodNames::TIME )    
        {
            double timeToNextWrite = ( std::floor( 
                                        ( time_ + Constants::VERY_SMALL )
                                        / writeInt_) + 1 ) * writeInt_ - time_;

            dt_ = std::max( timeToNextWrite / std::ceil(timeToNextWrite/dt_),
                            Constants::VERY_SMALL );
        }
    }
    return;
}


//-----------------------------------------------------------------------
//   updateCellStates_()
//-----------------------------------------------------------------------


void  Solver:: updateCellStates_()
{
    Matrix<double> fluxes(cellStates_.rowSize(), cellStates_.colSize() );

    for (int i = 0; i < cellStates_.rowSize(); i++)
    {
        for (int j = 1; j < mesh_->cells(i,0) +1; j++)
        {
            Vector<double> flux (cellStates_.colSize(), 0.0);

            faces_[i][j-1].getFlux( flux, false );
            fluxes[i] += flux;
        }
    }
    
    fluxes.elemwiseMultiply( -dt_/mesh_->surfaces );    
    cellStates_ += fluxes;

    return;
}


//-----------------------------------------------------------------------
//   calcMaxCourantNr()
//-----------------------------------------------------------------------


double  Solver::calcMaxCourantNr( ) const
{
    
    Vector<double> uMags(cellStates_.rowSize());

    for (int i = 0; i < uMags.rowSize(); i++)
    {
        uMags[i] = std::sqrt( 
                   ( std::pow(cellStates_[i][1],2) 
                 + std::pow(cellStates_[i][2],2) )
                 / std::pow(cellStates_[i][0],2) );
    }

    return ( dt_*(uMags/utils::sqrt(mesh_->surfaces)).max() );
}


//-----------------------------------------------------------------------
//   iterate()
//-----------------------------------------------------------------------


void  Solver::iterate( )
{
    calcTimeStep_();
    

    updateCellStates_();
    
    time_ += dt_;
    iter_ ++;

    return;
}


//-----------------------------------------------------------------------
//   getTime()
//-----------------------------------------------------------------------


double Solver::getTime() const
{
    return time_;
}


//-----------------------------------------------------------------------
//   getTimeStep()
//-----------------------------------------------------------------------


double Solver::getTimeStep() const
{
    return dt_;
}


//-----------------------------------------------------------------------
//   getIteration()
//-----------------------------------------------------------------------


int Solver::getIteration() const
{
    return iter_;
}


//-----------------------------------------------------------------------
//   hasFinished()
//-----------------------------------------------------------------------


bool Solver::hasFinished() const
{
    return (time_ >= config_.getEndTime() );
}


//-----------------------------------------------------------------------
//   getResults()
//-----------------------------------------------------------------------


void Solver::getResults( Matrix<double>& results ) const
{
    Matrix<double> temp (cellStates_.colSize(), cellStates_.rowSize());

    utils::transpose(cellStates_,temp);

    Vector<double> uMags(cellStates_.rowSize());    
    Vector<double> machs(cellStates_.rowSize());

    for (int i = 0; i < cellStates_.rowSize(); i++)
    {
        uMags[i] = std::sqrt( 
                   ( std::pow(cellStates_[i][1],2) 
                 + std::pow(cellStates_[i][2],2) )
                 / std::pow(cellStates_[i][0],2) );

        machs[i] = uMags[i]/std::sqrt(gamma_* Thermo::getPressure(cellStates_[i],gamma_)/temp[0][i]);
    }


    results =  
    { 
        temp[0],                        	                                       // density
        temp[1]/temp[0],                                                           // u velocity
        temp[2]/temp[0],                                                           // v velocity
        Thermo::getPressure(cellStates_,gamma_),                                   // pressure
        (Thermo::getPressure(cellStates_,gamma_)- refP_) / refDynP_,               // cp
        machs                                                                      // Mach
    };
    
}


//-----------------------------------------------------------------------
//   getWallCoeffs()
//-----------------------------------------------------------------------


void Solver::getWallCoeffs
( 
    Matrix<double>& results,
    double& cfx,
    double& cfy,
    std::string boundary ) const
{
    //force coefficients are calculated via Riemann sum
    cfx = 0.0; // force coefficient in x direction
    cfy = 0.0; // force coefficient in y direction
    
    //write out the pressure plus coordinates, i.e. 4 variables
    //initialise with total wall count, then later resize to only current boundary
    results.reinit(wallCount_, 4 );

    int wallidx = 0;
    for (int i = 0; i < cellStates_.rowSize(); i++)
    {
        for (int j = 1; j < mesh_->cells(i,0) +1; j++)
        {
            Vector<double> result (results.colSize(), 0.0);
            Face face = faces_[i][j-1];

            if (face.isWallBoundary() && face.getBoundary().name == boundary )
            {
                face.getWallPressure(result);
                results[wallidx] = result; 

                //convert pressure to cp
                results[wallidx][3] = (results[wallidx][3]- refP_) / refDynP_;

                //add contribution to cfx and cfy, cfx is multiplied by x normal
                // and cfy is multiplied by y normal
                cfx += results[wallidx][3] * face.getHeight()*face.getNormal()[0]; 
                cfy += results[wallidx][3] * face.getHeight()*face.getNormal()[1]; 

                wallidx++;
            }
        }
    }

    //make sure force coefficients are scaled by reference length
    cfx /= refLength_;
    cfy /= refLength_;

    if ( wallidx !=wallCount_ )
    {
        results.resize(wallidx);
    }
}
