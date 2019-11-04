/***************************************************************************
  iPIC3D.cpp  -  Main file for 3D simulation
  -------------------
 ************************************************************************** */


// Topology
#include "processtopology/VirtualTopology3D.h"
#include "processtopology/VCtopology3D.h"
// input
#include "inputoutput/CollectiveIO.h"
#include "inputoutput/Collective.h"
// grid
#include "grids/Grid.h"
#include "grids/Grid3DCU.h"
// fields
#include "fields/Field.h"
#include "fields/EMfields3D.h"
// particles
#include "particles/Particles.h"
#include "particles/Particles3Dcomm.h"
#include "particles/Particles3D.h"

// serial ASCII output
#include "inputoutput/SerialIO.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;
using std::ofstream;


int main(int argc, char **argv) {
  Collective *col = new Collective(argc, argv); // Every proc loads the parameters of simulation from class Collective
  bool verbose = col->getVerbose();
  int restart_cycle = col->getRestartOutputCycle();
  string SaveDirName = col->getSaveDirName();
  string RestartDirName = col->getRestartDirName();
  const int restart = col->getRestart_status();
  const int ns = col->getNs();  // get the number of particle species involved in simulation
  const int first_cycle = col->getLast_cycle() + 1; // get the last cycle from the restart
  // initialize the virtual cartesian topology 
  VCtopology3D *vct = new VCtopology3D();
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology();
  // Print the initial settings to stdout and a file
  col->Print();
  col->save();
  // Create the local grid
  Grid3DCU *grid = new Grid3DCU(col, vct);  // Create the local grid
  EMfields3D *EMf = new EMfields3D(col, grid);  // Create Electromagnetic Fields Object
  //EMf->init(vct,grid);
  EMf->initGEM(vct,grid);  // init simulation for reconnection
  // Allocation of particles
  Particles3D *part = new Particles3D[ns];
  for (int i = 0; i < ns; i++)
        part[i].allocate(i, col, vct, grid);
  // MAXWELLIAN
  for (int i = 0; i < ns; i++)
      part[i].maxwellian(grid, EMf, vct); // all the species have Maxwellian distribution in the velocity
  double Eenergy = 0.0, Benergy = 0.0, TOTenergy = 0.0;
  double *Ke = new double[ns];
  string cq = SaveDirName + "/ConservedQuantities.txt";  
  
    ofstream my_file(cq.c_str());
    my_file.close();
  
 
  // *******************************************//
  // **** Start the Simulation! ***//
  // *******************************************//
  for (int cycle = first_cycle; cycle < (col->getNcycles() + first_cycle); cycle++) {
      cout << endl;
      cout << endl;
      cout << "***********************" << endl;
      cout << "   cycle = " << cycle + 1 << endl;
      cout << "***********************" << endl;
    

    // set to zero the densities - needed for interpolation
    EMf->setZeroDensities();    
    
    // PARTICLE MOVER
    for (int i = 0; i < ns; i++)  // move each species
      //part[i].mover_relativistic(grid, vct, EMf);
      part[i].mover_PC(grid, vct, EMf); // use the classical Predictor Corrector scheme
    
    // interpolation particle to grid 
    cout << "***  INTERPOLATION P->G ***" << endl;
    for (int i = 0; i < ns; i++)
      part[i].interpP2G(EMf, grid, vct);  // interpolate Particles to Grid(Nodes)
    // rest of interpolation
    EMf->sumOverSpecies(vct);   // sum all over the species
    EMf->interpDensitiesN2C(vct, grid); // calculate densities on centers from nodes
    EMf->calculateHatFunctions(grid, vct);  // calculate the hat quantities for the implicit method
    
    // MAXWELL'S SOLVER
    EMf->calculateE(grid, vct); // calculate the E field
    // calculate the B field
    EMf->calculateB(grid, vct);

    // diagnostics
    if (cycle % col->getDiagnosticsOutputCycle() == 0) {
      VTK_Write_Scalars(cycle,grid,EMf);
      VTK_Write_Vectors(cycle,grid,EMf);
      // write the conserved quantities
      Eenergy = EMf->getEenergy(); Benergy = EMf->getBenergy();
      TOTenergy = 0.0;
      for (int is = 0; is < ns; is++) {
        Ke[is] = part[is].getKe();
        TOTenergy += Ke[is]; // store kinetic energies
      }
      
      ofstream my_file(cq.c_str(), fstream::app);
      // cycle, total energy, electric field energy, magnetic field energy, total kinetic energy, electron kinetic energy, ion kinetic energy
      my_file << cycle << "\t" << (Eenergy + Benergy + TOTenergy) << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTenergy << "\t" << Ke[0] << "\t" << Ke[1] << endl;
      my_file.close();
    }
   

  }
  // deallocate
  delete[]Ke;
  return (0);

}
