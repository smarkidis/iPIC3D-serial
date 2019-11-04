#ifndef Collective_H
#define Collective_H



#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "../ConfigFile/src/ConfigFile.h"
#include "../ConfigFile/src/input_array.h"


#include "CollectiveIO.h"

using std::cout;
using std::endl;
using std::ofstream;
using namespace std;

class Collective:public CollectiveIO {
public:
    /** constructor: initialize physical parameters with values */
    Collective(int argc, char **argv);
    /** destructor */
    ~Collective();
    /** read input file */
    void ReadInput(string inputfile);
    /** read the restart input file from HDF5 */
    int ReadRestart(string inputfile);
    /** Print physical parameters */
    void Print();
    /** save setting in a file */
    void save();
    /** get the physical space dimensions            */
    int getDim();
    /** Get length of the system - direction X */
    double getLx();
    /** Get length of the system - direction Y */
    double getLy();
    /** Get length of the system - direction Z */
    double getLz();
    /** Get object center - direction X */
    double getx_center();
    /** Get object center - direction Y */
    double gety_center();
    /** Get object center - direction Z */
    double getz_center();
    /** Get object size - cubic box */
    double getL_square();
    /** Get the number of cells - direction X*/
    int getNxc();
    /** Get the number of cells - direction Y*/
    int getNyc();
    /** Get the number of cells - direction Z*/
    int getNzc();
    /** Get the grid spacing - direction X*/
    double getDx();
    /** Get the grid spacing - direction Y*/
    double getDy();
    /** Get the grid spacing - direction Z*/
    double getDz();
    /** get the light speed */
    double getC();
    /** get the time step */
    double getDt();
    /** get the decentering parameter */
    double getTh();
    /** get the Smoothing value*/
    double getSmooth();
    /** get the number of time cycles */
    int getNcycles();
    /** get the number of species */
    int getNs();
    /** get the number of particles for different species */
    long getNp(int nspecies);
    /** get the number of particles per cell  */
    int getNpcel(int nspecies);
    /** get the number of particles per cell - direction X  */
    int getNpcelx(int nspecies);
    /** get the number of particles per cell - direction Y  */
    int getNpcely(int nspecies);
    /** get the number of particles per cell - direction Z  */
    int getNpcelz(int nspecies);
    /** get maximum number of particles for different species */
    long getNpMax(int nspecies);
    /** NpMax/Np is the ratio between the maximum number of particles allowed on a processor and the number of particles*/
    double getNpMaxNpRatio();
    /** get charge to mass ratio for different species */
    double getQOM(int nspecies);
    /** get background charge for GEM challenge */
    double getRHOinit(int nspecies);
    /** get thermal velocity  - X direction    */
    double getUth(int nspecies);
    /** get thermal velocity  - Y direction    */
    double getVth(int nspecies);
    /** get thermal velocity  - Z direction    */
    double getWth(int nspecies);
    /** get Drift velocity - Direction X         */
    double getU0(int nspecies);
    /** get Drift velocity - Direction Y         */
    double getV0(int nspecies);
    /** get Drift velocity - Direction Z         */
    double getW0(int nspecies);
    /** get the boolean value for TrackParticleID */
    bool getTrackParticleID(int nspecies);
    /** get SaveDirName  */
    string getSaveDirName();
    /** get last_cycle  */
    int getLast_cycle();
    /** get RestartDirName  */
    string getRestartDirName();
    
    /** get Boundary Condition Particles: FaceXright */
    int getBcPfaceXright();
    /** get Boundary Condition Particles: FaceXleft */
    int getBcPfaceXleft();
    /** get Boundary Condition Particles: FaceYright */
    int getBcPfaceYright();
    /** get Boundary Condition Particles: FaceYleft */
    int getBcPfaceYleft();
    /** get Boundary Condition Particles: FaceYright */
    int getBcPfaceZright();
    /** get Boundary Condition Particles: FaceYleft */
    int getBcPfaceZleft();
    
    /** get Boundary Condition Electrostatic Potential: FaceXright */
    int getBcPHIfaceXright();
    /** get Boundary Condition Electrostatic Potential:FaceXleft */
    int getBcPHIfaceXleft();
    /** get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceYright();
    /** get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceYleft();
    /** get Boundary Condition Electrostatic Potential:FaceYright */
    int getBcPHIfaceZright();
    /** get Boundary Condition Electrostatic Potential:FaceYleft */
    int getBcPHIfaceZleft();
    
    /** get Boundary ConditionElectric Field: FaceXright */
    int getBcEMfaceXright();
    /** get Boundary Condition Electric Field: FaceXleft */
    int getBcEMfaceXleft();
    /** get Boundary Condition Electric Field: FaceYright */
    int getBcEMfaceYright();
    /** get Boundary Condition Electric Field: FaceYleft */
    int getBcEMfaceYleft();
    /** get Boundary Condition Electric Field: FaceZright */
    int getBcEMfaceZright();
    /** get Boundary Condition Electric Field: FaceZleft */
    int getBcEMfaceZleft();
    
    /** get RESTART */
    int getRestart_status();
    
    /** get the sheet thickness */
    double getDelta();
    /** get the amplitude of the magnetic field along x */
    double getB0x();
    /** get the amplitude of the magnetic field along y */
    double getB0y();
    /** get the amplitude of the magnetic field along z */
    double getB0z();
    
    /** get the boolean value for verbose results */
    bool getVerbose();
    /** get the PoissonCorrection */
    bool getPoissonCorrection();
    
    /** get the velocity of injection of the plasma from the wall */
    double getVinj();
    
    /** get the converging tolerance for CG solver */
    double getCGtol();
    /** get the converging tolerance for GMRES solver */
    double getGMREStol();
    /** get the numbers of iteration for the PC mover */
    int getNiterMover();
    
    /** output of fields */
    int getFieldOutputCycle();
    /** output of particles */
    int getParticlesOutputCycle();
    /** output of restart */
    int getRestartOutputCycle();
    /** output of diagnostics */
    int getDiagnosticsOutputCycle();
    /** get four pi value */
    double getFourPI();
    /** get Nwaves */
    int getNwaves();
    /** get the Number of subcycles */
    int getN_sub_cycles();
    /** get dBoB0  */
    double getdBoB0();
    /** get Wave file */
    string getWaveFile();
    
    /** Test particles */
    double getEnergy();
    double getPitchAngle();
    
    
private:
    /** inputfile */
    string inputfile;
    /** light speed */
    double c;
    /** 4  pi */
    double fourpi;
    /** time step */
    double dt;
    /** decentering parameter */
    double th;
    /** Smoothing value*/
    double Smooth;
    /** number of time cycles */
    int ncycles;
    /** physical space dimensions            */
    int dim;
    /** simulation box length - X direction   */
    double Lx;
    /** simulation box length - Y direction   */
    double Ly;
    /** simulation box length - Z direction   */
    double Lz;
    /** object center - X direction   */
    double x_center;
    /** object center - Y direction   */
    double y_center;
    /** object center - Z direction   */
    double z_center;
    /** object size - assuming a cubic box  */
    double L_square;
    /** number of cells - X direction        */
    int nxc;
    /** number of cells - Y direction        */
    int nyc;
    /** number of cells - Z direction        */
    int nzc;
    /** grid spacing - X direction */
    double dx;
    /** grid spacing - Y direction */
    double dy;
    /** grid spacing - Z direction */
    double dz;
    /** number of species */
    int ns;
    /** number of particles per cell */
    int *npcel;
    /** number of particles per cell - X direction */
    int *npcelx;
    /** number of particles per cell - Y direction */
    int *npcely;
    /** number of particles per cell - Z direction */
    int *npcelz;
    /** number of particles array for different species */
    long *np;
    /** maximum number of particles array for different species */
    long *npMax;
    /** max number of particles */
    double NpMaxNpRatio;
    /** charge to mass ratio array for different species */
    double *qom;
    /** charge to mass ratio array for different species */
    double *rhoINIT;
    /** thermal velocity  - Direction X  */
    double *uth;
    /** thermal velocity  - Direction Y  */
    double *vth;
    /** thermal velocity  - Direction Z  */
    double *wth;
    /** Drift velocity - Direction X     */
    double *u0;
    /** Drift velocity - Direction Y    */
    double *v0;
    /** Drift velocity - Direction Z     */
    double *w0;
    /** number of particle of subcycles */
    int n_sub_cycles;
    
    /** TrackParticleID     */
    bool *TrackParticleID;
    /** SaveDirName     */
    string SaveDirName;
    /** RestartDirName     */
    string RestartDirName;
    /** get inputfile  */
    string getinputfile();
    /** restart_status  0 --> no restart; 1--> restart, create new; 2--> restart, append; */
    int restart_status;
    /** last cycle */
    int last_cycle;
    
    
    /** pitch angle */
    double pitch_angle;
    /** energy of the particle */
    double energy;
    
    
    /** Boundary condition on particles
     0 = exit
     1 = perfect mirror
     2 = riemission
     */
    /** Boundary Condition Particles: FaceXright */
    int bcPfaceXright;
    /** Boundary Condition Particles: FaceXleft */
    int bcPfaceXleft;
    /** Boundary Condition Particles: FaceYright */
    int bcPfaceYright;
    /** Boundary Condition Particles: FaceYleft */
    int bcPfaceYleft;
    /** Boundary Condition Particles: FaceYright */
    int bcPfaceZright;
    /** Boundary Condition Particles: FaceYleft */
    int bcPfaceZleft;
    
    
    /** Field Boundary Condition
     0 =  Dirichlet Boundary Condition: specifies the valueto take pn the boundary of the domain
     1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
     2 =  Periodic Condition
     */
    /** Boundary Condition Electrostatic Potential: FaceXright */
    int bcPHIfaceXright;
    /** Boundary Condition Electrostatic Potential:FaceXleft */
    int bcPHIfaceXleft;
    /** Boundary Condition Electrostatic Potential:FaceYright */
    int bcPHIfaceYright;
    /** Boundary Condition Electrostatic Potential:FaceYleft */
    int bcPHIfaceYleft;
    /** Boundary Condition Electrostatic Potential:FaceZright */
    int bcPHIfaceZright;
    /** Boundary Condition Electrostatic Potential:FaceZleft */
    int bcPHIfaceZleft;
    
    /** Boundary Condition EM Field: FaceXright */
    int bcEMfaceXright;
    /** Boundary Condition EM Field: FaceXleft */
    int bcEMfaceXleft;
    /** Boundary Condition EM Field: FaceYright */
    int bcEMfaceYright;
    /** Boundary Condition EM Field: FaceYleft */
    int bcEMfaceYleft;
    /** Boundary Condition EM Field: FaceZright */
    int bcEMfaceZright;
    /** Boundary Condition EM Field: FaceZleft */
    int bcEMfaceZleft;
    
    
    /** GEM Challenge parameters */
    /** current sheet thickness */
    double delta;
    /* Amplitude of the field */
    double B0x;
    double B0y;
    double B0z;
    
    
    /** boolean value for verbose results */
    bool verbose;
    /** RESTART */
    bool RESTART1;
    /** Poisson Correction */
    bool PoissonCorrection;
    
    /** velocity of the injection from the wall */
    double Vinj;
    
    /** CG solver stopping criterium tolerance */
    double CGtol;
    /** GMRES solver stopping criterium tolerance */
    double GMREStol;
    /** mover predictor correcto iteration */
    int NiterMover;
    
    /** Output for field */
    int FieldOutputCycle;
    /** Output for particles */
    int ParticlesOutputCycle;
    /** restart cycle */
    int RestartOutputCycle;
    /** Output for diagnostics */
    int DiagnosticsOutputCycle;
    /** Number of waves present in the system */
    int Nwaves;
    /** Perturbation amplitude */
    double dBoB0;
    /** name of the file with wave amplitude and phases */
    string WaveFile;
    
};
/** Read the input file from text file and put the data in a collective wrapper:
 if it's a restart read from input file basic sim data and load particles and EM field
 from restart file */
inline void Collective::ReadInput(string inputfile) {
    using namespace std;
    int test_verbose;
    // Loading the input file
    ConfigFile config(inputfile);
    // the following variables are ALWAYS taken from inputfile, even if restarting
    {
        dt = config.read < double >("dt");
        ncycles = config.read < int >("ncycles");
        th = config.read < double >("th",1);
        Smooth = config.read < double >("Smooth",0.2); // Default 0.2
        SaveDirName = config.read < string > ("SaveDirName");
        RestartDirName = config.read < string > ("RestartDirName");
        ns = config.read < int >("ns");
        NpMaxNpRatio = config.read < double >("NpMaxNpRatio",1);
        // GEM Challenge
        B0x = config.read < double >("B0x",0.0);
        B0y = config.read < double >("B0y",0.0);
        B0z = config.read < double >("B0z",0.0);
        delta = config.read < double >("delta",0.5);
        
        rhoINIT = new double[ns];
        array_double rhoINIT0 = config.read < array_double > ("rhoINIT");
        rhoINIT[0] = rhoINIT0.a;
        if (ns > 1)
            rhoINIT[1] = rhoINIT0.b;
        if (ns > 2)
            rhoINIT[2] = rhoINIT0.c;
        if (ns > 3)
            rhoINIT[3] = rhoINIT0.d;
        if (ns > 4)
            rhoINIT[4] = rhoINIT0.e;
        if (ns > 5)
            rhoINIT[5] = rhoINIT0.f;
        // take the tolerance of the solvers
        CGtol = config.read < double >("CGtol",1E-3);
        GMREStol = config.read < double >("GMREStol",1E-3);
        NiterMover = config.read < int >("NiterMover",3);
        // take the injection of the particless
        Vinj = config.read < double >("Vinj",0.0);
        
        // take the output cycles
        FieldOutputCycle = config.read < int >("FieldOutputCycle",10);
        ParticlesOutputCycle = config.read < int >("ParticlesOutputCycle",100);
        RestartOutputCycle = config.read < int >("RestartOutputCycle",100000);
        DiagnosticsOutputCycle = config.read < int >("DiagnosticsOutputCycle", FieldOutputCycle);
        fourpi = config.read <double>("fourpi", 12.5663706);
        // Waves info
        Nwaves = config.read <int>("Nwaves", 1);
        dBoB0 = config.read <double>("dBoB0", 0.0);
        WaveFile = config.read < string > ("WaveFile","WaveFile.txt");
        // Poisson Correction
        PoissonCorrection = config.read <bool>("PoissonCorrection",true);
        //
        c = config.read < double >("c",1);
        x_center = config.read < double >("x_center",2.5);
        y_center = config.read < double >("y_center",2.5);
        z_center = config.read < double >("z_center",2.5);
        L_square = config.read < double >("L_square",1.25);
        
        energy = config.read < double >("energy",0.018199864696222);
        pitch_angle = config.read < double >("pitch_angle",0.698131700797732); // 40 degree default
        
        n_sub_cycles = config.read <int>("n_sub_cycles", 1);
        
        
        verbose = config.read < bool > ("verbose",true);
        
        // PHI Electrostatic Potential
        bcPHIfaceXright = config.read < int >("bcPHIfaceXright",1);
        bcPHIfaceXleft = config.read < int >("bcPHIfaceXleft",1);
        bcPHIfaceYright = config.read < int >("bcPHIfaceYright",1);
        bcPHIfaceYleft = config.read < int >("bcPHIfaceYleft",1);
        bcPHIfaceZright = config.read < int >("bcPHIfaceZright",1);
        bcPHIfaceZleft = config.read < int >("bcPHIfaceZleft",1);
        
        // EM field boundary condition
        bcEMfaceXright = config.read < int >("bcEMfaceXright",1);
        bcEMfaceXleft = config.read < int >("bcEMfaceXleft",1);
        bcEMfaceYright = config.read < int >("bcEMfaceYright",1);
        bcEMfaceYleft = config.read < int >("bcEMfaceYleft",1);
        bcEMfaceZright = config.read < int >("bcEMfaceZright",1);
        bcEMfaceZleft = config.read < int >("bcEMfaceZleft",1);
        
        // Particles Boundary condition
        bcPfaceXright = config.read < int >("bcPfaceXright",1);
        bcPfaceXleft = config.read < int >("bcPfaceXleft",1);
        bcPfaceYright = config.read < int >("bcPfaceYright",1);
        bcPfaceYleft = config.read < int >("bcPfaceYleft",1);
        bcPfaceZright = config.read < int >("bcPfaceZright",1);
        bcPfaceZleft = config.read < int >("bcPfaceZleft",1);
        
    }
    
    if (RESTART1) {               // you are restarting
        RestartDirName = config.read < string > ("RestartDirName");
        ReadRestart(RestartDirName);
    }
    else {
        restart_status = 0;
        last_cycle = -1;
        
        Lx = config.read < double >("Lx",1);
        Ly = config.read < double >("Ly",1);
        Lz = config.read < double >("Lz",1);
        nxc = config.read < int >("nxc",1);
        nyc = config.read < int >("nyc",1);
        nzc = config.read < int >("nzc",1);
        
        n_sub_cycles = config.read <int>("n_sub_cycles", 1);
        
        npcelx = new int[ns];
        npcely = new int[ns];
        npcelz = new int[ns];
        qom = new double[ns];
        uth = new double[ns];
        vth = new double[ns];
        wth = new double[ns];
        u0 = new double[ns];
        v0 = new double[ns];
        w0 = new double[ns];
        
        array_int npcelx0 = config.read < array_int > ("npcelx");
        array_int npcely0 = config.read < array_int > ("npcely");
        array_int npcelz0 = config.read < array_int > ("npcelz");
        array_double qom0 = config.read < array_double > ("qom");
        array_double uth0 = config.read < array_double > ("uth");
        array_double vth0 = config.read < array_double > ("vth");
        array_double wth0 = config.read < array_double > ("wth");
        array_double u00 = config.read < array_double > ("u0");
        array_double v00 = config.read < array_double > ("v0");
        array_double w00 = config.read < array_double > ("w0");
        
        npcelx[0] = npcelx0.a;
        npcely[0] = npcely0.a;
        npcelz[0] = npcelz0.a;
        qom[0] = qom0.a;
        uth[0] = uth0.a;
        vth[0] = vth0.a;
        wth[0] = wth0.a;
        u0[0] = u00.a;
        v0[0] = v00.a;
        w0[0] = w00.a;
        
        if (ns > 1) {
            npcelx[1] = npcelx0.b;
            npcely[1] = npcely0.b;
            npcelz[1] = npcelz0.b;
            qom[1] = qom0.b;
            uth[1] = uth0.b;
            vth[1] = vth0.b;
            wth[1] = wth0.b;
            u0[1] = u00.b;
            v0[1] = v00.b;
            w0[1] = w00.b;
        }
        if (ns > 2) {
            npcelx[2] = npcelx0.c;
            npcely[2] = npcely0.c;
            npcelz[2] = npcelz0.c;
            qom[2] = qom0.c;
            uth[2] = uth0.c;
            vth[2] = vth0.c;
            wth[2] = wth0.c;
            u0[2] = u00.c;
            v0[2] = v00.c;
            w0[2] = w00.c;
        }
        if (ns > 3) {
            npcelx[3] = npcelx0.d;
            npcely[3] = npcely0.d;
            npcelz[3] = npcelz0.d;
            qom[3] = qom0.d;
            uth[3] = uth0.d;
            vth[3] = vth0.d;
            wth[3] = wth0.d;
            u0[3] = u00.d;
            v0[3] = v00.d;
            w0[3] = w00.d;
        }
        if (ns > 4) {
            npcelx[4] = npcelx0.e;
            npcely[4] = npcely0.e;
            npcelz[4] = npcelz0.e;
            qom[4] = qom0.e;
            uth[4] = uth0.e;
            vth[4] = vth0.e;
            wth[4] = wth0.e;
            u0[4] = u00.e;
            v0[4] = v00.e;
            w0[4] = w00.e;
        }
        if (ns > 5) {
            npcelx[5] = npcelx0.f;
            npcely[5] = npcely0.f;
            npcelz[5] = npcelz0.f;
            qom[5] = qom0.f;
            uth[5] = uth0.f;
            vth[5] = vth0.f;
            wth[5] = wth0.f;
            u0[5] = u00.f;
            v0[5] = v00.f;
            w0[1] = w00.f;
        }
        
        
        
        
        
        
        
    }
    TrackParticleID = new bool[ns];
    array_bool TrackParticleID0 = config.read < array_bool > ("TrackParticleID");
    TrackParticleID[0] = TrackParticleID0.a;
    if (ns > 1)
        TrackParticleID[1] = TrackParticleID0.b;
    if (ns > 2)
        TrackParticleID[2] = TrackParticleID0.c;
    if (ns > 3)
        TrackParticleID[3] = TrackParticleID0.d;
    if (ns > 4)
        TrackParticleID[4] = TrackParticleID0.e;
    if (ns > 5)
        TrackParticleID[5] = TrackParticleID0.f;
    
    
    
    
}
/**
 *  Read the collective information from the RESTART file in HDF5 format
 *
 *  There are three restart status:
 *  restart_status = 0 --->  new inputfile
 *  restart_status = 1 --->  RESTART and restart and result directories does not coincide
 *  restart_status = 2 --->  RESTART and restart and result directories coincide
 *
 */
inline int Collective::ReadRestart(string inputfile) {
    restart_status = 1;
    return (0);
}
/** constructor */
inline Collective::Collective(int argc, char **argv) {
    if (argc < 2) {
        inputfile = "inputfile";
        RESTART1 = false;
    }
    else if (argc < 3) {
        inputfile = argv[1];
        RESTART1 = false;
    }
    else {
        if (strcmp(argv[1], "restart") == 0) {
            inputfile = argv[2];
            RESTART1 = true;
        }
        else if (strcmp(argv[2], "restart") == 0) {
            inputfile = argv[1];
            RESTART1 = true;
        }
        else {
            cout << "Error: syntax error in mpirun arguments. Did you mean to 'restart' ?" << endl;
            return;
        }
    }
    
    ReadInput(inputfile);
    
    /** dx = space step - X direction */
    dx = Lx / (double) nxc;
    /** dy = space step - Y direction */
    dy = Ly / (double) nyc;
    /** dz = space step - Z direction */
    dz = Lz / (double) nzc;
    /** npcel = number of particles per cell */
    npcel = new int[ns];
    /** np = number of particles of different species */
    np = new long[ns];
    /** npMax = maximum number of particles of different species */
    npMax = new long[ns];
    
    for (int i = 0; i < ns; i++) {
        npcel[i] = npcelx[i] * npcely[i] * npcelz[i];
        np[i] = npcel[i] * nxc * nyc * nzc;
        npMax[i] = (long) (NpMaxNpRatio * np[i]);
    }
    
    
    
}

/** destructor */
inline Collective::~Collective() {
    delete[]np;
    delete[]npcel;
    delete[]npcelx;
    delete[]npcely;
    delete[]npcelz;
    delete[]npMax;
    delete[]qom;
    
    delete[]uth;
    delete[]vth;
    delete[]wth;
    
    delete[]u0;
    delete[]v0;
    delete[]w0;
    
    delete[]TrackParticleID;
    
    delete[]rhoINIT;
    
}
/** Print Simulation Parameters */
inline void Collective::Print() {
    cout << endl;
    cout << "Simulation Parameters" << endl;
    cout << "---------------------" << endl;
    cout << "Number of species    = " << ns << endl;
    for (int i = 0; i < ns; i++)
        cout << "Number of particles of species " << i << " = " << np[i] << "\t (MAX = " << npMax[i] << ")" << "  QOM = " << qom[i] << endl;
    cout << "x-Length                 = " << Lx << endl;
    cout << "y-Length                 = " << Ly << endl;
    cout << "z-Length                 = " << Lz << endl;
    cout << "Number of cells (x)      = " << nxc << endl;
    cout << "Number of cells (y)      = " << nyc << endl;
    cout << "Number of cells (z)      = " << nzc << endl;
    cout << "Time step                = " << dt << endl;
    cout << "Number of cycles         = " << ncycles << endl;
    cout << "Results saved in: " << SaveDirName << endl;
    
    
}
/** Print Simulation Parameters */
inline void Collective::save() {
    string temp;
    temp = SaveDirName + "/SimulationData.txt";
    ofstream my_file(temp.c_str());
    my_file << "---------------------------" << endl;
    my_file << "-  Simulation Parameters  -" << endl;
    my_file << "---------------------------" << endl;
    
    my_file << "Number of species    = " << ns << endl;
    for (int i = 0; i < ns; i++)
        my_file << "Number of particles of species " << i << " = " << np[i] << "\t (MAX = " << npMax[i] << ")" << "  QOM = " << qom[i] << endl;
    my_file << "---------------------------" << endl;
    my_file << "x-Length                 = " << Lx << endl;
    my_file << "y-Length                 = " << Ly << endl;
    my_file << "z-Length                 = " << Lz << endl;
    my_file << "Number of cells (x)      = " << nxc << endl;
    my_file << "Number of cells (y)      = " << nyc << endl;
    my_file << "Number of cells (z)      = " << nzc << endl;
    my_file << "---------------------------" << endl;
    my_file << "Time step                = " << dt << endl;
    my_file << "Number of cycles         = " << ncycles << endl;
    my_file << "---------------------------" << endl;
    for (int is = 0; is < ns; is++)
        my_file << "rho init species " << is << " = " << rhoINIT[is] << endl;
    my_file << "current sheet thickness  = " << delta << endl;
    my_file << "B0x                      = " << B0x << endl;
    my_file << "BOy                      = " << B0y << endl;
    my_file << "B0z                      = " << B0z << endl;
    my_file << "---------------------------" << endl;
    my_file << "Smooth                   = " << Smooth << endl;
    my_file << "GMRES error tolerance    = " << GMREStol << endl;
    my_file << "CG error tolerance       = " << CGtol << endl;
    my_file << "Mover error tolerance    = " << NiterMover << endl;
    my_file << "---------------------------" << endl;
    my_file << "Results saved in: " << SaveDirName << endl;
    my_file << "Restart saved in: " << RestartDirName << endl;
    my_file << "---------------------" << endl;
    my_file.close();
    
}

/** get the physical space dimensions            */
inline int Collective::getDim() {
    return (dim);
}
/** get Lx */
inline double Collective::getLx() {
    return (Lx);
}
/** get Ly */
inline double Collective::getLy() {
    return (Ly);
}
/** get Lz */
inline double Collective::getLz() {
    return (Lz);
}
/** get x_center */
inline double Collective::getx_center() {
    return (x_center);
}
/** get y_center */
inline double Collective::gety_center() {
    return (y_center);
}
/** get z_center */
inline double Collective::getz_center() {
    return (z_center);
}
/** get L_square */
inline double Collective::getL_square() {
    return (L_square);
}
/** get nxc */
inline int Collective::getNxc() {
    return (nxc);
}
/** get nyx */
inline int Collective::getNyc() {
    return (nyc);
}
/** get nzc */
inline int Collective::getNzc() {
    return (nzc);
}
/** get dx */
inline double Collective::getDx() {
    return (dx);
}
/** get dy */
inline double Collective::getDy() {
    return (dy);
}
/** get dz */
inline double Collective::getDz() {
    return (dz);
}
/** get the light speed */
inline double Collective::getC() {
    return (c);
}
/** get the time step */
inline double Collective::getDt() {
    return (dt);
}
/** get the decentering parameter */
inline double Collective::getTh() {
    return (th);
}
/** get the smooth parameter */
inline double Collective::getSmooth() {
    return (Smooth);
}
/** get energy **/
inline double Collective::getEnergy() {
    return (energy);
}
/** get pitch angle **/
inline double Collective::getPitchAngle() {
    return (pitch_angle);
}

/** get the number of mover subcycling **/
inline int Collective::getN_sub_cycles(){
    return (n_sub_cycles);
}

/** get the number of time cycles */
inline int Collective::getNcycles() {
    return (ncycles);
}
/** get the number of species */
inline int Collective::getNs() {
    return (ns);
}
/** get the number of particles per cell for species nspecies */
inline int Collective::getNpcel(int nspecies) {
    return (npcel[nspecies]);
}
/** get the number of particles per cell for species nspecies - direction X */
inline int Collective::getNpcelx(int nspecies) {
    return (npcelx[nspecies]);
}
/** get the number of particles per cell for species nspecies - direction Y */
inline int Collective::getNpcely(int nspecies) {
    return (npcely[nspecies]);
}
/** get the number of particles per cell for species nspecies - direction Z */
inline int Collective::getNpcelz(int nspecies) {
    return (npcelz[nspecies]);
}
/** get the number of particles for different species */
inline long Collective::getNp(int nspecies) {
    return (np[nspecies]);
}
/** get maximum number of particles for different species */
inline long Collective::getNpMax(int nspecies) {
    return (npMax[nspecies]);
}
inline double Collective::getNpMaxNpRatio() {
    return (NpMaxNpRatio);
}
/** get charge to mass ratio for different species */
inline double Collective::getQOM(int nspecies) {
    return (qom[nspecies]);
}
/** get the background density for GEM challenge */
inline double Collective::getRHOinit(int nspecies) {
    return (rhoINIT[nspecies]);
}
/** get thermal velocity  - Direction X     */
inline double Collective::getUth(int nspecies) {
    return (uth[nspecies]);
}
/** get thermal velocity  - Direction Y     */
inline double Collective::getVth(int nspecies) {
    return (vth[nspecies]);
}
/** get thermal velocity  - Direction Z     */
inline double Collective::getWth(int nspecies) {
    return (wth[nspecies]);
}
/** get beam velocity - Direction X        */
inline double Collective::getU0(int nspecies) {
    return (u0[nspecies]);
}
/** get beam velocity - Direction Y        */
inline double Collective::getV0(int nspecies) {
    return (v0[nspecies]);
}
/** get beam velocity - Direction Z        */
inline double Collective::getW0(int nspecies) {
    return (w0[nspecies]);
}
/** get Boundary Condition Particles: FaceXright */
inline int Collective::getBcPfaceXright() {
    return (bcPfaceXright);
}
/** get Boundary Condition Particles: FaceXleft */
inline int Collective::getBcPfaceXleft() {
    return (bcPfaceXleft);
}
/** get Boundary Condition Particles: FaceYright */
inline int Collective::getBcPfaceYright() {
    return (bcPfaceYright);
}
/** get Boundary Condition Particles: FaceYleft */
inline int Collective::getBcPfaceYleft() {
    return (bcPfaceYleft);
}
/** get Boundary Condition Particles: FaceZright */
inline int Collective::getBcPfaceZright() {
    return (bcPfaceZright);
}
/** get Boundary Condition Particles: FaceZleft */
inline int Collective::getBcPfaceZleft() {
    return (bcPfaceZleft);
}
/** get Boundary Condition Electrostatic Potential: FaceXright */
inline int Collective::getBcPHIfaceXright() {
    return (bcPHIfaceXright);
}
/** get Boundary Condition Electrostatic Potential:FaceXleft */
inline int Collective::getBcPHIfaceXleft() {
    return (bcPHIfaceXleft);
}
/** get Boundary Condition Electrostatic Potential:FaceYright */
inline int Collective::getBcPHIfaceYright() {
    return (bcPHIfaceYright);
}
/** get Boundary Condition Electrostatic Potential:FaceYleft */
inline int Collective::getBcPHIfaceYleft() {
    return (bcPHIfaceYleft);
}
/** get Boundary Condition Electrostatic Potential:FaceZright */
inline int Collective::getBcPHIfaceZright() {
    return (bcPHIfaceZright);
}
/** get Boundary Condition Electrostatic Potential:FaceZleft */
inline int Collective::getBcPHIfaceZleft() {
    return (bcPHIfaceZleft);
}
/** get Boundary Condition EM Field: FaceXright */
inline int Collective::getBcEMfaceXright() {
    return (bcEMfaceXright);
}
/** get Boundary Condition EM Field: FaceXleft */
inline int Collective::getBcEMfaceXleft() {
    return (bcEMfaceXleft);
}
/** get Boundary Condition EM Field: FaceYright */
inline int Collective::getBcEMfaceYright() {
    return (bcEMfaceYright);
}
/** get Boundary Condition EM Field: FaceYleft */
inline int Collective::getBcEMfaceYleft() {
    return (bcEMfaceYleft);
}
/** get Boundary Condition EM Field: FaceZright */
inline int Collective::getBcEMfaceZright() {
    return (bcEMfaceZright);
}
/** get Boundary Condition EM Field: FaceZleft */
inline int Collective::getBcEMfaceZleft() {
    return (bcEMfaceZleft);
}
/** Get GEM Challenge parameters */
inline double Collective::getDelta() {
    return (delta);
}
inline double Collective::getB0x() {
    return (B0x);
}
inline double Collective::getB0y() {
    return (B0y);
}
inline double Collective::getB0z() {
    return (B0z);
}
/** get the boolean value for verbose results */
inline bool Collective::getVerbose() {
    return (verbose);
}
/** get the boolean value for TrackParticleID */
inline bool Collective::getTrackParticleID(int nspecies) {
    return (TrackParticleID[nspecies]);
}
inline int Collective::getRestart_status() {
    return (restart_status);
}
/** get SaveDirName  */
inline string Collective::getSaveDirName() {
    return (SaveDirName);
}
/** get RestartDirName  */
inline string Collective::getRestartDirName() {
    return (RestartDirName);
}
/** get inputfile  */
inline string Collective::getinputfile() {
    return (inputfile);
}
/** get last_cycle  */
inline int Collective::getLast_cycle() {
    return (last_cycle);
}
/** get the velocity of injection of the plasma from the wall */
inline double Collective::getVinj() {
    return (Vinj);
}
/** get the converging tolerance for CG solver */
inline double Collective::getCGtol() {
    return (CGtol);
}
/** get the converging tolerance for GMRES solver */
inline double Collective::getGMREStol() {
    return (GMREStol);
}
/** get the numbers of iteration for the PC mover */
inline int Collective::getNiterMover() {
    return (NiterMover);
}
/** output of fields */
inline int Collective::getFieldOutputCycle() {
    return (FieldOutputCycle);
}
/** output of particles */
inline int Collective::getParticlesOutputCycle() {
    return (ParticlesOutputCycle);
}
/** restart cycle */
inline int Collective::getRestartOutputCycle() {
    return (RestartOutputCycle);
}
/** output of fields */
inline int Collective::getDiagnosticsOutputCycle() {
    return (DiagnosticsOutputCycle);
}
/** get the four pi value */
inline double Collective::getFourPI() {
    return (fourpi);
}
/** get the number of waves */
inline int Collective::getNwaves() {
    return (Nwaves);
}
/** get the number of amplitude of perturbation*/
inline double Collective::getdBoB0() {
    return (dBoB0);
}
/** get the name of the Wave files*/
inline string Collective::getWaveFile() {
    return (WaveFile);
}
/** get the name of the Wave files*/
inline bool Collective::getPoissonCorrection() {
    return (PoissonCorrection);
}
#endif
