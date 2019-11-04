#include <iostream>
#include <math.h>
#include "../processtopology/VirtualTopology3D.h"
#include "../processtopology/VCtopology3D.h"
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../grids/Grid.h"
#include "../grids/Grid3DCU.h"
#include "../fields/Field.h"
#include "Particles3Dcomm.h"

#include <vector>
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));



/** constructor */
Particles3Dcomm::Particles3Dcomm() {
  // see allocate(int species, CollectiveIO* col, VirtualTopology3D* vct, Grid* grid)

}
/** deallocate particles */
Particles3Dcomm::~Particles3Dcomm() {
  delete[]x;
  delete[]y;
  delete[]z;
  delete[]u;
  delete[]v;
  delete[]w;
  delete[]q;
}
/** constructors fo a single species*/
void Particles3Dcomm::allocate(int species, CollectiveIO * col, VirtualTopology3D * vct, Grid * grid) {
  
  // subcycling is set here
  
    
  // info from collectiveIO
  ns = species;
  npcel = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  npcelz = col->getNpcelz(species);
  nop = col->getNp(species) / (vct->getNprocs());
  np_tot = col->getNp(species);
  npmax = col->getNpMax(species) / (vct->getNprocs());
  qom = col->getQOM(species);
  uth = col->getUth(species);
  vth = col->getVth(species);
  wth = col->getWth(species);
  u0 = col->getU0(species);
  v0 = col->getV0(species);
  w0 = col->getW0(species);
  dt = col->getDt();
  Lx = col->getLx();
  Ly = col->getLy();
  Lz = col->getLz();
  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();
  delta = col->getDelta();
  TrackParticleID = col->getTrackParticleID(species);
  c = col->getC();
  // info for mover
  NiterMover = col->getNiterMover();
  if (qom > 0)
      NiterMover = 1;
  // number of subcycling
  n_sub_cycles = col->getN_sub_cycles();
  if (qom > 0)
    n_sub_cycles = 1;
  // velocity of the injection from the wall
  Vinj = col->getVinj();
  // info from Grid
  xstart = grid->getXstart();
  xend = grid->getXend();
  ystart = grid->getYstart();
  yend = grid->getYend();
  zstart = grid->getZstart();
  zend = grid->getZend();

  dx = grid->getDX();
  dy = grid->getDY();
  dz = grid->getDZ();

  nxn = grid->getNXN();
  nyn = grid->getNYN();
  nzn = grid->getNZN();
  invVOL = grid->getInvVOL();
  // info from VirtualTopology3D
  cVERBOSE = vct->getcVERBOSE();

  // boundary condition for particles
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
  bcPfaceZright = col->getBcPfaceZright();
  bcPfaceZleft = col->getBcPfaceZleft();
  // //////////////////////////////////////////////////////////////
  // ////////////// ALLOCATE ARRAYS /////////////////////////
  // //////////////////////////////////////////////////////////////
  // positions
  x = new double[npmax];
  y = new double[npmax];
  z = new double[npmax];
  // velocities
  u = new double[npmax];
  v = new double[npmax];
  w = new double[npmax];
  // charge
  q = new double[npmax];
  // ID
  if (TrackParticleID)
    ParticleID = new unsigned long[npmax];
  // BUFFERS
  // the buffer size should be decided depending on number of particles
  // the buffer size should be decided depending on number of particles
  if (TrackParticleID)
    nVar = 8;
  else
    nVar = 7;

  // if RESTART is true initialize the particle in allocate method
  restart = col->getRestart_status();

}
/** calculate the weights given the position of particles 0,0,0 is the left,left, left node */
void Particles3Dcomm::calculateWeights(double weight[][2][2], double xp, double yp, double zp, int ix, int iy, int iz, Grid * grid) {
  double xi[2], eta[2], zeta[2];
  xi[0] = xp - grid->getXN(ix - 1, iy, iz);
  eta[0] = yp - grid->getYN(ix, iy - 1, iz);
  zeta[0] = zp - grid->getZN(ix, iy, iz - 1);
  xi[1] = grid->getXN(ix, iy, iz) - xp;
  eta[1] = grid->getYN(ix, iy, iz) - yp;
  zeta[1] = grid->getZN(ix, iy, iz) - zp;
  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      for (int k = 0; k < 2; k++)
        weight[i][j][k] = xi[i] * eta[j] * zeta[k] * invVOL;
}

/** Interpolation Particle --> Grid */
void Particles3Dcomm::interpP2G(Field * EMf, Grid * grid, VirtualTopology3D * vct) {
  double weight[2][2][2];
  double temp[2][2][2];
  double xi[2], eta[2], zeta[2];
  int ix, iy, iz;
  double inv_dx = 1.0 / dx, inv_dy = 1.0 / dy, inv_dz = 1.0 / dz;
  for (register long long i = 0; i < nop; i++) {
    ix = 2 + int (floor((x[i] - xstart) * inv_dx));
    iy = 2 + int (floor((y[i] - ystart) * inv_dy));
    iz = 2 + int (floor((z[i] - zstart) * inv_dz));
    xi[0] = x[i] - grid->getXN(ix - 1, iy, iz);
    eta[0] = y[i] - grid->getYN(ix, iy - 1, iz);
    zeta[0] = z[i] - grid->getZN(ix, iy, iz - 1);
    xi[1] = grid->getXN(ix, iy, iz) - x[i];
    eta[1] = grid->getYN(ix, iy, iz) - y[i];
    zeta[1] = grid->getZN(ix, iy, iz) - z[i];
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          weight[ii][jj][kk] = q[i] * xi[ii] * eta[jj] * zeta[kk] * invVOL;
    // add charge density
    EMf->addRho(weight, ix, iy, iz, ns);
    // add current density - X
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = u[i] * weight[ii][jj][kk];
    EMf->addJx(temp, ix, iy, iz, ns);
    // add current density - Y
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = v[i] * weight[ii][jj][kk];
    EMf->addJy(temp, ix, iy, iz, ns);
    // add current density - Z
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = w[i] * weight[ii][jj][kk];
    EMf->addJz(temp, ix, iy, iz, ns);
    // Pxx - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = u[i] * u[i] * weight[ii][jj][kk];
    EMf->addPxx(temp, ix, iy, iz, ns);
    // Pxy - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = u[i] * v[i] * weight[ii][jj][kk];
    EMf->addPxy(temp, ix, iy, iz, ns);
    // Pxz - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = u[i] * w[i] * weight[ii][jj][kk];
    EMf->addPxz(temp, ix, iy, iz, ns);
    // Pyy - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = v[i] * v[i] * weight[ii][jj][kk];
    EMf->addPyy(temp, ix, iy, iz, ns);
    // Pyz - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = v[i] * w[i] * weight[ii][jj][kk];
    EMf->addPyz(temp, ix, iy, iz, ns);
    // Pzz - add pressure tensor
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++)
        for (int kk = 0; kk < 2; kk++)
          temp[ii][jj][kk] = w[i] * w[i] * weight[ii][jj][kk];
    EMf->addPzz(temp, ix, iy, iz, ns);
  }
  // communicate contribution from ghost cells 
  EMf->communicateGhostP2G(ns, 0, 0, 0, 0, vct);
}


/** communicate buffers */
int Particles3Dcomm::communicate(VirtualTopology3D * ptVCT) {
  
  return (0);                   // everything was fine


}
/** resize the buffers */
void Particles3Dcomm::resize_buffers(int new_buffer_size) {
  
}
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
void Particles3Dcomm::bufferXright(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferYright(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** put a particle exiting to Z-LEFT in the bufferZLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZleft(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** put a particle exiting to Z-RIGHT in the bufferZRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles3Dcomm::bufferZright(double *b_, long long np_current, VirtualTopology3D * vct) {
  
}
/** This unbuffer the last communication */
int Particles3Dcomm::unbuffer(double *b_) {
    return (0);                   // everything was fine
}
/** Delete the a particle from the array and pack the the array, update the number of 
 * particles that are exiting
 * For deleting the particle from the array take the last particle and put it
 * in the position of the particle you want to delete
 * @param np = the index of the particle that must be deleted
 * @param nplast = the index of the last particle in the array
 */
void Particles3Dcomm::del_pack(long long np_current, long long *nplast) {
  x[np_current] = x[*nplast];
  y[np_current] = y[*nplast];
  z[np_current] = z[*nplast];
  u[np_current] = u[*nplast];
  v[np_current] = v[*nplast];
  w[np_current] = w[*nplast];
  q[np_current] = q[*nplast];
  if (TrackParticleID)
    ParticleID[np_current] = ParticleID[*nplast];
  npExit++;
  (*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles3Dcomm::isMessagingDone(VirtualTopology3D * ptVCT) {
    return(1);

}
/** calculate the maximum number exiting from this domain */
int Particles3Dcomm::maxNpExiting() {
  int maxNp = 0;
  if (npExitXright > maxNp)
    maxNp = npExitXright;
  if (npExitXleft > maxNp)
    maxNp = npExitXleft;
  if (npExitYright > maxNp)
    maxNp = npExitYright;
  if (npExitYleft > maxNp)
    maxNp = npExitYleft;
  if (npExitZright > maxNp)
    maxNp = npExitZright;
  if (npExitZleft > maxNp)
    maxNp = npExitZleft;
  return (maxNp);
}
/** return X-coordinate of particle array */
double *Particles3Dcomm::getXall() const {
  return (x);
}
/** return Y-coordinate  of particle array */
double *Particles3Dcomm::getYall() const {
  return (y);
}
/** return Z-coordinate  of particle array*/
double *Particles3Dcomm::getZall() const {
  return (z);
}
/** get X-velocity of particle with label indexPart */
double *Particles3Dcomm::getUall() const {
  return (u);
}
/** get Y-velocity of particle with label indexPart */
double *Particles3Dcomm::getVall() const {
  return (v);
}
/**get Z-velocity of particle with label indexPart */
double *Particles3Dcomm::getWall() const {
  return (w);
}
/**get ID of particle with label indexPart */
unsigned long *Particles3Dcomm::getParticleIDall() const {
  return (ParticleID);
}
/**get charge of particle with label indexPart */
double *Particles3Dcomm::getQall() const {
  return (q);
}
/** return X-coordinate of particle with index indexPart */
double Particles3Dcomm::getX(long long indexPart) const {
  return (x[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getY(long long indexPart) const {
  return (y[indexPart]);
}
/** return Y-coordinate  of particle with index indexPart */
double Particles3Dcomm::getZ(long long indexPart) const {
  return (z[indexPart]);
}
/** get u (X-velocity) of particle with label indexPart */
double Particles3Dcomm::getU(long long indexPart) const {
  return (u[indexPart]);
}
/** get v (Y-velocity) of particle with label indexPart */
double Particles3Dcomm::getV(long long indexPart) const {
  return (v[indexPart]);
}
/**get w (Z-velocity) of particle with label indexPart */
double Particles3Dcomm::getW(long long indexPart) const {
  return (w[indexPart]);
}
/**get ID of particle with label indexPart */
unsigned long Particles3Dcomm::getParticleID(long long indexPart) const {
  return (ParticleID[indexPart]);
}
/**get charge of particle with label indexPart */
double Particles3Dcomm::getQ(long long indexPart) const {
  return (q[indexPart]);
}
/** return the number of particles */
long long Particles3Dcomm::getNOP() const {
  return (nop);
}
/** return the Kinetic energy */
double Particles3Dcomm::getKe() {
  double localKe = 0.0;
  double totalKe = 0.0;
  for (register long long i = 0; i < nop; i++)
    localKe += .5 * (q[i] / qom) * (u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    totalKe = localKe;
  return (totalKe);
}
/** return the total momentum */
double Particles3Dcomm::getP() {
  double localP = 0.0;
  double totalP = 0.0;
  for (register long long i = 0; i < nop; i++)
    localP += (q[i] / qom) * sqrt(u[i] * u[i] + v[i] * v[i] + w[i] * w[i]);
    totalP = localP;
  return (totalP);
}

/** return the highest kinetic energy */
double Particles3Dcomm::getMaxVelocity(){  
  double localVel = 0.0;
  double maxVel = 0.0;
  for (long long i=0; i < nop; i++)
    localVel = max(localVel, sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]));
    maxVel = localVel;
  return(maxVel);
}


/** get energy spectrum */
unsigned long* Particles3Dcomm::getVelocityDistribution(int nBins, double maxVel){  
  unsigned long* f = new unsigned long [nBins];
  for (int i=0; i < nBins; i++)
      f[i] = 0;
  double Vel = 0.0;
  double dv = maxVel / nBins;
  int bin = 0;
  for (long long i=0; i < nop; i++) {
    Vel = sqrt(u[i]*u[i] + v[i]*v[i] + w[i]*w[i]);
    bin = int(floor(Vel/dv));
    if (bin >= nBins) 
      f[nBins-1] += 1;
    else
      f[bin] += 1;
  }
  unsigned long localN = 0;
  unsigned long totalN = 0;
  for (int i=0; i < nBins; i++) {
    localN = f[i];
    totalN = localN;
    f[i] = totalN;
  }
  return f;
}


/** print particles info */
void Particles3Dcomm::Print(VirtualTopology3D * ptVCT) const {
  cout << endl;
  cout << "Number of Particles: " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << "Xin = " << xstart << "; Xfin = " << xend << endl;
  cout << "Yin = " << ystart << "; Yfin = " << yend << endl;
  cout << "Zin = " << zstart << "; Zfin = " << zend << endl;
  cout << "Number of species = " << ns << endl;
  for (long long i = 0; i < nop; i++)
    cout << "Particles #" << i << " x=" << x[i] << " y=" << y[i] << " z=" << z[i] << " u=" << u[i] << " v=" << v[i] << " w=" << w[i] << endl;
  cout << endl;
}
/** print just the number of particles */
void Particles3Dcomm::PrintNp(VirtualTopology3D * ptVCT) const {
  cout << endl;
  cout << "Number of Particles of species " << ns << ": " << nop << endl;
  cout << "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << "," << ptVCT->getCoordinates(2) << ")" << endl;
  cout << endl;
}
