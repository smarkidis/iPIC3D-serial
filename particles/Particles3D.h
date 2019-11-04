#ifndef Part2D_H
#define Part2D_H

#include "Particles3Dcomm.h"

class Particles3D:public Particles3Dcomm {
public:
  /** constructor */
  Particles3D();
  /** destructor */
  ~Particles3D();
  /** Initial condition: uniform in space and motionless */
  void uniform_background(Grid * grid, Field * EMf);
  /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
          <ul>
          <li> dim = 0 --> constant velocity on X direction </li>
          <li> dim = 1 --> constant velocity on Y direction </li>
          </ul>
          */
  void constantVelocity(double vel, int dim, Grid * grid, Field * EMf);
  /** Initial condition: uniform in space and maxwellian in velocity */
  void maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
  /** Initialize particles with a given pitch angle and energy (B field assumed to be in z direction) */
  void pitch_angle_energy(Grid * grid, Field * EMf, VirtualTopology3D * vct);
  /** Force Free initialization (JxB=0) for particles */
  void force_free(Grid * grid, Field * EMf, VirtualTopology3D * vct);
  /** Initial condition: uniform in space and maxwellian in velocity */
  void alt_maxwellian(Grid * grid, Field * EMf, VirtualTopology3D * vct);
  /** Linear_perturbation */
  void linear_perturbation(double deltaBX, double kx, double ky, double theta, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid * grid, Field * EMf, VirtualTopology3D * vct);
  /**Add a periodic perturbation in velocity exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
  void AddPerturbationJ(double deltaBoB, double kx, double ky, double kz, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid * grid);
  /** Linear delta f for bi-maxwellian plasma */
  double delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_ampl, double Ex_phase, double Ey_ampl, double Ey_phase, double Ez_ampl, double Ez_phase, double theta, Field * EMf);
  /** Derivative of f0 wrt vpar */
  double df0_dvpar(double vpar, double vperp);
  /** Derivative of f0 wrt vperp */
  double df0_dvperp(double vpar, double vperp);
  /** Equilibrium bi-maxwellian f0 */
  double f0(double vpar, double vperp);
  /** Rotate velocities in plane XY of angle theta */
  void RotatePlaneXY(double theta);
  /** mover with the esplicit non relativistic scheme */
  void mover_explicit(Grid * grid, VirtualTopology3D * vct, Field * EMf);
  /** mover with a Predictor-Corrector Scheme */
  int mover_PC(Grid * grid, VirtualTopology3D * vct, Field * EMf);
  /** relativistic mover with a Predictor-Corrector scheme */
  int mover_relativistic(Grid * grid, VirtualTopology3D * vct, Field * EMf);
  /** interpolation Particle->Grid only charge density, current */
  void interpP2G_notP(Field * EMf, Grid * grid, VirtualTopology3D * vct);
  /** interpolation Particle->Grid only for pressure tensor */
  void interpP2G_onlyP(Field * EMf, Grid * grid, VirtualTopology3D * vct);

};


#endif
