/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NVE_SPIN_SYNC_H
#define LMP_FIX_NVE_SPIN_SYNC_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
   common base of the spin integrators that advance all spins of the group
   together with the effective field of one common configuration
   (nve/spin/sib, nve/spin/pc), in contrast to the sequential sweep of
   fix nve/spin
------------------------------------------------------------------------- */

class FixNVESpinSync : public Fix {
 public:
  FixNVESpinSync(class LAMMPS *, int, char **);
  ~FixNVESpinSync() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double compute_vector(int) override;
  void *extract(const char *, int &) override;

 protected:
  enum { FIELD_AUTO, FIELD_PERSPIN, FIELD_GLOBAL };

  int lattice_flag;         // 1 if atom positions are integrated
  int field_mode;           // choice of the field keyword
  int field_global;         // 1 if the pair field is evaluated by re-running the pair style
  int field_tol_loose;      // 1 if the pair style may not reproduce its field exactly
  int check_parallel;       // 1 if the field check includes the components along the spins
  int field_checked;        // 1 once the field of the first step has been verified
  double field_check[3];    // largest deviation, field scale, re-run spread of the check
  double dtv, dtf;          // timesteps for positions and velocities
  double hbar;              // reduced Planck constant (energy*time/rad)

  // implicit midpoint iteration and its diagnostics

  int maxiter;              // largest number of field evaluations of one midpoint solve
  double tolerance;         // convergence criterion of a midpoint solve
  bigint nsteps, nevals;    // timesteps and field evaluations in this run
  bigint nsolve, nfail;     // midpoint solves in this run, and those without convergence
  double residual_max;      // largest final change of a midpoint solve in this run

  // spin bath: transverse damping and thermal field

  int bathbit;        // group bit of the bath, 0 if there is none
  double alpha;       // transverse (Gilbert) damping
  double temp;        // temperature of the bath
  double sigma;       // standard deviation of the thermal field components
  int seed_offset;    // added to the seed of the bath for the random numbers
  int random_seed;    // seed the random number generator was created with
  class RanMars *random;
  std::vector<double> restart_state;    // random number states read from a restart file

  std::vector<class PairSpin *> spin_pairs;
  std::vector<class FixPrecessionSpin *> precession_fixes;
  class FixSetForceSpin *setforce_fix;

  std::vector<double> fsave;                       // atom forces kept across a pair re-run
  std::vector<double> sold, swork, snew, noise;    // work arrays of the transverse step

  void parse_args(int, int, char **);
  virtual int parse_extra(int, int, char **) { return 0; }
  virtual void init_bath();
  void create_random(int);

  // advance the spins of the group by one timestep

  virtual void advance_spins() = 0;
  void transverse_sib();
  void transverse_pc();

  void refresh_field();
  void spin_field(int, double *);
  virtual void rebuild_field(int, double *);
  void check_field();
  void draw_noise(std::vector<double> &, int, double);
  double damping(int) const;
  void evaluate_pair();

  static void llg_omega(const double *, const double *, double, double *);
  static void llg_rate(const double *, const double *, double, double *);
  static void cayley(const double *, const double *, double, double *);
  static void normalize(double *);
};

}    // namespace LAMMPS_NS

#endif
