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

#ifdef FIX_CLASS
// clang-format off
FixStyle(mc/spin,FixMCSpin);
// clang-format on
#else

#ifndef LMP_FIX_MC_SPIN_H
#define LMP_FIX_MC_SPIN_H

#include "fix.h"

#include <vector>

namespace LAMMPS_NS {

class FixMCSpin : public Fix {
 public:
  FixMCSpin(class LAMMPS *, int, char **);
  ~FixMCSpin() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void setup_pre_neighbor() override;
  void pre_neighbor() override;
  void pre_force(int) override;
  double compute_vector(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void *extract(const char *, int &) override;

 protected:
  enum { EXCHANGE, BIQUADRATIC, NEEL };

  double temperature, beta;
  int seed, nsweep, every;
  double cone;           // half-angle of the cone of trial directions (radians)
  double cone_target;    // acceptance rate the cone adaptation aims at
  int adapt_flag;
  int modulus_flag;     // 1 if trial moves also change the spin magnitude
  int jacobian_flag;    // 1 if the |S|^2 weight of the magnitude is included
  double mod_sigma, mod_min, mod_max;
  int energy_global;    // 1 if every trial is scored by re-running the pair style
  int lattice_flag;     // -1: the atoms are moved, or not, by other fixes

  double naccept, nattempt, mod_naccept, mod_nattempt;
  double accept_rate, mod_accept_rate;
  double epair_cur;    // pair energy of the current configuration (energy global)

  std::vector<class Pair *> spin_pairs;
  std::vector<int> pair_kind;
  std::vector<std::vector<double>> pair_cutsq;    // per pair of types, < 0 where not applied
  std::vector<class FixPrecessionSpin *> precession_fixes;
  std::vector<class SpinFieldFix *> modulus_fixes;

  class NeighList *list;
  class RanMars *random;    // with energy global, the same numbers on all processors
  std::vector<int> order;
  std::vector<tagint> trial_tags;    // atom IDs in the order of a sweep (energy global)

  // periodic images of the local atoms on this processor (energy local, serial)

  int image_flag;
  std::vector<int> image_head, image_next;

  // sectoring of the subdomains for parallel sweeps, as in fix nve/spin

  int sector_flag, nsectors, nlocal_max;
  double rsec[3];
  int *stack_head, *stack_foot, *forward_stacks, *backward_stacks;

  double bond_energy(int, const double *, int);
  double single_spin_energy(int, const double *, double);
  double single_site_energy(int, const double *);
  double magnetic_energy();
  double pair_energy_local();
  double model_energy();
  void clear_forces();
  void set_spin(int, const double *, double);
  void propose_direction(const double *, double *);
  void mc_single_spin(int);
  void mc_global_spin(int);
  void local_sweep();
  void global_sweep();
  void shuffle_and_try(int);
  void link_images();
  void mc_diagnostics();
  void verify_energy_local();
  void verify_energy_global();
  void sectoring();
  int coords2sector(double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
