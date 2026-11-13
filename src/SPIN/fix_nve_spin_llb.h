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
FixStyle(nve/spin/llb,FixNVESpinLLB);
// clang-format on
#else

#ifndef LMP_FIX_NVE_SPIN_LLB_H
#define LMP_FIX_NVE_SPIN_LLB_H

#include "fix_nve_spin_sync.h"

#include <vector>

namespace LAMMPS_NS {

class FixNVESpinLLB : public FixNVESpinSync {
 public:
  FixNVESpinLLB(class LAMMPS *, int, char **);

  void init() override;

 protected:
  enum { SIB, PC };
  enum { PAIR_NONE, PAIR_PERSPIN, PAIR_GLOBAL };

  int scheme;             // SIB or PC
  double alpha_par;       // longitudinal damping
  double kt_long;         // kB T of the longitudinal noise, 0 without it
  int user_alpha_perp;    // 1 if alpha_perp was given on this fix
  int user_alpha_par;     // 1 if alpha_par was given on this fix
  double alpha_perp_user, alpha_par_user;
  int pair_long;    // how the pair styles contribute to the longitudinal drive

  std::vector<class SpinFieldFix *> modulus_fixes;
  std::vector<class PairSpin *> full_pairs;    // spin pair styles with a longitudinal force

  std::vector<double> u0, u1, k1, lnoise;
  std::vector<int> active;

  int parse_extra(int, int, char **) override;
  void init_bath() override;
  void advance_spins() override;
  void rebuild_field(int, double *) override;

  double longitudinal_drive(int);
  double longitudinal_rate(int);
  void refresh_longitudinal();
  void start_longitudinal(double);
  void longitudinal_sib(double);
  void longitudinal_pc(double);
};

}    // namespace LAMMPS_NS

#endif
#endif
