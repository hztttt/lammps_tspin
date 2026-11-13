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
FixStyle(langevin/spin/llb,FixLangevinSpinLLB);
// clang-format on
#else

#ifndef LMP_FIX_LANGEVIN_SPIN_LLB_H
#define LMP_FIX_LANGEVIN_SPIN_LLB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLangevinSpinLLB : public Fix {
 public:
  FixLangevinSpinLLB(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void *extract(const char *, int &) override;

 protected:
  double temp;          // temperature of the bath
  double alpha_perp;    // transverse (Gilbert) damping
  double alpha_par;     // longitudinal damping
  int seed;
  int tnoise, lnoise;    // 1 if the transverse / longitudinal thermal noise is applied
};

}    // namespace LAMMPS_NS

#endif
#endif
