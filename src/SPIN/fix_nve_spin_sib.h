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
FixStyle(nve/spin/sib,FixNVESpinSIB);
// clang-format on
#else

#ifndef LMP_FIX_NVE_SPIN_SIB_H
#define LMP_FIX_NVE_SPIN_SIB_H

#include "fix_nve_spin_sync.h"

namespace LAMMPS_NS {

class FixNVESpinSIB : public FixNVESpinSync {
 public:
  FixNVESpinSIB(class LAMMPS *, int, char **);

 protected:
  int parse_extra(int, int, char **) override;
  void advance_spins() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
