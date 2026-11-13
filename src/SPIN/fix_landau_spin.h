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
FixStyle(landau/spin,FixLandauSpin);
// clang-format on
#else

#ifndef LMP_FIX_LANDAU_SPIN_H
#define LMP_FIX_LANDAU_SPIN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLandauSpin : public Fix {
 public:
  FixLandauSpin(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  double compute_scalar() override;

  // Landau energy of one site with spin modulus smag, in energy units
  double site_energy(double smag) const;
  // -(1/hbar) dE/d|S| of one site with spin modulus smag
  double site_drive(double smag) const;

 protected:
  double coeff_a, coeff_b, coeff_c;    // coefficients of |S|^2, |S|^4, |S|^6
  double elandau;                      // Landau energy of the local atoms
  double hbar;                         // reduced Planck constant (energy*time/rad)
  int ilevel_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
