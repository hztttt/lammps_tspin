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

   Interface of the fixes that apply an energy to the spin magnitudes,
   E_i(|S_i|), such as fix landau/spin and fix spring/tspin.

   In their post_force() these fixes add the radial magnetic force
   fm_i += |S_i| * L_i * s_hat_i, with the longitudinal drive

       L_i = (1/hbar) (-dE_i / d|S_i|) ,

   following the convention of the inertial spin dynamics styles.  A spin
   integrator that evolves the spin magnitudes evaluates L_i at
   intermediate configurations of a timestep through this interface, since
   it cannot call post_force() again: that would apply the forces twice.
------------------------------------------------------------------------- */

#ifndef LMP_SPIN_FIELD_FIX_H
#define LMP_SPIN_FIELD_FIX_H

namespace LAMMPS_NS {

class SpinFieldFix {
 public:
  // longitudinal drive (1/hbar)(-dE/d|S|) on local atom i at its current spin
  // magnitude, zero for atoms outside the group of the fix

  virtual double longitudinal_drive(int i) = 0;

 protected:
  // the fixes are never deleted through this interface; a virtual destructor
  // would conflict with the exception specification of the Fix destructor

  ~SpinFieldFix() = default;
};

}    // namespace LAMMPS_NS

#endif
