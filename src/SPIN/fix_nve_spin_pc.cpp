// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
// clang-format on

/* ------------------------------------------------------------------------
   Contributing author: Zhengtao Huang (The University of Hong Kong)
                        hzt990224@gmail.com

   Predictor-corrector (Heun) integration of the spin directions, applied
   to all spins of the group together:

       k1 = ds/dt( s(t), H(s(t)) ),        s_p = normalized( s(t) + dt k1 )
       k2 = ds/dt( s_p,  H(s_p) ),
       s(t+dt) = normalized( s(t) + dt (k1 + k2) / 2 )

   with two evaluations of the field of all spins per step.  With a thermal
   field held fixed over the step this is the Heun scheme for the
   stochastic equation in the Stratonovich interpretation.  The scheme is
   explicit, so it does not conserve the energy of an undamped system.
------------------------------------------------------------------------- */

#include "fix_nve_spin_pc.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixNVESpinPC::FixNVESpinPC(LAMMPS *lmp, int narg, char **arg) : FixNVESpinSync(lmp, narg, arg)
{
  parse_args(3, narg, arg);
}

/* ---------------------------------------------------------------------- */

void FixNVESpinPC::advance_spins()
{
  transverse_pc();
}
