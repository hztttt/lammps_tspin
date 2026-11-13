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

   Self-consistent implicit midpoint integration of the spin directions
   (semi-implicit B scheme, SIB, of Mentink et al., J. Phys.: Condens.
   Matter 22, 176001 (2010)), applied to all spins of the group together:

       s(t+dt) = Cayley[ omega( s_mid, H(s_mid) ), dt ] s(t),
       s_mid   = (s(t) + s(t+dt)) / |s(t) + s(t+dt)|,

   solved by fixed-point iteration, one evaluation of the field of all spins
   per iteration.  The Cayley rotation keeps every spin length exactly, and
   the converged step is symmetric in time.
------------------------------------------------------------------------- */

#include "fix_nve_spin_sib.h"

#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixNVESpinSIB::FixNVESpinSIB(LAMMPS *lmp, int narg, char **arg) : FixNVESpinSync(lmp, narg, arg)
{
  size_vector = 6;
  parse_args(3, narg, arg);
}

/* ---------------------------------------------------------------------- */

int FixNVESpinSIB::parse_extra(int iarg, int narg, char **arg)
{
  if (strcmp(arg[iarg], "maxiter") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/sib maxiter", error);
    maxiter = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if (maxiter < 1) error->all(FLERR, iarg + 1, "Fix nve/spin/sib maxiter must be > 0");
    return 2;
  } else if (strcmp(arg[iarg], "tol") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/sib tol", error);
    tolerance = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    if (tolerance < 0.0) error->all(FLERR, iarg + 1, "Fix nve/spin/sib tol must be >= 0");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixNVESpinSIB::advance_spins()
{
  transverse_sib();
}
