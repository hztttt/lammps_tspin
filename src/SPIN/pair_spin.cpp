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

/* ------------------------------------------------------------------------
   Contributing authors: Julien Tranchida (SNL)
                         Aidan Thompson (SNL)

   Please cite the related publication:
   Tranchida, J., Plimpton, S. J., Thibaudeau, P., & Thompson, A. P. (2018).
   Massively parallel symplectic algorithm for coupled magnetic spin dynamics
   and molecular dynamics. Journal of Computational Physics.
------------------------------------------------------------------------- */

#include "pair_spin.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairSpin::PairSpin(LAMMPS *lmp) : Pair(lmp), emag(nullptr)
{
  hbar = force->hplanck/MY_2PI;
  single_enable = 0;
  respa_enable = 0;
  no_virial_fdotr_compute = 1;
  lattice_flag = 0;
  full_spin_gradient = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSpin::settings(int narg, char **/*arg*/)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect number of args in pair_style pair/spin command");

  // pair spin/* need the metal unit style

  if (strcmp(update->unit_style,"metal") != 0)
    error->all(FLERR,"Spin pair styles require metal units");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSpin::init_style()
{
  if (!atom->sp_flag)
    error->all(FLERR,"Pair spin requires atom/spin style");

  const int ntspin = modify->get_fix_by_style("^nve/tspin").size() +
      modify->get_fix_by_style("^nvt/tspin").size() + modify->get_fix_by_style("^npt/tspin").size() +
      modify->get_fix_by_style("^nph/tspin").size();
  const auto spin_fixes = modify->get_fix_by_style("^nve/spin");

  if (!spin_fixes.empty() && (ntspin > 0))
    error->all(FLERR, Error::NOLASTLINE,
               "Cannot combine fix {} with tspin integrators when using SPIN pair styles",
               spin_fixes.front()->style);

  // the tspin integrators need the derivative of the energy with respect to
  // the full spin vector, magnitude as well as direction

  if ((ntspin > 0) && !full_spin_gradient)
    error->all(FLERR, Error::NOLASTLINE,
               "Pair style {} does not provide the derivative of the energy with respect "
               "to the full spin vector and cannot be used with the tspin integrators",
               force->pair_style);

  // the spin integrators report through Fix::extract() whether they move the
  // atoms (1) or not (0); the mechanical forces of the magnetic interactions
  // are computed when any of them does.  A spin sampler that leaves the atoms
  // to other fixes, fix mc/spin, reports -1: the forces are then computed
  // when another fix integrates the atom positions.

  int nintegrators = 0;
  int moving = 0, frozen = 0, deferred = 0, other = 0;
  for (const auto &ifix : modify->get_fix_list()) {
    int dim = -1;
    auto *flag = static_cast<int *>(ifix->extract("lattice_flag", dim));
    if (flag && (dim == 0)) {
      nintegrators++;
      if (*flag > 0)
        moving = 1;
      else if (*flag == 0)
        frozen = 1;
      else
        deferred = 1;
    } else if (ifix->time_integrate) {
      other = 1;
    }
  }
  lattice_flag = (moving || (deferred && !frozen && other)) ? 1 : 0;

  // checking if a spin integrator or neb/spin is a listed fix

  if ((comm->me == 0) && ((nintegrators + modify->get_fix_by_style("^neb/spin").size()) == 0))
    error->warning(FLERR,"Using spin pair style without a spin integrator or neb/spin");
  if ((comm->me == 0) && (nintegrators > 1))
    error->warning(FLERR,"Using multiple spin integrators");

  // check if newton pair is on

  if (force->newton_pair == 0)
    error->all(FLERR, Error::NOLASTLINE, "Pair style spin requires newton pair on");

  // need a full neighbor list

  neighbor->add_request(this, NeighConst::REQ_FULL);

  // init. size of energy stacking lists

  nlocal_max = atom->nlocal;
  memory->grow(emag,nlocal_max,"pair/spin:emag");
}

/* ---------------------------------------------------------------------- */

double PairSpin::memory_usage()
{
  double bytes = Pair::memory_usage();
  bytes += (double) nlocal_max * sizeof(double);    // emag[nlocal_max]
  return bytes;
}
