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

   Single-site Landau energy of the spin modulus,

       E_i = A |S_i|^2 + B |S_i|^4 + C |S_i|^6 ,

   the longitudinal part of a Landau-Heisenberg Hamiltonian.  With A < 0 and
   B > 0 it has a minimum at a finite modulus, with A > 0 at zero modulus.
   The energy depends only on the modulus, so it exerts no torque: it leaves
   fixed-modulus Landau-Lifshitz dynamics unchanged and acts only in styles
   in which the modulus is a dynamical variable.

   The force is stored in fm with the convention of the inertial spin
   dynamics styles (fix nve/tspin and friends) and of pair style deepspin,
   the full gradient scaled by |S|/hbar:

       fm_i += |S_i|/hbar * (-dE_i/dS_i)
             = -(2A |S_i|^2 + 4B |S_i|^4 + 6C |S_i|^6)/hbar * s_hat_i .

   The polynomial is smooth at zero modulus, so a spin passing through
   |S| = 0 needs no special treatment.
------------------------------------------------------------------------- */

#include "fix_landau_spin.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

/* ---------------------------------------------------------------------- */

FixLandauSpin::FixLandauSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), coeff_a(0.0), coeff_b(0.0), coeff_c(0.0), elandau(0.0), hbar(0.0),
    ilevel_respa(0)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix landau/spin", error);

  if (!atom->sp_flag) error->all(FLERR, "Fix landau/spin requires atom style spin");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  energy_global_flag = 1;
  respa_level_support = 1;

  coeff_a = utils::numeric(FLERR, arg[3], false, lmp);
  coeff_b = utils::numeric(FLERR, arg[4], false, lmp);

  int iarg = 5;
  int iarg_c = Error::NOPOINTER;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "c") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix landau/spin c", error);
      coeff_c = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg_c = iarg + 1;
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix landau/spin keyword: {}", arg[iarg]);
    }
  }

  // the energy must be bounded from below, otherwise the modulus runs away:
  // the highest-order nonzero coefficient has to be positive

  if (coeff_c < 0.0) error->all(FLERR, iarg_c, "Fix landau/spin coefficient C must be >= 0");
  if ((coeff_c == 0.0) && (coeff_b < 0.0))
    error->all(FLERR, 4, "Fix landau/spin coefficient B must be >= 0 when C is zero");
  if ((coeff_c == 0.0) && (coeff_b == 0.0) && (coeff_a <= 0.0))
    error->all(FLERR, 3, "Fix landau/spin coefficient A must be > 0 when B and C are zero");

  hbar = force->hplanck / MY_2PI;
}

/* ---------------------------------------------------------------------- */

int FixLandauSpin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::init()
{
  hbar = force->hplanck / MY_2PI;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    auto *respa = dynamic_cast<Respa *>(update->integrate);
    respa->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    respa->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   Landau energy of a site with spin modulus smag
------------------------------------------------------------------------- */

double FixLandauSpin::site_energy(double smag) const
{
  const double s2 = smag * smag;
  return s2 * (coeff_a + s2 * (coeff_b + s2 * coeff_c));
}

/* ----------------------------------------------------------------------
   longitudinal drive -(1/hbar) dE/d|S| of a site with spin modulus smag
------------------------------------------------------------------------- */

double FixLandauSpin::site_drive(double smag) const
{
  const double s2 = smag * smag;
  return -2.0 * smag * (coeff_a + s2 * (2.0 * coeff_b + 3.0 * s2 * coeff_c)) / hbar;
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::post_force(int /*vflag*/)
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  elandau = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double smag = sp[i][3];
    elandau += site_energy(smag);

    // fm = |S|/hbar * (-dE/dS) = |S| * site_drive(|S|) * s_hat; purely radial

    const double pref = smag * site_drive(smag);
    fm[i][0] += pref * sp[i][0];
    fm[i][1] += pref * sp[i][1];
    fm[i][2] += pref * sp[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLandauSpin::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   Landau energy of the group
------------------------------------------------------------------------- */

double FixLandauSpin::compute_scalar()
{
  double all = 0.0;
  MPI_Allreduce(&elandau, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}
