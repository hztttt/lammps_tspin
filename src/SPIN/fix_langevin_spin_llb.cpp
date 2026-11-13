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

   Spin bath of the Landau-Lifshitz-Bloch integrator fix nve/spin/llb:
   the temperature, the transverse damping of the spin directions, the
   longitudinal damping of the spin magnitudes, and the random number
   seed.  The fix holds parameters only.  fix nve/spin/llb draws the
   thermal noise and applies the damping itself, since the noise has to be
   kept fixed during the iterations of a timestep.

   The two damping constants are parameters of the atomistic model.  No
   parameter depends on the temperature: the bath samples exp(-E/kT) d^3S
   for the energy E that the pair styles and spin fixes define, and the
   temperature dependence of the magnetization and of the longitudinal
   susceptibility follow from that sampling.
------------------------------------------------------------------------- */

#include "fix_langevin_spin_llb.h"

#include "error.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixLangevinSpinLLB::FixLangevinSpinLLB(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), temp(-1.0), alpha_perp(-1.0), alpha_par(-1.0), seed(0), tnoise(1),
    lnoise(1)
{
  if (narg < 9) utils::missing_cmd_args(FLERR, "fix langevin/spin/llb", error);

  int iarg = 3;
  while (iarg < narg) {
    if (iarg + 2 > narg)
      utils::missing_cmd_args(FLERR, std::string("fix langevin/spin/llb ") + arg[iarg], error);
    if (strcmp(arg[iarg], "temp") == 0) {
      temp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (temp < 0.0) error->all(FLERR, iarg + 1, "Fix langevin/spin/llb temp must be >= 0");
    } else if (strcmp(arg[iarg], "alpha_perp") == 0) {
      alpha_perp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (alpha_perp < 0.0)
        error->all(FLERR, iarg + 1, "Fix langevin/spin/llb alpha_perp must be >= 0");
    } else if (strcmp(arg[iarg], "alpha_par") == 0) {
      alpha_par = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (alpha_par < 0.0)
        error->all(FLERR, iarg + 1, "Fix langevin/spin/llb alpha_par must be >= 0");
    } else if (strcmp(arg[iarg], "seed") == 0) {
      seed = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (seed <= 0) error->all(FLERR, iarg + 1, "Fix langevin/spin/llb seed must be > 0");
    } else if (strcmp(arg[iarg], "noise") == 0) {
      const std::string value = arg[iarg + 1];
      if (value == "both") {
        tnoise = lnoise = 1;
      } else if (value == "transverse") {
        tnoise = 1;
        lnoise = 0;
      } else if (value == "longitudinal") {
        tnoise = 0;
        lnoise = 1;
      } else if (value == "none") {
        tnoise = lnoise = 0;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix langevin/spin/llb noise option: {}", value);
      }
    } else {
      error->all(FLERR, iarg, "Unknown fix langevin/spin/llb keyword: {}", arg[iarg]);
    }
    iarg += 2;
  }

  if ((temp < 0.0) || (alpha_perp < 0.0) || (alpha_par < 0.0))
    error->all(FLERR,
               "Fix langevin/spin/llb requires the temp, alpha_perp, and alpha_par keywords");
  if ((temp > 0.0) && (tnoise || lnoise) && (seed <= 0))
    error->all(FLERR, "Fix langevin/spin/llb requires the seed keyword for thermal noise");
}

/* ---------------------------------------------------------------------- */

int FixLangevinSpinLLB::setmask()
{
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixLangevinSpinLLB::init()
{
  if (modify->get_fix_by_style("^nve/spin/llb$").empty())
    error->all(FLERR, Error::NOLASTLINE,
               "Fix langevin/spin/llb requires fix nve/spin/llb, which applies it");
}

/* ----------------------------------------------------------------------
   parameters of the bath, for fix nve/spin/llb
------------------------------------------------------------------------- */

void *FixLangevinSpinLLB::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "temp") == 0) return &temp;
  if (strcmp(str, "alpha_perp") == 0) return &alpha_perp;
  if (strcmp(str, "alpha_par") == 0) return &alpha_par;
  if (strcmp(str, "seed") == 0) return &seed;
  if (strcmp(str, "tnoise") == 0) return &tnoise;
  if (strcmp(str, "lnoise") == 0) return &lnoise;
  return nullptr;
}
