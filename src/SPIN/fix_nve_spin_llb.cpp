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

   Stochastic Landau-Lifshitz-Bloch (LLB) dynamics of atomic spins whose
   magnitude is a dynamical variable, S_i = R_i s_i with |s_i| = 1, for the
   energy E({S_i}) that the pair styles and spin fixes define.  At
   temperature T the spins sample

       P({S_i}) prod_i d^3S_i  ~  exp(-E / kT) prod_i R_i^2 dR_i dOmega_i .

   No parameter depends on the temperature; the temperature dependence of
   the magnetization and of the longitudinal susceptibility follows from the
   sampling.

   Transverse channel: the Landau-Lifshitz-Gilbert equation of the spin
   direction on the sphere of the current radius, as integrated by fix
   nve/spin/sib and fix nve/spin/pc, with damping alpha_perp and thermal
   field variance 2 alpha_perp kT / (hbar dt).

   Longitudinal channel: Langevin dynamics of u = ln R,

       du = alpha_par [ R L + 3 kT/hbar ] dt + sqrt(2 alpha_par kT dt/hbar) xi,

   with the longitudinal drive L = (1/hbar)(-dE/dR) at fixed direction.
   Its stationary density exp(3u - E/kT) du equals R^2 exp(-E/kT) dR, the
   3 kT/hbar term being the measure of d^3S.  R = exp(u) stays positive, and
   the noise is additive, so there is no ambiguity of interpretation.  The
   measure term belongs to the noise and is omitted without it.  With
   alpha_par = 0 the magnitudes are fixed and the fix integrates the same
   equation as fix nve/spin/sib or nve/spin/pc.

   Each timestep is split symmetrically: half a longitudinal step, a full
   transverse step, and half a longitudinal step.  With scheme sib both
   channels use the implicit midpoint rule, solved by iteration with the
   thermal noise held fixed; with scheme pc both use the predictor-corrector
   (Heun) method.

   The longitudinal drive is collected from the contributions whose energy
   depends on the spin magnitudes: spin pair styles that provide the
   derivative with respect to the full spin vector (fm.s/|S| of their
   field), a pair style re-run as a whole (keyword field global, for
   example a machine-learning potential, whose fm is that full derivative),
   the Zeeman term of fix precession/spin, and fixes implementing
   SpinFieldFix such as fix landau/spin.  The energies of the other spin
   pair styles depend on the spin directions only, so their longitudinal
   drive is zero.
------------------------------------------------------------------------- */

#include "fix_nve_spin_llb.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix_precession_spin.h"
#include "force.h"
#include "modify.h"
#include "pair_spin.h"
#include "random_mars.h"
#include "spin_field_fix.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixNVESpinLLB::FixNVESpinLLB(LAMMPS *lmp, int narg, char **arg) :
    FixNVESpinSync(lmp, narg, arg), scheme(SIB), alpha_par(0.0), kt_long(0.0), user_alpha_perp(0),
    user_alpha_par(0), alpha_perp_user(0.0), alpha_par_user(0.0), pair_long(PAIR_NONE)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix nve/spin/llb", error);

  if (strcmp(arg[3], "sib") == 0)
    scheme = SIB;
  else if (strcmp(arg[3], "pc") == 0)
    scheme = PC;
  else
    error->all(FLERR, 3, "Unknown fix nve/spin/llb scheme: {}", arg[3]);

  // the random numbers of the bath are drawn here only, so no offset is needed

  seed_offset = 0;
  check_parallel = 1;
  size_vector = 6;

  parse_args(4, narg, arg);
}

/* ---------------------------------------------------------------------- */

int FixNVESpinLLB::parse_extra(int iarg, int narg, char **arg)
{
  if (strcmp(arg[iarg], "alpha_perp") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/llb alpha_perp", error);
    alpha_perp_user = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    if (alpha_perp_user < 0.0)
      error->all(FLERR, iarg + 1, "Fix nve/spin/llb alpha_perp must be >= 0");
    user_alpha_perp = 1;
    return 2;
  } else if (strcmp(arg[iarg], "alpha_par") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/llb alpha_par", error);
    alpha_par_user = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    if (alpha_par_user < 0.0)
      error->all(FLERR, iarg + 1, "Fix nve/spin/llb alpha_par must be >= 0");
    user_alpha_par = 1;
    return 2;
  } else if (strcmp(arg[iarg], "maxiter") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/llb maxiter", error);
    maxiter = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
    if (maxiter < 1) error->all(FLERR, iarg + 1, "Fix nve/spin/llb maxiter must be > 0");
    return 2;
  } else if (strcmp(arg[iarg], "tol") == 0) {
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix nve/spin/llb tol", error);
    tolerance = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    if (tolerance < 0.0) error->all(FLERR, iarg + 1, "Fix nve/spin/llb tol must be >= 0");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void FixNVESpinLLB::init()
{
  FixNVESpinSync::init();

  // fix setforce/spin overwrites the field, which has no meaning for the
  // spin magnitudes

  if (setforce_fix)
    error->all(FLERR, Error::NOLASTLINE, "Fix nve/spin/llb cannot be used with fix setforce/spin");

  modulus_fixes.clear();
  for (const auto &ifix : modify->get_fix_list())
    if (auto *mfix = dynamic_cast<SpinFieldFix *>(ifix)) modulus_fixes.push_back(mfix);

  // the longitudinal drive of the pair styles: spin pair styles whose
  // energy depends on the directions only contribute nothing

  full_pairs.clear();
  for (auto &ipair : spin_pairs)
    if (ipair->has_full_spin_gradient()) full_pairs.push_back(ipair);

  if (field_global) {
    if (spin_pairs.empty() || (full_pairs.size() == spin_pairs.size()))
      pair_long = PAIR_GLOBAL;
    else if (full_pairs.empty())
      pair_long = PAIR_NONE;
    else
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix nve/spin/llb field global cannot separate the longitudinal force of spin "
                 "pair styles with and without the full spin gradient; use field perspin");
  } else {
    pair_long = full_pairs.empty() ? PAIR_NONE : PAIR_PERSPIN;
  }
}

/* ----------------------------------------------------------------------
   damping, temperature, and random numbers of fix langevin/spin/llb, or
   the damping constants given on this fix for dynamics without a bath
------------------------------------------------------------------------- */

void FixNVESpinLLB::init_bath()
{
  bathbit = 0;
  alpha = alpha_par = temp = sigma = kt_long = 0.0;

  for (const auto &ifix : modify->get_fix_list())
    if (utils::strmatch(ifix->style, "^langevin/spin") && strcmp(ifix->style, "langevin/spin/llb"))
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix nve/spin/llb cannot be used with fix {}; use fix langevin/spin/llb",
                 ifix->style);

  const auto baths = modify->get_fix_by_style("^langevin/spin/llb$");
  if (baths.size() > 1)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix nve/spin/llb supports only one fix langevin/spin/llb");

  if (baths.empty()) {
    alpha = alpha_perp_user;
    alpha_par = alpha_par_user;
    bathbit = groupbit;
    return;
  }

  Fix *bath = baths.front();
  if (user_alpha_perp || user_alpha_par)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix nve/spin/llb takes alpha_perp and alpha_par from fix {} when it is defined",
               bath->style);
  if (bath->igroup != igroup)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix nve/spin/llb and fix langevin/spin/llb must use the same group");

  int dim = 0;
  auto *palpha_perp = (double *) bath->extract("alpha_perp", dim);
  auto *palpha_par = (double *) bath->extract("alpha_par", dim);
  auto *ptemp = (double *) bath->extract("temp", dim);
  auto *pseed = (int *) bath->extract("seed", dim);
  auto *ptnoise = (int *) bath->extract("tnoise", dim);
  auto *plnoise = (int *) bath->extract("lnoise", dim);
  if (!palpha_perp || !palpha_par || !ptemp || !pseed || !ptnoise || !plnoise)
    error->all(FLERR, Error::NOLASTLINE, "Fix nve/spin/llb cannot access the parameters of fix {}",
               bath->style);

  bathbit = groupbit;
  alpha = *palpha_perp;
  alpha_par = *palpha_par;
  temp = *ptemp;

  const double kt = force->boltz * temp;
  const bool use_tnoise = *ptnoise && (kt > 0.0) && (alpha > 0.0);
  const bool use_lnoise = *plnoise && (kt > 0.0) && (alpha_par > 0.0);
  if (use_tnoise) sigma = sqrt(2.0 * alpha * kt / (hbar * dtv));
  if (use_lnoise) kt_long = kt;
  if (use_tnoise || use_lnoise) create_random(*pseed);
}

/* ----------------------------------------------------------------------
   one timestep of the spins: half a longitudinal step, a full transverse
   step, half a longitudinal step
------------------------------------------------------------------------- */

void FixNVESpinLLB::advance_spins()
{
  if (scheme == SIB) {
    longitudinal_sib(0.5 * dtv);
    transverse_sib();
    longitudinal_sib(0.5 * dtv);
  } else {
    longitudinal_pc(0.5 * dtv);
    transverse_pc();
    longitudinal_pc(0.5 * dtv);
  }
}

/* ----------------------------------------------------------------------
   longitudinal drive L = (1/hbar)(-dE/d|S|) of local atom i at the
   current spins
------------------------------------------------------------------------- */

double FixNVESpinLLB::longitudinal_drive(int i)
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  const double smag = sp[i][3];
  if (smag <= 0.0) return 0.0;

  double drive = 0.0;

  // pair styles: the radial component of a field that is the full gradient,
  // fm = |S|/hbar (-dE/dS), is |S| L

  if (pair_long == PAIR_GLOBAL) {
    double **fm = atom->fm;
    drive += (fm[i][0] * sp[i][0] + fm[i][1] * sp[i][1] + fm[i][2] * sp[i][2]) / smag;
  } else if (pair_long == PAIR_PERSPIN) {
    for (auto &ipair : full_pairs) {
      double fmi[3] = {0.0, 0.0, 0.0};
      ipair->compute_single_pair(i, fmi);
      drive += (fmi[0] * sp[i][0] + fmi[1] * sp[i][1] + fmi[2] * sp[i][2]) / smag;
    }
  }

  // Zeeman energy of fix precession/spin, -hbar |S| (h.s); its anisotropy
  // energies depend on the spin direction only

  for (auto &ifix : precession_fixes) {
    if (!ifix->zeeman_flag || !(mask[i] & ifix->groupbit)) continue;
    double zee[3] = {0.0, 0.0, 0.0};
    ifix->compute_zeeman(i, zee);
    drive += (zee[0] * sp[i][0] + zee[1] * sp[i][1] + zee[2] * sp[i][2]) / smag;
  }

  for (auto &ifix : modulus_fixes) drive += ifix->longitudinal_drive(i);

  return drive;
}

/* ----------------------------------------------------------------------
   du/dt of local atom i, u = ln |S|, without the noise
------------------------------------------------------------------------- */

double FixNVESpinLLB::longitudinal_rate(int i)
{
  return alpha_par * (atom->sp[i][3] * longitudinal_drive(i) + 3.0 * kt_long / hbar);
}

/* ----------------------------------------------------------------------
   make the current spin magnitudes known to the ghost atoms, and re-run the
   pair style if its longitudinal drive is needed
------------------------------------------------------------------------- */

void FixNVESpinLLB::refresh_longitudinal()
{
  comm->forward_comm();
  if (pair_long == PAIR_GLOBAL) evaluate_pair();
  nevals++;
}

/* ----------------------------------------------------------------------
   start a longitudinal step of length h: the active spins (in the group,
   with a finite magnitude), their u = ln |S|, and the noise of the step
------------------------------------------------------------------------- */

void FixNVESpinLLB::start_longitudinal(double h)
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  const double sigma_long = sqrt(2.0 * alpha_par * kt_long * h / hbar);

  active.assign(nlocal, 0);
  u0.assign(nlocal, 0.0);
  u1.assign(nlocal, 0.0);
  k1.assign(nlocal, 0.0);
  lnoise.assign(nlocal, 0.0);

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || (sp[i][3] <= 0.0)) continue;
    active[i] = 1;
    u0[i] = u1[i] = log(sp[i][3]);
    if (random && (sigma_long > 0.0)) lnoise[i] = sigma_long * random->gaussian();
  }
}

/* ----------------------------------------------------------------------
   implicit midpoint step of u = ln |S| over time h, solved by iteration:
     u(t+h) = u(t) + h du/dt( exp[(u(t) + u(t+h))/2] ) + noise
------------------------------------------------------------------------- */

void FixNVESpinLLB::longitudinal_sib(double h)
{
  if (alpha_par <= 0.0) return;

  double **sp = atom->sp;
  start_longitudinal(h);
  const int n = (int) active.size();

  double residual = 0.0;
  for (int iter = 0; iter < maxiter; iter++) {
    for (int i = 0; i < n; i++)
      if (active[i]) sp[i][3] = exp(0.5 * (u0[i] + u1[i]));
    refresh_longitudinal();

    double local = 0.0;
    for (int i = 0; i < n; i++) {
      if (!active[i]) continue;
      const double unew = u0[i] + h * longitudinal_rate(i) + lnoise[i];
      local = std::max(local, fabs(unew - u1[i]));
      u1[i] = unew;
    }

    MPI_Allreduce(&local, &residual, 1, MPI_DOUBLE, MPI_MAX, world);
    if (residual <= tolerance) break;
  }

  for (int i = 0; i < n; i++)
    if (active[i]) sp[i][3] = exp(u1[i]);
  comm->forward_comm();

  nsolve++;
  if (residual > tolerance) nfail++;
  residual_max = std::max(residual_max, residual);
}

/* ----------------------------------------------------------------------
   predictor-corrector (Heun) step of u = ln |S| over time h
------------------------------------------------------------------------- */

void FixNVESpinLLB::longitudinal_pc(double h)
{
  if (alpha_par <= 0.0) return;

  double **sp = atom->sp;
  start_longitudinal(h);
  const int n = (int) active.size();

  // predictor, with every rate computed before any magnitude is changed

  refresh_longitudinal();
  for (int i = 0; i < n; i++) {
    if (!active[i]) continue;
    k1[i] = longitudinal_rate(i);
    u1[i] = u0[i] + h * k1[i] + lnoise[i];
  }
  for (int i = 0; i < n; i++)
    if (active[i]) sp[i][3] = exp(u1[i]);

  // corrector

  refresh_longitudinal();
  for (int i = 0; i < n; i++)
    if (active[i]) u1[i] = u0[i] + 0.5 * h * (k1[i] + longitudinal_rate(i)) + lnoise[i];
  for (int i = 0; i < n; i++)
    if (active[i]) sp[i][3] = exp(u1[i]);
  comm->forward_comm();
}

/* ----------------------------------------------------------------------
   field of local atom i as the force computation builds it, for the check:
   the precession field plus the radial fields of the fixes acting on the
   spin magnitudes, fm = |S| L s
------------------------------------------------------------------------- */

void FixNVESpinLLB::rebuild_field(int i, double *fmi)
{
  double **sp = atom->sp;

  spin_field(i, fmi);

  double drive = 0.0;
  for (auto &ifix : modulus_fixes) drive += ifix->longitudinal_drive(i);
  const double pref = sp[i][3] * drive;
  fmi[0] += pref * sp[i][0];
  fmi[1] += pref * sp[i][1];
  fmi[2] += pref * sp[i][2];
}
