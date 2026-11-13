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

   Common base of the spin integrators that advance all spins of the group
   together, every spin seeing the effective field of the same configuration
   (a synchronous update), as opposed to the sequential sweep of fix
   nve/spin.  A synchronous update does not depend on the order of the atoms
   or on the domain decomposition, and needs a fixed number of field
   evaluations per step instead of four per spin, which matters when the
   field comes from a pair style that has to be re-run as a whole.

   The spin directions obey the Landau-Lifshitz-Gilbert equation

       ds/dt = omega x s,   omega = [H - alpha (H x s)] / (1 + alpha^2),

   with H the precession field held in fm (rad/time) plus, at finite
   temperature, a thermal field of variance 2 alpha kB T / (hbar dt) per
   component.  The thermal field is drawn once per step and held fixed while
   the field is evaluated at several intermediate configurations, as the
   Stratonovich interpretation of the stochastic equation requires.

   The atom positions are integrated with velocity Verlet: a half-step
   velocity update and a full position update, then the spin update at the
   new positions, and the second half-step velocity update with the forces
   of the new configuration.

   The effective field at an intermediate configuration is obtained either
   spin by spin, from PairSpin::compute_single_pair() and the single-spin
   functions of fix precession/spin and fix setforce/spin, or by re-running
   the whole pair style, for pair styles that have no per-spin entry point
   (keyword field).  At the first step of each run the field is rebuilt the
   same way at the current configuration and compared with the one of the
   force computation, so that a contribution the integrator cannot
   re-evaluate stops the run instead of being dropped silently.
------------------------------------------------------------------------- */

#include "fix_nve_spin_sync.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "error.h"
#include "fix_precession_spin.h"
#include "fix_setforce_spin.h"
#include "force.h"
#include "math_const.h"
#include "modify.h"
#include "pair.h"
#include "pair_spin.h"
#include "random_mars.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_2PI;

// offset added to the seed of fix langevin/spin, so that the thermal field of
// the integrator does not repeat the random numbers of the bath itself

static constexpr int SEED_OFFSET = 571;

// largest seed accepted by the RanMars random number generator

static constexpr int MAXSEED = 900000000;

/* ---------------------------------------------------------------------- */

FixNVESpinSync::FixNVESpinSync(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), lattice_flag(1), field_mode(FIELD_AUTO), field_global(0),
    field_tol_loose(0), check_parallel(0), field_checked(0), field_check{0.0, 0.0, 0.0}, dtv(0.0),
    dtf(0.0), hbar(0.0), maxiter(10), tolerance(1.0e-12), nsteps(0), nevals(0), nsolve(0), nfail(0),
    residual_max(0.0), bathbit(0), alpha(0.0), temp(0.0), sigma(0.0), seed_offset(SEED_OFFSET),
    random_seed(0), random(nullptr), setforce_fix(nullptr)
{
  if (!atom->sp_flag) error->all(FLERR, "Fix {} requires atom style spin", style);

  time_integrate = 1;
  restart_global = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 0;
}

/* ---------------------------------------------------------------------- */

FixNVESpinSync::~FixNVESpinSync()
{
  delete random;
}

/* ----------------------------------------------------------------------
   parse the keywords common to all synchronous spin integrators, starting
   at argument iarg, and pass the others to the derived class; called by
   the derived constructors
------------------------------------------------------------------------- */

void FixNVESpinSync::parse_args(int iarg, int narg, char **arg)
{
  while (iarg < narg) {
    if (strcmp(arg[iarg], "lattice") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " lattice", error);
      const std::string value = arg[iarg + 1];
      if (value == "moving")
        lattice_flag = 1;
      else if (value == "frozen")
        lattice_flag = 0;
      else
        lattice_flag = utils::logical(FLERR, value, false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "field") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, std::string("fix ") + style + " field", error);
      const std::string value = arg[iarg + 1];
      if (value == "auto")
        field_mode = FIELD_AUTO;
      else if (value == "perspin")
        field_mode = FIELD_PERSPIN;
      else if (value == "global")
        field_mode = FIELD_GLOBAL;
      else
        error->all(FLERR, iarg + 1, "Unknown fix {} field option: {}", style, value);
      iarg += 2;
    } else {
      const int nused = parse_extra(iarg, narg, arg);
      if (nused <= 0) error->all(FLERR, iarg, "Unknown fix {} keyword: {}", style, arg[iarg]);
      iarg += nused;
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixNVESpinSync::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESpinSync::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  hbar = force->hplanck / MY_2PI;

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Fix {} is not compatible with run_style respa", style);

  // the spins can be advanced by only one integrator, and the inertial spin
  // dynamics styles integrate a different equation of motion

  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix == this) continue;
    const std::string other = ifix->style;
    if (utils::strmatch(other, "^nve/spin") || (other == "nve/tspin") || (other == "nvt/tspin") ||
        (other == "npt/tspin") || (other == "nph/tspin") || (other == "langevin/tspin"))
      error->all(FLERR, Error::NOLASTLINE, "Fix {} cannot be combined with fix {}", style, other);
  }

  // every SPIN pair style, including all sub-styles of a hybrid pair style

  spin_pairs.clear();
  for (int k = 1; force->pair; k++) {
    auto *ipair = dynamic_cast<PairSpin *>(force->pair_match("^spin", 0, k));
    if (!ipair) break;
    if (std::find(spin_pairs.begin(), spin_pairs.end(), ipair) != spin_pairs.end()) break;
    spin_pairs.push_back(ipair);
  }

  precession_fixes.clear();
  for (const auto &ifix : modify->get_fix_by_style("^precession/spin$"))
    precession_fixes.push_back(dynamic_cast<FixPrecessionSpin *>(ifix));

  setforce_fix = nullptr;
  const auto setforce = modify->get_fix_by_style("^setforce/spin$");
  if (setforce.size() > 1)
    error->all(FLERR, Error::NOLASTLINE, "Fix {} supports only one fix setforce/spin", style);
  if (!setforce.empty()) setforce_fix = dynamic_cast<FixSetForceSpin *>(setforce.front());

  init_bath();

  // how the pair field is evaluated at intermediate configurations

  if ((field_mode == FIELD_GLOBAL) && !force->pair)
    error->all(FLERR, Error::NOLASTLINE, "Fix {} field global requires a pair style", style);
  field_global = (field_mode == FIELD_GLOBAL) ? 1 : 0;

  // a pair style without per-spin field, such as a machine-learning model
  // that may run on a GPU, need not reproduce its field to the last digit

  field_tol_loose = (field_global && spin_pairs.empty()) ? 1 : 0;
  field_checked = 0;

  nsteps = nevals = nsolve = nfail = 0;
  residual_max = 0.0;
}

/* ----------------------------------------------------------------------
   transverse damping, temperature, and random numbers of fix langevin/spin
------------------------------------------------------------------------- */

void FixNVESpinSync::init_bath()
{
  bathbit = 0;
  alpha = temp = sigma = 0.0;

  for (const auto &ifix : modify->get_fix_by_style("^langevin/spin/"))
    error->all(FLERR, Error::NOLASTLINE, "Fix {} cannot be combined with fix {}", style,
               ifix->style);

  const auto baths = modify->get_fix_by_style("^langevin/spin$");
  if (baths.size() > 1)
    error->all(FLERR, Error::NOLASTLINE, "Fix {} supports only one fix langevin/spin", style);
  if (baths.empty()) return;

  Fix *bath = baths.front();
  int dim = 0;
  auto *palpha = (double *) bath->extract("alpha", dim);
  auto *ptemp = (double *) bath->extract("temp", dim);
  auto *pseed = (int *) bath->extract("seed", dim);
  if (!palpha || !ptemp || !pseed)
    error->all(FLERR, Error::NOLASTLINE, "Fix {} cannot access the parameters of fix {}", style,
               bath->style);

  bathbit = bath->groupbit;
  alpha = *palpha;
  temp = *ptemp;
  if ((alpha > 0.0) && (temp > 0.0)) {
    sigma = sqrt(2.0 * alpha * force->boltz * temp / (hbar * dtv));
    create_random(*pseed);
  }
}

/* ----------------------------------------------------------------------
   create the random number generator for a seed, unless it exists already,
   so that consecutive runs continue the same sequence; apply the state of
   a restart file the first time
------------------------------------------------------------------------- */

void FixNVESpinSync::create_random(int seed)
{
  if (random && (random_seed == seed)) return;

  if ((seed <= 0) || (seed > MAXSEED - seed_offset - comm->nprocs))
    error->all(FLERR, Error::NOLASTLINE, "Fix {} needs a random number seed between 1 and {}",
               style, MAXSEED - seed_offset - comm->nprocs);

  delete random;
  random = new RanMars(lmp, seed + seed_offset + comm->me);
  random_seed = seed;

  if (!restart_state.empty()) {
    if (((int) restart_state[0] == comm->nprocs) && ((int) restart_state[1] == seed)) {
      const int stride = RanMars::state_size(&restart_state[2]);
      random->set_state(&restart_state[2 + (size_t) comm->me * stride]);
    } else if (comm->me == 0) {
      error->warning(FLERR,
                     "Fix {} cannot restore the random number state from the restart file: "
                     "the number of processors or the seed differs",
                     style);
    }
    restart_state.clear();
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESpinSync::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  if (sigma > 0.0) sigma = sqrt(2.0 * alpha * force->boltz * temp / (hbar * dtv));
}

/* ---------------------------------------------------------------------- */

void FixNVESpinSync::initial_integrate(int /*vflag*/)
{
  // fm still holds the field of the force computation at the configuration
  // this step starts from: verify once per run that it can be rebuilt

  if (!field_checked) {
    check_field();
    field_checked = 1;
  }

  if (lattice_flag) {
    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        const double dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }

  advance_spins();
  nsteps++;
}

/* ---------------------------------------------------------------------- */

void FixNVESpinSync::final_integrate()
{
  if (!lattice_flag) return;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const double dtfm = dtf / (rmass ? rmass[i] : mass[type[i]]);
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
}

/* ----------------------------------------------------------------------
   one timestep of the spin directions with the semi-implicit midpoint
   method (Mentink et al., J. Phys.: Condens. Matter 22, 176001 (2010)):

       s(t+dt) = Cayley[ omega( s_mid, H(s_mid) ), dt ] s(t),
       s_mid   = (s(t) + s(t+dt)) / |s(t) + s(t+dt)|,

   solved by repeated iteration, one evaluation of the field of all spins
   per iteration.  The Cayley rotation keeps the spin lengths exactly.
------------------------------------------------------------------------- */

void FixNVESpinSync::transverse_sib()
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  sold.resize(3 * (size_t) nlocal);
  swork.resize(3 * (size_t) nlocal);
  draw_noise(noise, 3, sigma);

  for (int i = 0; i < nlocal; i++)
    for (int d = 0; d < 3; d++) sold[3 * i + d] = swork[3 * i + d] = sp[i][d];

  double residual = 0.0;
  for (int iter = 0; iter < maxiter; iter++) {

    // all spins of the group move to the current midpoint estimate before
    // any field is evaluated, so every spin sees the same configuration

    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      double mid[3] = {sold[3 * i] + swork[3 * i], sold[3 * i + 1] + swork[3 * i + 1],
                       sold[3 * i + 2] + swork[3 * i + 2]};
      normalize(mid);
      sp[i][0] = mid[0];
      sp[i][1] = mid[1];
      sp[i][2] = mid[2];
    }
    refresh_field();

    double local = 0.0;
    for (int i = 0; i < nlocal; i++) {
      if (!(mask[i] & groupbit)) continue;
      double fmi[3], omega[3], next[3];
      spin_field(i, fmi);
      fmi[0] += noise[3 * i];
      fmi[1] += noise[3 * i + 1];
      fmi[2] += noise[3 * i + 2];
      llg_omega(sp[i], fmi, damping(i), omega);
      cayley(&sold[3 * i], omega, dtv, next);
      normalize(next);
      for (int d = 0; d < 3; d++) {
        local = std::max(local, fabs(next[d] - swork[3 * i + d]));
        swork[3 * i + d] = next[d];
      }
    }

    MPI_Allreduce(&local, &residual, 1, MPI_DOUBLE, MPI_MAX, world);
    if (residual <= tolerance) break;
  }

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    sp[i][0] = swork[3 * i];
    sp[i][1] = swork[3 * i + 1];
    sp[i][2] = swork[3 * i + 2];
  }
  comm->forward_comm();

  nsolve++;
  if (residual > tolerance) nfail++;
  residual_max = std::max(residual_max, residual);
}

/* ----------------------------------------------------------------------
   one timestep of the spin directions with the predictor-corrector (Heun)
   method, two evaluations of the field of all spins:

       k1 = ds/dt( s(t), H(s(t)) ),   s_p = normalized( s(t) + dt k1 )
       k2 = ds/dt( s_p, H(s_p) ),     s(t+dt) = normalized( s(t) + dt (k1+k2)/2 )

   Every slope is computed before any spin is moved, so that the update is
   independent of the order of the atoms.
------------------------------------------------------------------------- */

void FixNVESpinSync::transverse_pc()
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  sold.resize(3 * (size_t) nlocal);
  swork.resize(3 * (size_t) nlocal);
  snew.resize(3 * (size_t) nlocal);
  draw_noise(noise, 3, sigma);

  for (int i = 0; i < nlocal; i++)
    for (int d = 0; d < 3; d++) sold[3 * i + d] = sp[i][d];

  // predictor: slope k1 in swork, predicted spin in snew

  refresh_field();
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    double fmi[3];
    spin_field(i, fmi);
    fmi[0] += noise[3 * i];
    fmi[1] += noise[3 * i + 1];
    fmi[2] += noise[3 * i + 2];
    llg_rate(&sold[3 * i], fmi, damping(i), &swork[3 * i]);
    for (int d = 0; d < 3; d++) snew[3 * i + d] = sold[3 * i + d] + dtv * swork[3 * i + d];
    normalize(&snew[3 * i]);
  }

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    sp[i][0] = snew[3 * i];
    sp[i][1] = snew[3 * i + 1];
    sp[i][2] = snew[3 * i + 2];
  }

  // corrector: average of the two slopes, applied to s(t)

  refresh_field();
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    double fmi[3], k2[3];
    spin_field(i, fmi);
    fmi[0] += noise[3 * i];
    fmi[1] += noise[3 * i + 1];
    fmi[2] += noise[3 * i + 2];
    llg_rate(&snew[3 * i], fmi, damping(i), k2);
    for (int d = 0; d < 3; d++)
      swork[3 * i + d] = sold[3 * i + d] + 0.5 * dtv * (swork[3 * i + d] + k2[d]);
    normalize(&swork[3 * i]);
  }

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    sp[i][0] = swork[3 * i];
    sp[i][1] = swork[3 * i + 1];
    sp[i][2] = swork[3 * i + 2];
  }
  comm->forward_comm();
}

/* ----------------------------------------------------------------------
   make the current spins (and positions) known to the ghost atoms, and
   re-run the pair style if its field is evaluated as a whole
------------------------------------------------------------------------- */

void FixNVESpinSync::refresh_field()
{
  comm->forward_comm();
  if (field_global) evaluate_pair();
  nevals++;
}

/* ----------------------------------------------------------------------
   re-run the pair style at the current configuration, which must have been
   communicated to the ghost atoms; afterwards fm holds the field of the pair
   style alone.  The atom forces are restored, since other fixes may still
   use them during this step.
------------------------------------------------------------------------- */

void FixNVESpinSync::evaluate_pair()
{
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  double **f = atom->f;

  fsave.resize(3 * (size_t) nlocal);
  for (int i = 0; i < nlocal; i++) {
    fsave[3 * i] = f[i][0];
    fsave[3 * i + 1] = f[i][1];
    fsave[3 * i + 2] = f[i][2];
  }

  const size_t nbytes = sizeof(double) * nall;
  if (nbytes) {
    memset(&f[0][0], 0, 3 * nbytes);
    if (atom->avec->forceclearflag) atom->avec->force_clear(0, nbytes);
  }

  force->pair->compute(0, 0);
  if (force->newton_pair) comm->reverse_comm();

  for (int i = 0; i < nlocal; i++) {
    f[i][0] = fsave[3 * i];
    f[i][1] = fsave[3 * i + 1];
    f[i][2] = fsave[3 * i + 2];
  }
}

/* ----------------------------------------------------------------------
   precession field on local atom i at the current spins: the field of the
   pair styles, fix precession/spin, and fix setforce/spin.  Fixes acting on
   the spin magnitudes only add radial fields, which do not move the spin
   directions and are not included.
------------------------------------------------------------------------- */

void FixNVESpinSync::spin_field(int i, double *fmi)
{
  double **sp = atom->sp;

  fmi[0] = fmi[1] = fmi[2] = 0.0;

  if (field_global) {
    double **fm = atom->fm;
    fmi[0] = fm[i][0];
    fmi[1] = fm[i][1];
    fmi[2] = fm[i][2];
  } else {
    for (auto &ipair : spin_pairs) ipair->compute_single_pair(i, fmi);
  }

  double spi[3] = {sp[i][0], sp[i][1], sp[i][2]};
  for (auto &ifix : precession_fixes) ifix->compute_single_precession(i, spi, fmi);

  if (setforce_fix) setforce_fix->single_setforce_spin(i, fmi);
}

/* ----------------------------------------------------------------------
   field of local atom i as the force computation builds it, for the check
------------------------------------------------------------------------- */

void FixNVESpinSync::rebuild_field(int i, double *fmi)
{
  spin_field(i, fmi);
}

/* ----------------------------------------------------------------------
   rebuild the field at the configuration of the force computation and
   compare it with atom->fm.  Unless check_parallel is set, only the
   components perpendicular to the spins are compared, the only ones that
   act on the spin directions.
------------------------------------------------------------------------- */

void FixNVESpinSync::check_field()
{
  double **sp = atom->sp;
  double **fm = atom->fm;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  std::vector<double> fmpost(3 * (size_t) nlocal);
  for (int i = 0; i < nlocal; i++)
    for (int d = 0; d < 3; d++) fmpost[3 * i + d] = fm[i][d];

  // a pair style evaluated on a GPU may not reproduce its own field exactly:
  // run it twice and take the difference as the floor of the comparison

  double local[3] = {0.0, 0.0, 0.0};    // deviation, field scale, re-run spread
  if (field_global) {
    evaluate_pair();
    std::vector<double> fmfirst(3 * (size_t) nlocal);
    for (int i = 0; i < nlocal; i++)
      for (int d = 0; d < 3; d++) fmfirst[3 * i + d] = fm[i][d];
    evaluate_pair();
    for (int i = 0; i < nlocal; i++)
      for (int d = 0; d < 3; d++)
        local[2] = std::max(local[2], fabs(fm[i][d] - fmfirst[3 * i + d]));
  }

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    double fmi[3], diff[3], post[3];
    rebuild_field(i, fmi);
    for (int d = 0; d < 3; d++) {
      post[d] = fmpost[3 * i + d];
      diff[d] = fmi[d] - post[d];
    }
    if (!check_parallel) {
      const double ddot = diff[0] * sp[i][0] + diff[1] * sp[i][1] + diff[2] * sp[i][2];
      const double pdot = post[0] * sp[i][0] + post[1] * sp[i][1] + post[2] * sp[i][2];
      for (int d = 0; d < 3; d++) {
        diff[d] -= ddot * sp[i][d];
        post[d] -= pdot * sp[i][d];
      }
    }
    for (int d = 0; d < 3; d++) {
      local[0] = std::max(local[0], fabs(diff[d]));
      local[1] = std::max(local[1], fabs(post[d]));
    }
  }

  for (int i = 0; i < nlocal; i++)
    for (int d = 0; d < 3; d++) fm[i][d] = fmpost[3 * i + d];

  MPI_Allreduce(local, field_check, 3, MPI_DOUBLE, MPI_MAX, world);

  double tol = std::max(1.0e-8 * field_check[1], 10.0 * field_check[2]);
  if (field_tol_loose) tol = std::max(tol, 1.0e-4 * field_check[1]);

  if (field_check[0] > tol)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix {} cannot rebuild the magnetic field of the force computation: largest "
               "deviation {:.6g} for a field scale of {:.6g} (tolerance {:.3g}).  A pair style, "
               "fix, or long-range solver contributes to the magnetic force in a way this "
               "integrator cannot re-evaluate{}",
               style, field_check[0], field_check[1], tol,
               field_global ? "" : "; for a pair style without a per-spin field use field global");
}

/* ----------------------------------------------------------------------
   ncomp random numbers of standard deviation sig for each local atom of
   the group in the bath
------------------------------------------------------------------------- */

void FixNVESpinSync::draw_noise(std::vector<double> &out, int ncomp, double sig)
{
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  out.assign((size_t) ncomp * nlocal, 0.0);
  if (!random || (sig <= 0.0)) return;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit) || !(mask[i] & bathbit)) continue;
    for (int k = 0; k < ncomp; k++) out[(size_t) ncomp * i + k] = sig * random->gaussian();
  }
}

/* ---------------------------------------------------------------------- */

double FixNVESpinSync::damping(int i) const
{
  return (atom->mask[i] & bathbit) ? alpha : 0.0;
}

/* ----------------------------------------------------------------------
   Landau-Lifshitz-Gilbert precession vector of spin s in field h
     omega = [h - a (h x s)] / (1 + a^2)
------------------------------------------------------------------------- */

void FixNVESpinSync::llg_omega(const double *s, const double *h, double a, double *omega)
{
  const double hxs[3] = {h[1] * s[2] - h[2] * s[1], h[2] * s[0] - h[0] * s[2],
                         h[0] * s[1] - h[1] * s[0]};
  const double pref = 1.0 / (1.0 + a * a);
  omega[0] = (h[0] - a * hxs[0]) * pref;
  omega[1] = (h[1] - a * hxs[1]) * pref;
  omega[2] = (h[2] - a * hxs[2]) * pref;
}

/* ----------------------------------------------------------------------
   time derivative of spin s in field h:  ds/dt = omega x s
------------------------------------------------------------------------- */

void FixNVESpinSync::llg_rate(const double *s, const double *h, double a, double *rate)
{
  double omega[3];
  llg_omega(s, h, a, omega);
  rate[0] = omega[1] * s[2] - omega[2] * s[1];
  rate[1] = omega[2] * s[0] - omega[0] * s[2];
  rate[2] = omega[0] * s[1] - omega[1] * s[0];
}

/* ----------------------------------------------------------------------
   rotate s about omega by the Cayley transform over time dt: the implicit
   midpoint step of ds/dt = omega x s for a fixed omega, which keeps |s|
------------------------------------------------------------------------- */

void FixNVESpinSync::cayley(const double *s, const double *omega, double dt, double *out)
{
  const double a[3] = {0.5 * dt * omega[0], 0.5 * dt * omega[1], 0.5 * dt * omega[2]};
  const double a2 = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  const double axs[3] = {a[1] * s[2] - a[2] * s[1], a[2] * s[0] - a[0] * s[2],
                         a[0] * s[1] - a[1] * s[0]};
  const double ads = a[0] * s[0] + a[1] * s[1] + a[2] * s[2];
  const double inv = 1.0 / (1.0 + a2);
  out[0] = ((1.0 - a2) * s[0] + 2.0 * axs[0] + 2.0 * ads * a[0]) * inv;
  out[1] = ((1.0 - a2) * s[1] + 2.0 * axs[1] + 2.0 * ads * a[1]) * inv;
  out[2] = ((1.0 - a2) * s[2] + 2.0 * axs[2] + 2.0 * ads * a[2]) * inv;
}

/* ----------------------------------------------------------------------
   scale v to unit length; a zero vector, as for a non-magnetic atom, is
   left unchanged
------------------------------------------------------------------------- */

void FixNVESpinSync::normalize(double *v)
{
  const double n2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  if (n2 <= 1.0e-30) return;
  const double inv = 1.0 / sqrt(n2);
  v[0] *= inv;
  v[1] *= inv;
  v[2] *= inv;
}

/* ----------------------------------------------------------------------
   write the states of the random number generators of all processors
   layout: [0] = number of stored states (0 without random numbers),
           [1] = seed, then one state per processor
------------------------------------------------------------------------- */

void FixNVESpinSync::write_restart(FILE *fp)
{
  const int nstates = random ? comm->nprocs : 0;
  const int nsize = 2 + RanMars::STATE_SIZE * nstates;

  std::vector<double> list;
  if (comm->me == 0) {
    list.resize(nsize);
    list[0] = nstates;
    list[1] = random_seed;
  }

  if (random) {
    double state[RanMars::STATE_SIZE];
    random->get_state(state);
    MPI_Gather(state, RanMars::STATE_SIZE, MPI_DOUBLE, (comm->me == 0) ? &list[2] : nullptr,
               RanMars::STATE_SIZE, MPI_DOUBLE, 0, world);
  }

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list.data(), sizeof(double), nsize, fp);
  }
}

/* ----------------------------------------------------------------------
   keep the random number states of a restart file until the generator is
   created in init()
------------------------------------------------------------------------- */

void FixNVESpinSync::restart(char *buf)
{
  auto *list = (double *) buf;
  const int nstates = (int) list[0];

  restart_state.clear();
  if (nstates > 0) {
    const int stride = RanMars::state_size(list + 2);
    restart_state.assign(list, list + 2 + (size_t) nstates * stride);
  }
}

/* ----------------------------------------------------------------------
   1-3 = field check at the start of the run: largest deviation, field
         scale, spread of two re-runs of the pair style
   4 = average number of field evaluations per timestep in this run
   5 = fraction of the midpoint solves of this run without convergence
   6 = largest final change of a midpoint solve in this run
------------------------------------------------------------------------- */

double FixNVESpinSync::compute_vector(int n)
{
  if (n < 3) return field_check[n];
  if (n == 3) return nsteps ? (double) nevals / (double) nsteps : 0.0;
  if (n == 4) return nsolve ? (double) nfail / (double) nsolve : 0.0;
  return residual_max;
}

/* ---------------------------------------------------------------------- */

void *FixNVESpinSync::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "lattice_flag") == 0) return &lattice_flag;
  return nullptr;
}
