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

   Metropolis Monte Carlo sampling of the atomic spins at temperature T.

   Each sweep attempts one move per spin of the group, in a random order.
   A move rotates the spin direction to a point drawn uniformly inside a
   cone around the current direction; the cone has the same solid angle
   seen from either end, so the proposal is symmetric.  With the modulus
   keyword the spin magnitude also takes a Gaussian random-walk step, which
   is symmetric as well; since the spin then lives in three dimensions, the
   weight |S|^2 of d^3S = |S|^2 d|S| dOmega enters the acceptance
   probability.  The cone half-angle can be adapted towards a target
   acceptance rate (Alzate-Cardona et al., J. Phys.: Condens. Matter 31,
   095802 (2019)).

   The fix moves only the spins.  The atoms can be integrated by any other
   fix, such as fix nve or fix nvt, for coupled spin-lattice sampling: the
   sweeps are done in pre_force(), after the atoms have moved and before
   the forces are computed, so the forces always correspond to the current
   spins.

   Energies:
     energy local    the energy change of a move is evaluated from the bonds
                     of the moved spin only, with the energy functions and
                     the cutoffs of each pair of atom types of the spin pair
                     styles spin/exchange, spin/exchange/biquadratic and
                     spin/neel, plus the single-site terms of fix
                     precession/spin and of fixes implementing SpinFieldFix
                     (fix landau/spin, fix spring/tspin).  In serial, the
                     periodic images of a moved spin are updated with it.  In
                     parallel, the subdomains are divided into sectors as in
                     fix nve/spin, so that spins attempted at the same time
                     on different processors do not interact.
     energy global   every move is scored by re-running the whole pair
                     style, for pair styles without a local energy, such as
                     machine-learning potentials.  The moves are attempted
                     one after another for the whole system, each pair
                     evaluation is shared by all processors, and all
                     processors draw the same random numbers, so the
                     trajectory does not depend on the number of processors.
------------------------------------------------------------------------- */

#include "fix_mc_spin.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_precession_spin.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "pair_spin.h"
#include "pair_spin_exchange.h"
#include "pair_spin_exchange_biquadratic.h"
#include "pair_spin_neel.h"
#include "random_mars.h"
#include "spin_field_fix.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <unordered_map>

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

FixMCSpin::FixMCSpin(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), temperature(0.0), beta(0.0), seed(0), nsweep(1), every(1),
    cone(30.0 * MY_PI / 180.0), cone_target(0.5), adapt_flag(1), modulus_flag(0), jacobian_flag(1),
    mod_sigma(0.05), mod_min(0.1), mod_max(5.0), energy_global(0), lattice_flag(-1), naccept(0.0),
    nattempt(0.0), mod_naccept(0.0), mod_nattempt(0.0), accept_rate(0.0), mod_accept_rate(0.0),
    epair_cur(0.0), list(nullptr), random(nullptr), image_flag(0), sector_flag(0), nsectors(1),
    nlocal_max(0), rsec{0.0, 0.0, 0.0}, stack_head(nullptr), stack_foot(nullptr),
    forward_stacks(nullptr), backward_stacks(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix mc/spin", error);
  if (!atom->sp_flag) error->all(FLERR, "Fix mc/spin requires atom style spin");
  if (!atom->tag_enable) error->all(FLERR, "Fix mc/spin requires atom IDs");

  temperature = utils::numeric(FLERR, arg[3], false, lmp);
  if (temperature <= 0.0) error->all(FLERR, 3, "Fix mc/spin temperature must be > 0");
  seed = utils::inumeric(FLERR, arg[4], false, lmp);
  if (seed <= 0) error->all(FLERR, 4, "Fix mc/spin seed must be > 0");

  int iarg = 5;
  while (iarg < narg) {
    const std::string key = arg[iarg];
    if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix mc/spin " + key, error);
    const std::string value = arg[iarg + 1];
    if (key == "cone") {
      cone = utils::numeric(FLERR, value, false, lmp) * MY_PI / 180.0;
      if ((cone <= 0.0) || (cone > MY_PI))
        error->all(FLERR, iarg + 1, "Fix mc/spin cone must be > 0 and <= 180 degrees");
    } else if (key == "adapt") {
      adapt_flag = utils::logical(FLERR, value, false, lmp);
    } else if (key == "target") {
      cone_target = utils::numeric(FLERR, value, false, lmp);
      if ((cone_target <= 0.0) || (cone_target >= 1.0))
        error->all(FLERR, iarg + 1, "Fix mc/spin target must be > 0 and < 1");
    } else if (key == "nsweep") {
      nsweep = utils::inumeric(FLERR, value, false, lmp);
      if (nsweep < 1) error->all(FLERR, iarg + 1, "Fix mc/spin nsweep must be > 0");
    } else if (key == "every") {
      every = utils::inumeric(FLERR, value, false, lmp);
      if (every < 1) error->all(FLERR, iarg + 1, "Fix mc/spin every must be > 0");
    } else if (key == "modulus") {
      if (value == "off")
        modulus_flag = 0;
      else if (value == "walk")
        modulus_flag = 1;
      else
        error->all(FLERR, iarg + 1, "Unknown fix mc/spin modulus option: {}", value);
    } else if (key == "modsigma") {
      mod_sigma = utils::numeric(FLERR, value, false, lmp);
      if (mod_sigma <= 0.0) error->all(FLERR, iarg + 1, "Fix mc/spin modsigma must be > 0");
    } else if (key == "modmin") {
      mod_min = utils::numeric(FLERR, value, false, lmp);
    } else if (key == "modmax") {
      mod_max = utils::numeric(FLERR, value, false, lmp);
    } else if (key == "jacobian") {
      jacobian_flag = utils::logical(FLERR, value, false, lmp);
    } else if (key == "energy") {
      if (value == "local")
        energy_global = 0;
      else if (value == "global")
        energy_global = 1;
      else
        error->all(FLERR, iarg + 1, "Unknown fix mc/spin energy option: {}", value);
    } else {
      error->all(FLERR, iarg, "Unknown fix mc/spin keyword: {}", key);
    }
    iarg += 2;
  }

  if ((mod_min < 0.0) || (mod_max <= mod_min))
    error->all(FLERR, "Fix mc/spin requires 0 <= modmin < modmax");

  vector_flag = 1;
  size_vector = 4;
  global_freq = 1;
  extvector = 0;
  restart_global = 1;

  // with energy global, all processors decide every move together and draw
  // the same random numbers; with energy local, each draws its own

  random = new RanMars(lmp, energy_global ? seed : seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixMCSpin::~FixMCSpin()
{
  delete random;
  memory->destroy(stack_head);
  memory->destroy(stack_foot);
  memory->destroy(forward_stacks);
  memory->destroy(backward_stacks);
}

/* ---------------------------------------------------------------------- */

int FixMCSpin::setmask()
{
  int mask = 0;
  mask |= PRE_NEIGHBOR;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::init()
{
  beta = 1.0 / (force->boltz * temperature);

  if (utils::strmatch(update->integrate_style, "^respa"))
    error->all(FLERR, Error::NOLASTLINE, "Fix mc/spin is not compatible with run_style respa");

  // a magnitude outside the allowed interval is never moved back into it

  if (modulus_flag) {
    double **sp = atom->sp;
    int *mask = atom->mask;
    int outside = 0, outside_all = 0;
    for (int i = 0; i < atom->nlocal; i++)
      if ((mask[i] & groupbit) && (sp[i][3] > 0.0) &&
          ((sp[i][3] < mod_min) || (sp[i][3] > mod_max)))
        outside = 1;
    MPI_Allreduce(&outside, &outside_all, 1, MPI_INT, MPI_MAX, world);
    if (outside_all && (comm->me == 0))
      error->warning(FLERR,
                     "Fix mc/spin: some spin magnitudes lie outside modmin..modmax and "
                     "will not change");
  }

  // the spins may be moved by one method only, and fields that are not
  // energies have no place in a Metropolis test

  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix == this) continue;
    const std::string other = ifix->style;
    if (utils::strmatch(other, "^nve/spin") || utils::strmatch(other, "^langevin/spin") ||
        utils::strmatch(other, "^mc/spin") || (other == "setforce/spin") ||
        (other == "nve/tspin") || (other == "nvt/tspin") || (other == "npt/tspin") ||
        (other == "nph/tspin") || (other == "langevin/tspin"))
      error->all(FLERR, Error::NOLASTLINE, "Fix mc/spin cannot be combined with fix {}", other);
  }

  precession_fixes.clear();
  for (const auto &ifix : modify->get_fix_by_style("^precession/spin$")) {
    auto *pfix = dynamic_cast<FixPrecessionSpin *>(ifix);
    if (pfix->stt_flag)
      error->all(FLERR, Error::NOLASTLINE,
                 "Fix mc/spin cannot sample the spin transfer torque of fix precession/spin");
    precession_fixes.push_back(pfix);
  }

  modulus_fixes.clear();
  for (const auto &ifix : modify->get_fix_list())
    if (auto *mfix = dynamic_cast<SpinFieldFix *>(ifix)) modulus_fixes.push_back(mfix);

  // spin pair styles, their local energy functions, and their cutoffs for
  // each pair of atom types; a sub-style of pair_style hybrid applies only
  // to the pairs of types assigned to it, the others get a negative cutoff.
  // The coefficients of a sub-style stay set when a later pair_coeff command
  // assigns its pair of types to another sub-style or to none, so the
  // assignment is taken from pair_style hybrid itself.

  spin_pairs.clear();
  pair_kind.clear();
  pair_cutsq.clear();
  const int ntypes = atom->ntypes;
  auto *hybrid = dynamic_cast<PairHybrid *>(force->pair);
  double cutmax = 0.0;
  for (int k = 1; force->pair; k++) {
    auto *ipair = dynamic_cast<PairSpin *>(force->pair_match("^spin", 0, k));
    if (!ipair) break;
    if (std::find(spin_pairs.begin(), spin_pairs.end(), ipair) != spin_pairs.end()) break;
    int kind = -1;
    const char *cutname = "cut_spin_exchange";
    std::string substyle = "spin/exchange";
    if (dynamic_cast<PairSpinExchangeBiquadratic *>(ipair)) {
      kind = BIQUADRATIC;
      substyle = "spin/exchange/biquadratic";
    } else if (dynamic_cast<PairSpinExchange *>(ipair)) {
      kind = EXCHANGE;
    } else if (dynamic_cast<PairSpinNeel *>(ipair)) {
      kind = NEEL;
      cutname = "cut_spin_neel";
      substyle = "spin/neel";
    }

    std::vector<double> cutsq((size_t) (ntypes + 1) * (ntypes + 1), -1.0);
    int dim = 0;
    auto **cut = (kind < 0) ? nullptr : (double **) ipair->extract(cutname, dim);
    if (cut && (dim == 2)) {
      for (int itype = 1; itype <= ntypes; itype++) {
        for (int jtype = 1; jtype <= ntypes; jtype++) {
          if (!ipair->setflag[MIN(itype, jtype)][MAX(itype, jtype)]) continue;
          if (hybrid && !hybrid->check_ijtype(itype, jtype, substyle.data())) continue;
          cutsq[itype * (ntypes + 1) + jtype] = cut[itype][jtype] * cut[itype][jtype];
          cutmax = MAX(cutmax, cut[itype][jtype]);
        }
      }
    }
    spin_pairs.push_back(ipair);
    pair_kind.push_back(kind);
    pair_cutsq.push_back(cutsq);
  }

  image_flag = sector_flag = 0;
  nsectors = 1;
  if (!energy_global) {
    for (std::size_t k = 0; k < spin_pairs.size(); k++)
      if (pair_kind[k] < 0)
        error->all(FLERR, Error::NOLASTLINE,
                   "Fix mc/spin energy local supports the pair styles spin/exchange, "
                   "spin/exchange/biquadratic, and spin/neel; use energy global for others");
    if (force->pair && spin_pairs.empty() && strcmp(force->pair_style, "zero") && (comm->me == 0))
      error->warning(FLERR,
                     "Fix mc/spin energy local found no spin pair style; the energy of "
                     "pair style {} is not sampled",
                     force->pair_style);

    if (!spin_pairs.empty()) neighbor->add_request(this, NeighConst::REQ_FULL);

    // in serial, the periodic images of a moved spin are updated with it; in
    // parallel, interacting spins are swept sector by sector, and spins
    // without interactions all at once

    if (comm->nprocs == 1) {
      image_flag = 1;
    } else if (!spin_pairs.empty()) {
      sector_flag = 1;
      rsec[0] = cutmax + neighbor->skin;
      sectoring();
    }
    memory->grow(stack_head, nsectors, "mc/spin:stack_head");
    memory->grow(stack_foot, nsectors, "mc/spin:stack_foot");
  } else {
    if (!force->pair)
      error->all(FLERR, Error::NOLASTLINE, "Fix mc/spin energy global requires a pair style");
    if (atom->natoms > MAXSMALLINT)
      error->all(FLERR, Error::NOLASTLINE, "Fix mc/spin energy global supports at most {} atoms",
                 MAXSMALLINT);
  }
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ----------------------------------------------------------------------
   verify the energy the Metropolis test uses against the energy LAMMPS
   tallies, at the configuration of the setup force computation
------------------------------------------------------------------------- */

void FixMCSpin::setup(int /*vflag*/)
{
  if (energy_global)
    verify_energy_global();
  else
    verify_energy_local();
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::verify_energy_local()
{
  // the pair energies are tallied only if the setup computed energies; the
  // single-site terms use the energy functions of the fixes themselves

  if (update->eflag_global != update->ntimestep) return;

  double **sp = atom->sp;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double local[2] = {0.0, 0.0};
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) local[0] += bond_energy(i, sp[i], 0);
  for (auto &ipair : spin_pairs) {
    auto *pspin = dynamic_cast<PairSpin *>(ipair);
    if (pspin && pspin->emag)
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) local[1] += pspin->emag[i];
  }

  double all[2];
  MPI_Allreduce(local, all, 2, MPI_DOUBLE, MPI_SUM, world);

  if (fabs(all[0] - all[1]) > 1.0e-8 * MAX(fabs(all[1]), 1.0))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix mc/spin: the pair energy from local sums {:.12g} differs from the energy of "
               "the spin pair styles {:.12g}; use energy global",
               all[0], all[1]);
}

/* ----------------------------------------------------------------------
   the energy change of a move is the difference of two totals with energy
   global, so the pair style must reproduce its energy: evaluate the setup
   configuration twice and compare with the energy of the setup
------------------------------------------------------------------------- */

void FixMCSpin::verify_energy_global()
{
  const bool have_ref = (update->eflag_global == update->ntimestep);
  const double eng_vdwl = force->pair->eng_vdwl;
  const double eng_coul = force->pair->eng_coul;
  const double eref_local = eng_vdwl + eng_coul;
  double eref = 0.0;
  MPI_Allreduce(&eref_local, &eref, 1, MPI_DOUBLE, MPI_SUM, world);

  // the forces of the setup are still needed by the first timestep

  const int nall = atom->nlocal + atom->nghost;
  std::vector<double> fsave(&atom->f[0][0], &atom->f[0][0] + 3 * (size_t) nall);
  std::vector<double> fmsave(&atom->fm[0][0], &atom->fm[0][0] + 3 * (size_t) nall);
  std::vector<double> fmlsave(&atom->fm_long[0][0], &atom->fm_long[0][0] + 3 * (size_t) nall);

  const double e1 = model_energy();
  const double e2 = model_energy();

  std::copy(fsave.begin(), fsave.end(), &atom->f[0][0]);
  std::copy(fmsave.begin(), fmsave.end(), &atom->fm[0][0]);
  std::copy(fmlsave.begin(), fmlsave.end(), &atom->fm_long[0][0]);
  force->pair->eng_vdwl = eng_vdwl;
  force->pair->eng_coul = eng_coul;

  const double kt = 1.0 / beta;
  const double spread = fabs(e1 - e2);
  if (!std::isfinite(e1) || !std::isfinite(e2))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix mc/spin: the pair style returned an energy that is not a number");
  if (spread > 0.01 * kt)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix mc/spin energy global: the same configuration gives the energies {:.12g} and "
               "{:.12g}, a difference of {:.3g} kT",
               e1, e2, spread / kt);
  if (have_ref && (fabs(e1 - eref) > MAX(10.0 * spread, 0.01 * kt)))
    error->all(FLERR, Error::NOLASTLINE,
               "Fix mc/spin energy global: the pair energy {:.12g} differs from the energy {:.12g} "
               "of the force computation by {:.3g} kT",
               e1, eref, fabs(e1 - eref) / kt);
  epair_cur = e1;
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::setup_pre_neighbor()
{
  pre_neighbor();
}

/* ----------------------------------------------------------------------
   periodic images of the local atoms (serial), or the order of the local
   atoms in each sector (parallel), as in fix nve/spin
------------------------------------------------------------------------- */

void FixMCSpin::pre_neighbor()
{
  if (image_flag) link_images();
  if (!sector_flag) return;

  double **x = atom->x;
  const int nlocal = atom->nlocal;

  if (nlocal_max < nlocal) {
    nlocal_max = nlocal;
    memory->grow(backward_stacks, nlocal_max, "mc/spin:backward_stacks");
    memory->grow(forward_stacks, nlocal_max, "mc/spin:forward_stacks");
  }

  for (int j = 0; j < nsectors; j++) stack_head[j] = stack_foot[j] = -1;

  for (int j = 0; j < nsectors; j++)
    for (int i = 0; i < nlocal; i++) {
      if (coords2sector(x[i]) != j) continue;
      backward_stacks[i] = stack_head[j];
      stack_head[j] = i;
    }

  for (int j = nsectors - 1; j >= 0; j--)
    for (int i = nlocal - 1; i >= 0; i--) {
      if (coords2sector(x[i]) != j) continue;
      forward_stacks[i] = stack_foot[j];
      stack_foot[j] = i;
    }
}

/* ----------------------------------------------------------------------
   link every ghost atom to the local atom it is a periodic image of; in
   serial, all ghost atoms are such images
------------------------------------------------------------------------- */

void FixMCSpin::link_images()
{
  tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;

  std::unordered_map<tagint, int> owner;
  owner.reserve(nlocal);
  for (int i = 0; i < nlocal; i++) owner[tag[i]] = i;

  image_head.assign(nlocal, -1);
  image_next.assign(nall, -1);
  for (int j = nlocal; j < nall; j++) {
    const auto it = owner.find(tag[j]);
    if (it == owner.end()) continue;
    image_next[j] = image_head[it->second];
    image_head[it->second] = j;
  }
}

/* ----------------------------------------------------------------------
   Monte Carlo sweeps at the current atom positions, before the forces
------------------------------------------------------------------------- */

void FixMCSpin::pre_force(int /*vflag*/)
{
  if (update->ntimestep % every) return;

  // the image links are made at every neighbor list build; renew them if
  // the atoms changed without one

  if (image_flag &&
      ((image_head.size() != (size_t) atom->nlocal) ||
       (image_next.size() != (size_t) (atom->nlocal + atom->nghost))))
    link_images();

  naccept = nattempt = mod_naccept = mod_nattempt = 0.0;
  for (int k = 0; k < nsweep; k++) {
    if (energy_global)
      global_sweep();
    else
      local_sweep();
  }

  // re-running the pair style left forces behind; the force computation
  // that follows expects them cleared

  if (energy_global) clear_forces();
  mc_diagnostics();
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::local_sweep()
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  // sites without a magnetic moment are skipped, as with energy global

  if (sector_flag) {
    for (int j = 0; j < nsectors; j++) {
      comm->forward_comm();
      order.clear();
      for (int i = stack_foot[j]; i >= 0; i = forward_stacks[i])
        if ((mask[i] & groupbit) && (sp[i][3] != 0.0)) order.push_back(i);
      shuffle_and_try((int) order.size());
    }
    comm->forward_comm();
  } else {
    order.clear();
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (sp[i][3] != 0.0)) order.push_back(i);
    shuffle_and_try((int) order.size());

    // in parallel without spin pair styles, the spins do not interact and
    // their ghost copies are updated once at the end

    if (!image_flag) comm->forward_comm();
  }
}

/* ----------------------------------------------------------------------
   one sweep with energy global: the spins with a magnetic moment are tried
   one at a time, in a random order of their atom IDs that all processors
   draw alike
------------------------------------------------------------------------- */

void FixMCSpin::global_sweep()
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  const int nlocal = atom->nlocal;
  const int nprocs = comm->nprocs;

  std::vector<tagint> mytags;
  std::unordered_map<tagint, int> owned;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    if (sp[i][3] == 0.0) continue;    // a site without a magnetic moment
    mytags.push_back(tag[i]);
    owned[tag[i]] = i;
  }

  int nmine = (int) mytags.size();
  std::vector<int> counts(nprocs), displs(nprocs);
  MPI_Allgather(&nmine, 1, MPI_INT, counts.data(), 1, MPI_INT, world);
  int ntrial = 0;
  for (int p = 0; p < nprocs; p++) {
    displs[p] = ntrial;
    ntrial += counts[p];
  }
  trial_tags.resize(ntrial);
  MPI_Allgatherv(mytags.data(), nmine, MPI_LMP_TAGINT, trial_tags.data(), counts.data(),
                 displs.data(), MPI_LMP_TAGINT, world);
  std::sort(trial_tags.begin(), trial_tags.end());

  epair_cur = model_energy();

  for (int k = ntrial - 1; k > 0; k--) {
    int r = static_cast<int>(random->uniform() * (k + 1));
    if (r > k) r = k;
    std::swap(trial_tags[k], trial_tags[r]);
  }
  for (int k = 0; k < ntrial; k++) {
    const auto it = owned.find(trial_tags[k]);
    mc_global_spin((it == owned.end()) ? -1 : it->second);
  }

  // a rejected move restored the spin of its owner only

  comm->forward_comm();
}

/* ----------------------------------------------------------------------
   attempt one move for each of the n atoms in order, in a random order
------------------------------------------------------------------------- */

void FixMCSpin::shuffle_and_try(int n)
{
  for (int k = n - 1; k > 0; k--) {
    int r = static_cast<int>(random->uniform() * (k + 1));
    if (r > k) r = k;
    std::swap(order[k], order[r]);
  }
  for (int k = 0; k < n; k++) mc_single_spin(order[k]);
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::mc_single_spin(int i)
{
  double **sp = atom->sp;
  const double sold[3] = {sp[i][0], sp[i][1], sp[i][2]};
  const double mold = sp[i][3];

  double snew[3];
  propose_direction(sold, snew);

  double mnew = mold;
  int mod_move = 0;
  if (modulus_flag) {
    // reject rather than resample outside the bounds, which would make the
    // walk asymmetric

    mnew = mold + mod_sigma * random->gaussian();
    if ((mnew < mod_min) || (mnew > mod_max))
      mnew = mold;
    else
      mod_move = 1;
  }

  const double eold = single_spin_energy(i, sold, mold);
  const double enew = single_spin_energy(i, snew, mnew);

  double logacc = -beta * (enew - eold);
  if (mod_move && jacobian_flag) logacc += 2.0 * log(mnew / mold);

  nattempt += 1.0;
  if (mod_move) mod_nattempt += 1.0;

  if ((logacc < 0.0) && (random->uniform() >= exp(logacc))) return;

  naccept += 1.0;
  if (mod_move) mod_naccept += 1.0;
  set_spin(i, snew, mnew);
}

/* ----------------------------------------------------------------------
   one move with energy global, for the local atom i on the processor that
   owns the spin and i = -1 on all others.  Every processor draws the same
   random numbers and takes the same decision; the owner moves the spin and
   contributes its single-site energies to the sums over all processors.
------------------------------------------------------------------------- */

void FixMCSpin::mc_global_spin(int i)
{
  double **sp = atom->sp;
  const int owner = (i >= 0);

  double sold[3] = {0.0, 0.0, 1.0};
  double mold = 0.0;
  if (owner) {
    sold[0] = sp[i][0];
    sold[1] = sp[i][1];
    sold[2] = sp[i][2];
    mold = sp[i][3];
  }

  double snew[3];
  propose_direction(sold, snew);

  double mnew = mold;
  int mod_move = 0;
  if (modulus_flag) {
    const double gauss = random->gaussian();
    if (owner) {
      mnew = mold + mod_sigma * gauss;
      if ((mnew < mod_min) || (mnew > mod_max))
        mnew = mold;
      else
        mod_move = 1;
    }
  }

  // [0] pair energy, [1] and [2] single-site energies before and after the
  // move, [3] logarithm of the magnitude weight, [4] magnitude move flag

  double local[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  if (owner) {
    const double so4[4] = {sold[0], sold[1], sold[2], mold};
    const double sn4[4] = {snew[0], snew[1], snew[2], mnew};
    local[1] = single_site_energy(i, so4);
    local[2] = single_site_energy(i, sn4);
    if (mod_move) local[3] = 2.0 * log(mnew / mold);
    local[4] = mod_move;
    set_spin(i, snew, mnew);
  }
  local[0] = pair_energy_local();
  double all[5];
  MPI_Allreduce(local, all, 5, MPI_DOUBLE, MPI_SUM, world);

  const double epair_new = all[0];
  const double eold = epair_cur + all[1];
  const double enew = epair_new + all[2];
  double logacc = -beta * (enew - eold);
  if ((all[4] > 0.0) && jacobian_flag) logacc += all[3];

  if (owner) {
    nattempt += 1.0;
    if (mod_move) mod_nattempt += 1.0;
  }

  if ((logacc < 0.0) && (random->uniform() >= exp(logacc))) {
    if (owner) set_spin(i, sold, mold);
    return;
  }

  if (owner) {
    naccept += 1.0;
    if (mod_move) mod_naccept += 1.0;
  }
  epair_cur = epair_new;
}

/* ----------------------------------------------------------------------
   write spin i and, in serial with energy local, its periodic images on
   this processor; otherwise forward communication updates the ghost atoms
------------------------------------------------------------------------- */

void FixMCSpin::set_spin(int i, const double *s, double smag)
{
  double **sp = atom->sp;

  sp[i][0] = s[0];
  sp[i][1] = s[1];
  sp[i][2] = s[2];
  sp[i][3] = smag;
  if (!image_flag) return;
  for (int j = image_head[i]; j >= 0; j = image_next[j]) {
    sp[j][0] = s[0];
    sp[j][1] = s[1];
    sp[j][2] = s[2];
    sp[j][3] = smag;
  }
}

/* ----------------------------------------------------------------------
   uniform direction inside a cone of half-angle cone about s
------------------------------------------------------------------------- */

void FixMCSpin::propose_direction(const double *s, double *snew)
{
  const double cosmin = cos(cone);
  const double cz = cosmin + (1.0 - cosmin) * random->uniform();
  const double sz = sqrt(MAX(0.0, 1.0 - cz * cz));
  const double phi = 2.0 * MY_PI * random->uniform();

  // orthonormal frame around s, seeded with the axis least aligned with s

  double a[3] = {0.0, 0.0, 0.0};
  if ((fabs(s[0]) <= fabs(s[1])) && (fabs(s[0]) <= fabs(s[2])))
    a[0] = 1.0;
  else if (fabs(s[1]) <= fabs(s[2]))
    a[1] = 1.0;
  else
    a[2] = 1.0;
  double e1[3] = {a[1] * s[2] - a[2] * s[1], a[2] * s[0] - a[0] * s[2], a[0] * s[1] - a[1] * s[0]};
  const double n1 = 1.0 / sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
  e1[0] *= n1;
  e1[1] *= n1;
  e1[2] *= n1;
  const double e2[3] = {s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2],
                        s[0] * e1[1] - s[1] * e1[0]};

  const double cp = cos(phi), sp = sin(phi);
  snew[0] = cz * s[0] + sz * (cp * e1[0] + sp * e2[0]);
  snew[1] = cz * s[1] + sz * (cp * e1[1] + sp * e2[1]);
  snew[2] = cz * s[2] + sz * (cp * e1[2] + sp * e2[2]);
  const double nn = 1.0 / sqrt(snew[0] * snew[0] + snew[1] * snew[1] + snew[2] * snew[2]);
  snew[0] *= nn;
  snew[1] *= nn;
  snew[2] *= nn;
}

/* ----------------------------------------------------------------------
   energy of the bonds of spin i when it has the value s (direction and
   magnitude).  The energy functions return half a bond energy, which the
   pair styles subtract once from each end of a bond.  With trial = 1, the
   result is the energy that changes when spin i changes: the bond to
   another spin counts twice, since its other half changes as well, and the
   bond to a periodic image of spin i counts once, with the image taking
   the value s.  With trial = 0, it is the share of spin i in the total
   energy, where every bond counts once.
------------------------------------------------------------------------- */

double FixMCSpin::bond_energy(int i, const double *s, int trial)
{
  if (!list) return 0.0;

  double **x = atom->x;
  double **sp = atom->sp;
  int *type = atom->type;
  tagint *tag = atom->tag;
  const int ntp1 = atom->ntypes + 1;
  const int itype = type[i];

  double spi[4] = {s[0], s[1], s[2], s[3]};
  double energy = 0.0;

  int *jlist = list->firstneigh[i];
  const int jnum = list->numneigh[i];
  for (int jj = 0; jj < jnum; jj++) {
    const int j = jlist[jj] & NEIGHMASK;
    const double delx = x[j][0] - x[i][0];
    const double dely = x[j][1] - x[i][1];
    const double delz = x[j][2] - x[i][2];
    const double rsq = delx * delx + dely * dely + delz * delz;
    if (rsq == 0.0) continue;
    const double rinv = 1.0 / sqrt(rsq);
    double eij[3] = {-delx * rinv, -dely * rinv, -delz * rinv};

    const int self = (tag[j] == tag[i]);
    double spj[3] = {sp[j][0], sp[j][1], sp[j][2]};
    if (self) {
      spj[0] = spi[0];
      spj[1] = spi[1];
      spj[2] = spi[2];
    }
    const double weight = (trial && !self) ? 2.0 : 1.0;
    const int ij = itype * ntp1 + type[j];

    for (std::size_t k = 0; k < spin_pairs.size(); k++) {
      if (rsq > pair_cutsq[k][ij]) continue;
      if (pair_kind[k] == EXCHANGE)
        energy -= weight *
            static_cast<PairSpinExchange *>(spin_pairs[k])->compute_energy(i, j, rsq, spi, spj);
      else if (pair_kind[k] == BIQUADRATIC)
        energy -= weight *
            static_cast<PairSpinExchangeBiquadratic *>(spin_pairs[k])
                ->compute_energy(i, j, rsq, spi, spj);
      else if (pair_kind[k] == NEEL)
        energy -= weight *
            static_cast<PairSpinNeel *>(spin_pairs[k])
                ->compute_neel_energy(i, j, rsq, eij, spi, spj);
    }
  }

  return energy;
}

/* ----------------------------------------------------------------------
   energy that changes with spin i, for direction s and magnitude smag
------------------------------------------------------------------------- */

double FixMCSpin::single_spin_energy(int i, const double *s, double smag)
{
  const double spit[4] = {s[0], s[1], s[2], smag};
  return bond_energy(i, spit, 1) + single_site_energy(i, spit);
}

/* ----------------------------------------------------------------------
   energy of spin i alone: fix precession/spin and the fixes acting on the
   spin magnitude
------------------------------------------------------------------------- */

double FixMCSpin::single_site_energy(int i, const double *s)
{
  double spit[4] = {s[0], s[1], s[2], s[3]};
  double energy = 0.0;
  int *mask = atom->mask;

  for (auto &ifix : precession_fixes) {
    if (!(mask[i] & ifix->groupbit)) continue;
    if (ifix->zeeman_flag) energy -= ifix->compute_zeeman_energy(spit);
    if (ifix->aniso_flag) energy -= ifix->compute_anisotropy_energy(spit);
    if (ifix->cubic_flag) energy -= ifix->compute_cubic_energy(spit);
    if (ifix->hexaniso_flag) energy -= ifix->compute_hexaniso_energy(spit);
  }
  for (auto &ifix : modulus_fixes) energy += ifix->modulus_energy(i, spit[3]);

  return energy;
}

/* ----------------------------------------------------------------------
   magnetic energy of the group: the spin pair styles and single-site terms
   (energy local), or the pair style and single-site terms (energy global)
------------------------------------------------------------------------- */

double FixMCSpin::magnetic_energy()
{
  double **sp = atom->sp;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;

  double local = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    local += single_site_energy(i, sp[i]);
    if (!energy_global) local += bond_energy(i, sp[i], 0);
  }

  double all = 0.0;
  MPI_Allreduce(&local, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  if (energy_global) all += epair_cur;
  return all;
}

/* ----------------------------------------------------------------------
   pair energy of this processor from a full re-run of the pair style with
   the current spins; all processors must call it together
------------------------------------------------------------------------- */

double FixMCSpin::pair_energy_local()
{
  comm->forward_comm();
  clear_forces();
  force->pair->compute(1, 0);
  return force->pair->eng_vdwl + force->pair->eng_coul;
}

/* ----------------------------------------------------------------------
   pair energy of the whole system with the current spins
------------------------------------------------------------------------- */

double FixMCSpin::model_energy()
{
  const double local = pair_energy_local();
  double all = 0.0;
  MPI_Allreduce(&local, &all, 1, MPI_DOUBLE, MPI_SUM, world);
  return all;
}

/* ---------------------------------------------------------------------- */

void FixMCSpin::clear_forces()
{
  const size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) {
    memset(&atom->f[0][0], 0, 3 * nbytes);
    if (atom->avec->forceclearflag) atom->avec->force_clear(0, nbytes);
  }
}

/* ----------------------------------------------------------------------
   acceptance rates of the last sweeps, and the adaptation of the cone
------------------------------------------------------------------------- */

void FixMCSpin::mc_diagnostics()
{
  double local[4] = {naccept, nattempt, mod_naccept, mod_nattempt};
  double all[4];
  MPI_Allreduce(local, all, 4, MPI_DOUBLE, MPI_SUM, world);

  accept_rate = (all[1] > 0.0) ? all[0] / all[1] : 0.0;
  mod_accept_rate = (all[3] > 0.0) ? all[2] / all[3] : 0.0;

  if (adapt_flag) {
    cone *= pow(MAX(accept_rate, 1.0e-3) / MAX(cone_target, 1.0e-3), 0.1);
    cone = MIN(cone, MY_PI);
    cone = MAX(cone, 1.0e-4);
  }
}

/* ----------------------------------------------------------------------
   1 = acceptance rate of the last sweeps
   2 = current cone half-angle (degrees)
   3 = acceptance rate of the moves that changed the spin magnitude
   4 = magnetic energy of the group
------------------------------------------------------------------------- */

double FixMCSpin::compute_vector(int n)
{
  if (n == 0) return accept_rate;
  if (n == 1) return cone * 180.0 / MY_PI;
  if (n == 2) return mod_accept_rate;
  return magnetic_energy();
}

/* ----------------------------------------------------------------------
   divide the subdomain into 8 sectors whose half-widths exceed the spin
   cutoff, so that spins of one sector on different processors never
   interact; as in fix nve/spin
------------------------------------------------------------------------- */

void FixMCSpin::sectoring()
{
  const double rv = rsec[0];
  double rs[3], rcut[3];
  if (domain->triclinic) {
    const double *h_inv = domain->h_inv;
    for (int d = 0; d < 3; d++) rs[d] = domain->subhi_lamda[d] - domain->sublo_lamda[d];
    rcut[0] = rv * sqrt(h_inv[0] * h_inv[0] + h_inv[5] * h_inv[5] + h_inv[4] * h_inv[4]);
    rcut[1] = rv * sqrt(h_inv[1] * h_inv[1] + h_inv[3] * h_inv[3]);
    rcut[2] = rv * h_inv[2];
  } else {
    for (int d = 0; d < 3; d++) {
      rs[d] = domain->subhi[d] - domain->sublo[d];
      rcut[d] = rv;
    }
  }

  int sec[3];
  for (int d = 0; d < 3; d++) sec[d] = (rs[d] >= 2.0 * rcut[d]) ? 2 : 1;
  nsectors = sec[0] * sec[1] * sec[2];
  if (nsectors != 8)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix mc/spin needs each subdomain to be at least twice as wide as the spin "
               "cutoff plus the neighbor skin in every direction for the parallel sweeps");
  for (int d = 0; d < 3; d++) rsec[d] = rs[d] / sec[d];
}

/* ---------------------------------------------------------------------- */

int FixMCSpin::coords2sector(double *x)
{
  double lamda[3];
  double *xs = x;
  double *sublo = domain->sublo;
  if (domain->triclinic) {
    domain->x2lamda(x, lamda);
    xs = lamda;
    sublo = domain->sublo_lamda;
  }
  const int seci0 = xs[0] > (sublo[0] + rsec[0]);
  const int seci1 = xs[1] > (sublo[1] + rsec[1]);
  const int seci2 = xs[2] > (sublo[2] + rsec[2]);
  return seci0 + 2 * seci1 + 4 * seci2;
}

/* ----------------------------------------------------------------------
   write the cone angle and the random number states of all processors
   layout: [0] = number of processors, [1] = seed, [2] = cone,
   [3] = energy global flag, then the states
------------------------------------------------------------------------- */

void FixMCSpin::write_restart(FILE *fp)
{
  const int nprocs = comm->nprocs;
  const int nsize = 4 + RanMars::STATE_SIZE * nprocs;

  std::vector<double> list_out;
  if (comm->me == 0) {
    list_out.resize(nsize);
    list_out[0] = nprocs;
    list_out[1] = seed;
    list_out[2] = cone;
    list_out[3] = energy_global;
  }

  double state[RanMars::STATE_SIZE];
  random->get_state(state);
  MPI_Gather(state, RanMars::STATE_SIZE, MPI_DOUBLE, (comm->me == 0) ? &list_out[4] : nullptr,
             RanMars::STATE_SIZE, MPI_DOUBLE, 0, world);

  if (comm->me == 0) {
    int size = nsize * sizeof(double);
    fwrite(&size, sizeof(int), 1, fp);
    fwrite(list_out.data(), sizeof(double), nsize, fp);
  }
}

/* ----------------------------------------------------------------------
   with energy global, all processors continue with the state of processor
   0, which is the same for any number of processors
------------------------------------------------------------------------- */

void FixMCSpin::restart(char *buf)
{
  auto *list_in = (double *) buf;
  cone = list_in[2];

  const int same = ((int) list_in[1] == seed) && ((int) list_in[3] == energy_global);
  if (same && energy_global) {
    random->set_state(list_in + 4);
  } else if (same && ((int) list_in[0] == comm->nprocs)) {
    const int stride = RanMars::state_size(list_in + 4);
    random->set_state(list_in + 4 + (size_t) comm->me * stride);
  } else if (comm->me == 0) {
    error->warning(FLERR,
                   "Fix mc/spin cannot restore the random number state from the restart file: "
                   "the seed, the energy keyword, or with energy local the number of "
                   "processors differs");
  }
}

/* ---------------------------------------------------------------------- */

void *FixMCSpin::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "lattice_flag") == 0) return &lattice_flag;
  return nullptr;
}
