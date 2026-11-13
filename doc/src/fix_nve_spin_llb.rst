.. index:: fix nve/spin/llb

fix nve/spin/llb command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID nve/spin/llb scheme keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* nve/spin/llb = style name of this fix command
* scheme = *sib* or *pc*

  .. parsed-literal::

       *sib* = semi-implicit midpoint method
       *pc* = predictor-corrector method

* zero or more keyword/value pairs may be appended
* keyword = *lattice* or *field* or *alpha_perp* or *alpha_par* or *maxiter* or *tol*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate the spins and the atomic positions
         frozen = integrate the spins only
       *field* value = *auto* or *perspin* or *global*
         auto = same as perspin
         perspin = evaluate the magnetic field of each spin separately
         global = evaluate the magnetic field by re-running the pair style
       *alpha_perp* value = transverse damping (unitless), only without fix langevin/spin/llb
       *alpha_par* value = longitudinal damping (unitless), only without fix langevin/spin/llb
       *maxiter* value = N
         N = largest number of field evaluations of one implicit step (scheme *sib*)
       *tol* value = delta
         delta = convergence criterion of an implicit step (scheme *sib*)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all landau/spin -0.44 0.08
   fix 2 all langevin/spin/llb temp 300.0 alpha_perp 0.1 alpha_par 0.1 seed 48279
   fix 3 all nve/spin/llb pc lattice moving

   fix 3 all nve/spin/llb sib lattice frozen alpha_perp 0.1 alpha_par 0.05

Description
"""""""""""

.. versionadded:: TBD

This fix integrates the stochastic Landau-Lifshitz-Bloch (LLB) equation
of atomic spins whose magnitude is a dynamical variable, together with
the motion of the atoms.  Writing each spin as
:math:`\vec{S}_i = |\vec{S}_i|\,\vec{s}_i` with a unit direction
:math:`\vec{s}_i`, the dynamics has two channels:

* The direction obeys the Landau-Lifshitz-Gilbert equation with the
  transverse damping :math:`\alpha_\perp`, as integrated by :doc:`fix
  nve/spin/sib and nve/spin/pc <fix_nve_spin_sib>`.
* The magnitude relaxes towards lower energy and fluctuates thermally,
  following the Langevin equation of :math:`u_i = \ln|\vec{S}_i|`,

  .. math::

     du_i = \alpha_\parallel \left( |\vec{S}_i| L_i
     + \frac{3 k_B T}{\hbar} \right) dt
     + \sqrt{\frac{2\alpha_\parallel k_B T\,dt}{\hbar}}\,\xi_i,
     \qquad L_i = -\frac{1}{\hbar}\frac{\partial E}{\partial|\vec{S}_i|},

  with the longitudinal damping :math:`\alpha_\parallel` and a Gaussian
  random number :math:`\xi_i`.

At temperature *T*, the spins sample the distribution
:math:`\exp(-E/k_B T)\,d^3\vec{S}_1 \cdots d^3\vec{S}_N` of the energy *E*
that the pair styles and fixes define; the term :math:`3k_B T/\hbar` is
the weight :math:`|\vec{S}|^2` of the spin magnitude in :math:`d^3\vec{S}`.
No parameter depends on the temperature: the temperature dependence of
the magnetization, of the spin magnitudes, and of the longitudinal
response follows from the energy model, as in the atomistic approach of
:ref:`(Ma and Dudarev) <MaDudarev2012>`, rather than from the
temperature-dependent parameters of the macroscopic LLB equation
:ref:`(Garanin) <Garanin1997>`.  Unlike the inertial spin dynamics of
:doc:`fix nve/tspin <fix_nve_tspin>`, the magnitude has no mass and no
velocity.  With :math:`\alpha_\parallel = 0` the magnitudes stay fixed,
and the fix integrates the same equations as fix nve/spin/sib or
nve/spin/pc.

The longitudinal force *L* is collected from the terms whose energy
depends on the spin magnitudes:

* :doc:`fix landau/spin <fix_landau_spin>` and :doc:`fix spring/tspin
  <fix_spring_tspin>`,
* the Zeeman term of :doc:`fix precession/spin <fix_precession_spin>`;
  its anisotropy terms depend on the spin direction only,
* spin pair styles that provide the derivative of the energy with
  respect to the full spin vector, such as :doc:`pair_style
  spin/dipole/cut <pair_spin_dipole>`; the energies of the other spin pair
  styles depend on the spin directions only and contribute nothing,
* with *field global*, a pair style that is not part of the SPIN package,
  such as *pair_style deepspin* of DeePMD-kit.  Its magnetic force must
  follow the convention of :doc:`fix nve/tspin <fix_nve_tspin>`, i.e.
  include the derivative with respect to the spin magnitude.

Without a term that keeps the spin magnitudes finite, such as fix
landau/spin, the magnitudes grow without limit at finite temperature.
Spins with zero magnitude keep it, since the logarithm of the magnitude
is undefined; start them from a finite magnitude.

Each timestep consists of half a longitudinal step, a full step of the
directions, and half a longitudinal step.  With scheme *sib*, both
channels use the implicit midpoint method, solved by repeated iteration
until no change exceeds *delta* or *N* field evaluations have been done;
with scheme *pc*, both use the predictor-corrector (Heun) method.  The
thermal noise is drawn once per step and kept fixed during the
iterations.  When the atoms move, their positions are integrated with
velocity Verlet as in fix nve/spin/sib.  The *field* keyword works as for
:doc:`fix nve/spin/sib <fix_nve_spin_sib>`.  Each longitudinal half step
evaluates the longitudinal force once per iteration of the implicit step
with scheme *sib*, and twice with scheme *pc*; with *field global* and a
pair style whose energy depends on the spin magnitudes, each evaluation
re-runs the pair style.

The finite timestep slightly biases the sampled distributions, by an
amount that grows with the product of damping and timestep.  For a
single spin in a Landau energy with its minimum at 1.74 Bohr magnetons, at
300 K with :math:`\alpha_\parallel = 0.2` and a timestep of 0.1 fs, the
average spin magnitude was 0.06 percent too large with scheme *sib*
(0.04 percent at 0.05 fs) and within 0.01 percent with scheme *pc*.
Verify that the quantities of interest do not change when the timestep
is reduced.

The damping constants and the temperature are taken from :doc:`fix
langevin/spin/llb <fix_langevin_spin_llb>`, which must use the same group
as this fix.  Without it, the damping constants can be given with the
*alpha_perp* and *alpha_par* keywords for damped dynamics without thermal
noise; both default to zero.

At the first timestep of each run, the magnetic field is rebuilt in the
selected way and compared with the field of the force computation,
including the components along the spins, as described for :doc:`fix
nve/spin/sib <fix_nve_spin_sib>`.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the state of the random number generators of the
thermal noise to :doc:`binary restart files <restart>`, with the same
conditions for an exact restart as fix nve/spin/sib.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes the same global vector of 6 "intensive" values as
:doc:`fix nve/spin/sib <fix_nve_spin_sib>`, where the implicit steps of
both channels are counted.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the SPIN package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

This fix requires :doc:`atom_style spin <atom_style>`.  It cannot be used
with :doc:`run_style respa <run_style>`, with another spin integrator,
with the inertial spin dynamics commands, with :doc:`fix setforce/spin
<fix_setforce_spin>`, or with :doc:`fix langevin/spin
<fix_langevin_spin>`.  With *field global*, the spin pair styles must
either all provide the derivative with respect to the full spin vector
or none of them.

With thermal noise, each processor draws its own random numbers for its
atoms in the order in which it stores them, so the trajectories depend on
the number of processors and on the sorting of the atoms; use
:doc:`atom_modify sort 0 0.0 <atom_modify>` for trajectories that can be
reproduced exactly.

Related commands
""""""""""""""""

:doc:`fix langevin/spin/llb <fix_langevin_spin_llb>`,
:doc:`fix nve/spin/sib <fix_nve_spin_sib>`,
:doc:`fix landau/spin <fix_landau_spin>`,
:doc:`fix nve/tspin <fix_nve_tspin>`

Default
"""""""

The keyword defaults are lattice = moving, field = auto, alpha_perp = 0.0,
alpha_par = 0.0, maxiter = 10, and tol = 1.0e-12.

----------

.. _MaDudarev2012:

**(Ma and Dudarev)** P.-W. Ma and S. L. Dudarev, Phys. Rev. B, 86,
054416 (2012).

.. _Garanin1997:

**(Garanin)** D. A. Garanin, Phys. Rev. B, 55, 3050 (1997).
