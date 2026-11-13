.. index:: fix nve/spin/sib
.. index:: fix nve/spin/pc

fix nve/spin/sib command
========================

fix nve/spin/pc command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID style keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *nve/spin/sib* or *nve/spin/pc*
* zero or more keyword/value pairs may be appended
* keyword = *lattice* or *field* or *maxiter* or *tol*

  .. parsed-literal::

       *lattice* value = *moving* or *frozen*
         moving = integrate the spins and the atomic positions
         frozen = integrate the spins only
       *field* value = *auto* or *perspin* or *global*
         auto = same as perspin
         perspin = evaluate the magnetic field of each spin separately
         global = evaluate the magnetic field by re-running the pair style
       *maxiter* value = N (*nve/spin/sib* only)
         N = largest number of field evaluations per timestep
       *tol* value = delta (*nve/spin/sib* only)
         delta = convergence criterion for the change of the spin directions

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/spin/sib lattice moving
   fix 1 all nve/spin/sib lattice frozen maxiter 4 field global

   fix 1 all precession/spin zeeman 1.0 0.0 0.0 1.0
   fix 2 all langevin/spin 300.0 0.1 48279
   fix 3 all nve/spin/pc lattice moving

Description
"""""""""""

.. versionadded:: TBD

These fixes integrate the precession of the magnetic spins together
with the motion of the atoms, like :doc:`fix nve/spin <fix_nve_spin>`.
The spin directions :math:`\vec{s}_i` obey the Landau-Lifshitz-Gilbert
equation

.. math::

   \frac{d\vec{s}_i}{dt} = \vec{\omega}_i \times \vec{s}_i, \qquad
   \vec{\omega}_i = \frac{1}{1+\alpha^2}\left(\vec{H}_i - \alpha\,
   \vec{H}_i\times\vec{s}_i\right),

where :math:`\vec{H}_i` is the magnetic precession field of spin *i*
(in rad/ps, the *fm* quantity of the SPIN package) and :math:`\alpha` is
the transverse damping.  The spin magnitudes are not changed.

:doc:`Fix nve/spin <fix_nve_spin>` advances the spins one after the
other, each spin seeing the spins that have already been moved.  The
fixes on this page instead advance all spins of the group together: at
every stage of a timestep, the magnetic fields of all spins are evaluated
for the same configuration, and then all spins are updated.  The result
does not depend on the order of the atoms or on the number of
processors.  Each timestep needs a fixed number of field evaluations,
independent of the number of atoms, which makes these fixes suitable for
magnetic potentials whose field has to be computed for the whole system
at once, such as machine-learning potentials.

Fix *nve/spin/sib* uses the semi-implicit midpoint method
:ref:`(Mentink) <Mentink2010>`.  Each spin is rotated about the
precession vector evaluated at the midpoint between its old and its new
direction, with a rotation that keeps the spin magnitude exactly.  The
new directions are found by repeated iteration, which stops when no
component of any spin direction changes by more than *delta* (keyword
*tol*), or after *N* field evaluations (keyword *maxiter*).
The converged step is symmetric in time and shows no long-time drift of
the energy of an undamped system.  A small value of *maxiter*, such as 2
or 3, reduces the cost for expensive potentials at the price of a less
accurate step.

Fix *nve/spin/pc* uses the predictor-corrector (Heun) method
:ref:`(Garcia-Palacios) <GarciaPalacios1998>` with two field
evaluations per timestep.  This explicit method does not conserve the
energy of an undamped system, whose energy grows slowly.  It is intended
for simulations with a spin bath, for example for sampling at finite
temperature.

When the atoms move (*lattice moving*), their positions and velocities
are integrated with velocity Verlet.  The spins are advanced once per
timestep, at the new atomic positions, and the forces of the following
velocity update are computed with the new positions and spins.

To add damping and a thermal field, define :doc:`fix langevin/spin
<fix_langevin_spin>` after any :doc:`fix precession/spin
<fix_precession_spin>`: these fixes take the transverse damping, the
temperature, and the random number seed from it.  Only spins in the
groups of both fixes are damped and thermostatted.  The
thermal field, a random field with a variance of
:math:`2\alpha k_B T/(\hbar\,\Delta t)` per component, is drawn once per
timestep and kept fixed within the timestep.  This corresponds to the
Stratonovich interpretation of the stochastic equation of motion and
samples the canonical distribution of the spin directions.  Without
fix langevin/spin, the spins precess without damping.

The *field* keyword selects how the magnetic field is evaluated at the
intermediate configurations of a timestep.  With *perspin* (the default),
the field of each spin is computed from the spin pair styles, :doc:`fix
precession/spin <fix_precession_spin>`, and :doc:`fix setforce/spin
<fix_setforce_spin>`.  With *global*, the pair style is run as a whole
for every field evaluation, which is needed for magnetic pair styles
that do not belong to the SPIN package and cannot compute the field of
one spin alone, such as *pair_style deepspin* of DeePMD-kit.  In this
case, the cost of a timestep is dominated by the field evaluations: two
for *nve/spin/pc*, and up to *maxiter* for *nve/spin/sib*.

At the first timestep of each run, the magnetic field is rebuilt in the
selected way and compared with the field of the force computation.  If
a pair style, fix, or long-range solver contributes to the field in a way
that these fixes cannot re-evaluate, the run stops with an error, rather
than silently integrating a different magnetic energy.  Only the field
components perpendicular to the spins are compared, since only they act
on the spin directions.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These fixes write the state of the random number generators of the
thermal field to :doc:`binary restart files <restart>`.  A restarted
simulation with the same number of processors and the same random number
seed continues with the same sequence of random numbers; see the
Restrictions section for the other conditions of an exact restart.

None of the :doc:`fix_modify <fix_modify>` options are relevant to these
fixes.

These fixes compute a global vector, which can be accessed by various
:doc:`output commands <Howto_output>`.  The vector values are
"intensive".  Fix *nve/spin/pc* computes 3 values and fix *nve/spin/sib*
computes 6 values:

  #. largest deviation of the rebuilt field at the start of the run (rad/ps)
  #. largest field component perpendicular to a spin at that time (rad/ps)
  #. difference between two evaluations of the pair style with *field global* (rad/ps)
  #. average number of field evaluations per timestep in this run
  #. fraction of the timesteps of this run that ended without convergence
  #. largest final change of a spin direction component in this run

These fixes are not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

These fixes are part of the SPIN package.  They are only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

These fixes require :doc:`atom_style spin <atom_style>`.  They cannot be
used with :doc:`run_style respa <run_style>`, with another spin
integrator such as :doc:`fix nve/spin <fix_nve_spin>`, or with the
inertial spin dynamics commands such as :doc:`fix nve/tspin
<fix_nve_tspin>`.  Only one :doc:`fix langevin/spin <fix_langevin_spin>`
and one :doc:`fix setforce/spin <fix_setforce_spin>` are supported.

The magnetic field of the long-range dipolar solvers *ewald/dipole/spin*
and *pppm/dipole/spin* cannot be re-evaluated by these fixes.

With a spin bath, each processor draws its own random numbers for its
atoms in the order in which it stores them, so the trajectories depend
on the number of processors and on the sorting of the atoms; use
:doc:`atom_modify sort 0 0.0 <atom_modify>` for trajectories that can be
reproduced exactly, for example across a restart.

Related commands
""""""""""""""""

:doc:`fix nve/spin <fix_nve_spin>`,
:doc:`fix langevin/spin <fix_langevin_spin>`,
:doc:`fix precession/spin <fix_precession_spin>`

Default
"""""""

The keyword defaults are lattice = moving, field = auto, maxiter = 10,
and tol = 1.0e-12.

----------

.. _Mentink2010:

**(Mentink)** J. H. Mentink, M. V. Tretyakov, A. Fasolino,
M. I. Katsnelson, and Th. Rasing, J. Phys.: Condens. Matter, 22, 176001
(2010).

.. _GarciaPalacios1998:

**(Garcia-Palacios)** J. L. Garcia-Palacios and F. J. Lazaro,
Phys. Rev. B, 58, 14937 (1998).
