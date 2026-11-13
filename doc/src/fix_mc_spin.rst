.. index:: fix mc/spin

fix mc/spin command
===================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID mc/spin T seed keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mc/spin = style name of this fix command
* T = temperature (temperature units)
* seed = random number seed (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *cone* or *adapt* or *target* or *nsweep* or *every* or *modulus* or *modsigma* or *modmin* or *modmax* or *jacobian* or *energy*

  .. parsed-literal::

       *cone* value = angle
         angle = initial half-angle of the cone of trial directions (degrees)
       *adapt* value = *yes* or *no*
         yes = adjust the cone angle after every set of sweeps
         no = keep the cone angle fixed
       *target* value = acceptance rate the cone adjustment aims at (between 0 and 1)
       *nsweep* value = number of sweeps per invocation
       *every* value = N
         N = perform the sweeps every N timesteps
       *modulus* value = *off* or *walk*
         off = keep the spin magnitudes fixed
         walk = also change the spin magnitude in each trial move
       *modsigma* value = width of the trial change of the spin magnitude (Bohr magnetons)
       *modmin*, *modmax* values = smallest and largest allowed spin magnitude (Bohr magnetons)
       *jacobian* value = *yes* or *no*
         yes = include the weight of the spin magnitude in the acceptance probability
         no = omit it
       *energy* value = *local* or *global*
         local = evaluate the energy change from the bonds of the moved spin
         global = evaluate the energy change by re-running the pair style

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all mc/spin 300.0 48279
   fix 1 all mc/spin 300.0 48279 cone 20 adapt no nsweep 2

   fix 1 all landau/spin -0.44 0.08
   fix 2 all mc/spin 300.0 48279 modulus walk modsigma 0.1
   fix 3 all nve

   fix 1 all mc/spin 600.0 48279 energy global every 10

Description
"""""""""""

.. versionadded:: TBD

Sample the atomic spins of the group at temperature *T* with the
Metropolis Monte Carlo method.  Each sweep attempts one trial move per
spin, in a random order.  A move rotates the spin to a direction drawn
uniformly inside a cone around its current direction; with *modulus
walk*, it also changes the spin magnitude by a Gaussian random number of
standard deviation *modsigma*, and moves that would leave the interval
from *modmin* to *modmax* keep the old magnitude.  Both proposals are
symmetric, so a move is accepted with the probability
:math:`\min[1, \exp(-\Delta E/k_B T)\,(|\vec{S}_{new}|/|\vec{S}_{old}|)^2]`,
where the last factor, included with *jacobian yes*, is the weight of the
spin magnitude in :math:`d^3\vec{S}`.  The spins thus sample
:math:`\exp(-E/k_B T)` on the sphere of their fixed magnitude, or in
three dimensions with *modulus walk*.  No time step, mass, or damping
constant enters, which makes this fix a reference for the equilibrium
averages of spin dynamics.

With *adapt yes*, the cone half-angle is scaled after every set of
sweeps towards the acceptance rate *target* :ref:`(Alzate-Cardona)
<AlzateCardona2019>`.  An adapted proposal does not strictly preserve the
equilibrium distribution, so the cone should be fixed with *adapt no*
for the production part of a simulation, for example by redefining the
fix after equilibration with the angle reported in the output.

The fix moves only the spins.  The sweeps are performed every *every*
timesteps after the atoms have moved and before the forces are
computed, so the forces always correspond to the current spins.  The
atoms can be integrated by another fix, such as :doc:`fix nve <fix_nve>`
or :doc:`fix nvt <fix_nh>`, for a coupled sampling of spins and atoms;
without such a fix the lattice is frozen.  When such a fix integrates
the atom positions, the spin pair styles compute their magneto-elastic
forces.

With *energy local*, the energy change of a move is computed from the
bonds of the moved spin, for the spin pair styles :doc:`spin/exchange
<pair_spin_exchange>`, *spin/exchange/biquadratic*, and :doc:`spin/neel
<pair_spin_neel>`, plus its single-site energies from :doc:`fix
precession/spin <fix_precession_spin>` (Zeeman and anisotropy terms),
:doc:`fix landau/spin <fix_landau_spin>`, and :doc:`fix spring/tspin
<fix_spring_tspin>`.  The bonds use the coefficients and cutoff of each
pair of atom types; with :doc:`pair_style hybrid <pair_hybrid>` or
*hybrid/overlay*, a spin pair style contributes only for the pairs of
atom types assigned to it.  In parallel, every subdomain is divided into
sectors, as for :doc:`fix nve/spin <fix_nve_spin>`, so that spins that
are moved at the same time on different processors do not interact.
With *energy global*, every move is scored by re-running the whole pair
style, which is required for magnetic potentials without a local energy,
such as *pair_style deepspin* of DeePMD-kit.  Each sweep then costs one
evaluation of the pair style per spin.  In parallel, the moves are still
attempted one at a time for the whole system, and all processors share
each evaluation of the pair style.  All processors use the same random
numbers, so that a run follows the same sequence of moves on any number
of processors; small round-off differences of the energies can change
the outcome of a move only very rarely.

At the start of each run, the energy used in the Metropolis test is
compared with the energy computed by the pair styles and fixes, and the
run stops if they differ.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the current cone angle and the state of the random
number generators to :doc:`binary restart files <restart>`.  A restarted
simulation with the same seed and *energy* setting continues with the
same random numbers; with *energy local*, it must also run on the same
number of processors.

None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.

This fix computes a global vector of length 4, which can be accessed by
various :doc:`output commands <Howto_output>`.  The vector values are
"intensive":

  #. acceptance rate of the last set of sweeps
  #. current cone half-angle (degrees)
  #. acceptance rate of the moves that changed the spin magnitude
  #. energy sampled by the Metropolis test (energy units)

With *energy local*, the fourth value is the energy of the spin pair
styles and the single-site terms of the spins in the group.  With
*energy global*, it is the energy of the whole pair style, including
any non-magnetic terms such as an embedded-atom potential, plus the
single-site terms of the group.

This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the SPIN package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

This fix requires :doc:`atom_style spin <atom_style>` and atom IDs.  It
cannot be used together with a spin integrator, a spin thermostat,
:doc:`fix setforce/spin <fix_setforce_spin>`, or the spin transfer
torque of fix precession/spin.  The energies of other fixes that act on
the spins are not included.

*energy local* supports the pair styles listed above.  With spin pair
styles in parallel, each subdomain must be at least twice as wide as the
largest spin cutoff plus the neighbor skin in every direction.  With
*energy global* and *every* larger than 1, the fourth vector value is
the energy of the last sweep.

With *energy local*, each processor draws its own random numbers, so the
trajectories depend on the number of processors.

Related commands
""""""""""""""""

:doc:`fix landau/spin <fix_landau_spin>`,
:doc:`fix nve/spin/llb <fix_nve_spin_llb>`,
:doc:`fix nve/spin <fix_nve_spin>`

Default
"""""""

The keyword defaults are cone = 30, adapt = yes, target = 0.5, nsweep = 1,
every = 1, modulus = off, modsigma = 0.05, modmin = 0.1, modmax = 5.0,
jacobian = yes, and energy = local.

----------

.. _AlzateCardona2019:

**(Alzate-Cardona)** J. D. Alzate-Cardona, D. Sabogal-Suarez,
R. F. L. Evans, and E. Restrepo-Parra, J. Phys.: Condens. Matter, 31,
095802 (2019).
