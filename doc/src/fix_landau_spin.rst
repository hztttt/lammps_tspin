.. index:: fix landau/spin

fix landau/spin command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID landau/spin A B keyword value

* ID, group-ID are documented in :doc:`fix <fix>` command
* landau/spin = style name of this fix command
* A = coefficient of the second-order term (energy per squared Bohr
  magneton, :math:`\mathrm{eV}/\mu_B^2` in metal units)
* B = coefficient of the fourth-order term (energy per Bohr magneton to
  the fourth power, :math:`\mathrm{eV}/\mu_B^4` in metal units)
* zero or more keyword/value pairs may be appended
* keyword = *c*

  .. parsed-literal::

       *c* value = C
         C = coefficient of the sixth-order term (energy per Bohr magneton to
             the sixth power, :math:`\mathrm{eV}/\mu_B^6` in metal units)

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all landau/spin -0.44 0.08
   fix_modify 1 energy yes

   group fe type 1
   fix 2 fe landau/spin -0.44 0.08 c 0.001

Description
"""""""""""

.. versionadded:: TBD

Apply a single-site Landau energy to the magnitude of each magnetic
spin in the group,

.. math::

   E_{L,i} = A|\vec{S}_i|^2 + B|\vec{S}_i|^4 + C|\vec{S}_i|^6,

where :math:`|\vec{S}_i|` is the spin magnitude of atom *i* in Bohr
magnetons.  This is the longitudinal part of a Landau-Heisenberg
Hamiltonian :ref:`(Ma) <Ma2012>`: combined with a magnetic pair
interaction, it gives the spin magnitude an energy of its own, so that
the magnitude can fluctuate and respond to temperature, magnetic order,
and the local atomic environment.  With :math:`A < 0` and :math:`B > 0`
the energy has its minimum at the finite magnitude
:math:`|\vec{S}| = \sqrt{-A/(2B)}` (for :math:`C = 0`), with
:math:`A > 0` at zero magnitude.

The energy depends only on the spin magnitude, so the resulting
magnetic force is purely radial,

.. math::

   \vec{F}^{m}_{i,\mathrm{L}} = -\frac{\partial E_{L,i}}{\partial\vec{S}_i}
   = -\left(2A|\vec{S}_i| + 4B|\vec{S}_i|^3 + 6C|\vec{S}_i|^5\right)
   \hat{s}_i,

and exerts no torque on the spin.  It is stored in the same form as
for the :doc:`inertial spin dynamics (TSPIN) commands <fix_nve_tspin>`,
:math:`\vec{fm}_i = (|\vec{S}_i|/\hbar)\,\vec{F}^{m}_{i,\mathrm{L}}`,
so that it can be combined with integrators that treat the spin
magnitude as a dynamical variable, such as :doc:`fix nve/tspin
<fix_nve_tspin>`, :doc:`fix nvt/tspin <fix_nvt_tspin>`, and :doc:`fix
nve/spin/llb <fix_nve_spin_llb>`.  With
integrators that keep the spin magnitude fixed, such as :doc:`fix
nve/spin <fix_nve_spin>`, the Landau energy is a constant and the fix
has no effect on the dynamics.

The energy must not decrease without limit at large spin magnitude,
otherwise the magnitude grows indefinitely: *C* must not be negative,
*B* must be positive or zero when *C* is zero, and *A* must be positive
when both *B* and *C* are zero.

The coefficients apply to all atoms in the group.  To use different
coefficients for different atom types, define one group and one fix
per atom type, as in the second example above.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

Enable :doc:`fix_modify energy yes <fix_modify>` to include the Landau
energy in the global potential energy *pe*; it is excluded by default.
Since the energy changes when the spin magnitudes change, it should be
included whenever the magnitudes evolve.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix.  This allows one to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its forces.  Default is the
outermost level.

This fix computes a global scalar containing the total Landau energy of
the group, accessible as *f_ID*.  This scalar is "extensive" and has
energy units.

This fix can be used during :doc:`energy minimization <minimize>`.
The spin minimizers do not change the spin magnitudes, so only the
energy contribution of this fix matters there.

Restrictions
""""""""""""

This fix requires :doc:`atom_style spin <atom_style>` and the SPIN
package; see :doc:`Build package <Build_package>`.

Related commands
""""""""""""""""

:doc:`fix spring/tspin <fix_spring_tspin>`,
:doc:`fix nve/spin/llb <fix_nve_spin_llb>`,
:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`

Default
""""""""

The default value of the optional keyword is C = 0.0.

----------

.. _Ma2012:

**(Ma)** P.-W. Ma and S. L. Dudarev, Longitudinal magnetic fluctuations in
Langevin spin dynamics, Phys. Rev. B, 86, 054416 (2012).
