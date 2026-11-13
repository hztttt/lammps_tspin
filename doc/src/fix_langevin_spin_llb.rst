.. index:: fix langevin/spin/llb

fix langevin/spin/llb command
=============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID langevin/spin/llb keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin/spin/llb = style name of this fix command
* one or more keyword/value pairs; *temp*, *alpha_perp*, and *alpha_par*
  are required, and *seed* is required for thermal noise at a finite
  temperature
* keyword = *temp* or *alpha_perp* or *alpha_par* or *seed* or *noise*

  .. parsed-literal::

       *temp* value = T
         T = temperature of the bath (temperature units)
       *alpha_perp* value = transverse damping of the spin directions (unitless)
       *alpha_par* value = longitudinal damping of the spin magnitudes (unitless)
       *seed* value = random number seed (positive integer)
       *noise* value = *both* or *transverse* or *longitudinal* or *none*
         both = thermal noise on the spin directions and magnitudes
         transverse = thermal noise on the spin directions only
         longitudinal = thermal noise on the spin magnitudes only
         none = damping without thermal noise

Examples
""""""""

.. code-block:: LAMMPS

   fix 2 all langevin/spin/llb temp 300.0 alpha_perp 0.1 alpha_par 0.1 seed 48279
   fix 3 all nve/spin/llb pc lattice moving

Description
"""""""""""

.. versionadded:: TBD

Define the spin bath of the Landau-Lifshitz-Bloch integrator :doc:`fix
nve/spin/llb <fix_nve_spin_llb>`: the temperature, the transverse damping
:math:`\alpha_\perp` of the spin directions, the longitudinal damping
:math:`\alpha_\parallel` of the spin magnitudes, and the random number
seed.  This fix only holds these parameters; fix nve/spin/llb applies the
damping and draws the thermal noise, which has to stay fixed during the
iterations of a timestep.  See that page for the equations of motion.

The two damping constants are parameters of the model and do not depend
on the temperature.  :math:`\alpha_\perp` has the same meaning as the
damping of :doc:`fix langevin/spin <fix_langevin_spin>`, and
:math:`\alpha_\parallel = \alpha_\perp` damps all three components of
the spin vector at the same rate.

The *noise* keyword selects the channels that receive thermal noise.
Without the noise of a channel its damping still acts, so that channel is
relaxed towards lower energy at zero temperature.  The default *both*
samples the canonical distribution.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`; fix nve/spin/llb writes the state of its random numbers.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this
fix.  No global or per-atom quantities are stored by this fix for access
by various :doc:`output commands <Howto_output>`.  This fix is not
invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the SPIN package.  It is only enabled if LAMMPS was
built with that package.  See the :doc:`Build package <Build_package>`
page for more info.

This fix requires :doc:`fix nve/spin/llb <fix_nve_spin_llb>` with the
same group.

Related commands
""""""""""""""""

:doc:`fix nve/spin/llb <fix_nve_spin_llb>`,
:doc:`fix langevin/spin <fix_langevin_spin>`

Default
"""""""

noise = both
