.. index:: velocity/tspin

velocity/tspin command
======================

Syntax
""""""

.. code-block:: LAMMPS

   velocity/tspin group-ID create T seed keyword value ...

* group-ID = ID of group of atoms whose spin velocities will be set
* create = style name of this command
* T = desired spin temperature (temperature units, K in metal units)
* seed = random number seed (positive integer)
* zero or more keyword/value pairs may be appended
* keyword = *spinmass* or *mom*

  .. parsed-literal::

       *spinmass* value = tau
         tau = characteristic spin timescale that sets the spin mass (time units)
       *mom* value = *yes* or *no*
         yes = remove net spin momentum for each atom type
         no = retain net spin momentum

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all nve/tspin
   velocity/tspin all create 300.0 12345 spinmass 0.041

Description
"""""""""""

.. versionadded:: TBD

Initialize auxiliary spin velocities for TSPIN at the specified
temperature.  Unlike the ordinary :doc:`velocity <velocity>` command,
this command does not modify atomic velocities.
The temperature is defined by :doc:`compute temp/tspin <compute_temp_tspin>`.

Spin velocities are generated from a Gaussian distribution and rescaled
separately for each atom type to the specified temperature.  With
*mom yes*, the spin-mass-weighted mean velocity of each atom type is
subtracted before rescaling, so that its total auxiliary spin momentum
vanishes,

.. math::

   \sum_{i\in\mathrm{type}}m_{s,i}\vec{v}_i^s = 0.

The *spinmass* keyword sets the spin mass from the characteristic spin
timescale :math:`\tau_s` as described for :doc:`fix nve/tspin
<fix_nve_tspin>`, with the temperature *T* of this command as the
reference temperature; *T* must then be positive.  If omitted,
previously assigned spin masses are used.  When *spinmass* is also
specified in a TSPIN integrator, the integrator overrides these values
at the start of each run, using its own reference temperature.

Restrictions
""""""""""""

This command requires :doc:`atom_style spin <atom_style>` and the SPIN
package; see :doc:`Build package <Build_package>`.  The simulation box,
atoms, and a TSPIN integrator must be defined before using this command.

Related commands
""""""""""""""""

:doc:`velocity <velocity>`,
:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`compute temp/tspin <compute_temp_tspin>`

Default
""""""""

The default is *mom yes*.  Without *spinmass*, existing spin masses
are used.
