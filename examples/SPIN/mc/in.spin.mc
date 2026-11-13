# fcc cobalt: Metropolis Monte Carlo sampling of the spins on a frozen lattice
#
# Each timestep performs one sweep of fix mc/spin at 300 K.  The cone of
# trial directions is adapted during the first run and then kept fixed.
# f_mc[1] is the acceptance rate, f_mc[2] the cone half-angle in degrees.

units           metal
atom_style      spin
atom_modify     map array
boundary        p p p

lattice         fcc 3.54
region          box block 0.0 6.0 0.0 6.0 0.0 6.0
create_box      1 box
create_atoms    1 box

mass            1 58.93
set             group all spin/atom/random 31 1.72

pair_style      hybrid/overlay eam/alloy spin/exchange 4.0
pair_coeff      * * eam/alloy Co_PurjaPun_2012.eam.alloy Co
pair_coeff      * * spin/exchange exchange 4.0 0.0446928 0.003496 1.4885

neighbor        0.1 bin
neigh_modify    every 10 check yes delay 20

fix             prec all precession/spin zeeman 1.0 0.0 0.0 1.0
fix_modify      prec energy yes
fix             mc all mc/spin 300.0 48279 cone 30

compute         mag all spin
variable        magnorm equal c_mag[4]

thermo          50
thermo_style    custom step v_magnorm pe f_mc[1] f_mc[2]
thermo_modify   format float %20.15g

run             200

# production with the adapted cone kept fixed

variable        angle equal f_mc[2]
variable        cone equal ${angle}
unfix           mc
fix             mc all mc/spin 300.0 48280 cone ${cone} adapt no
run             300
