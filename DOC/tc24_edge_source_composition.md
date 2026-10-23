# Focused TC24 edge source composition

This branch starts at owner main fa116ee5 plus the already isolated Arnoldi
projection repair9935a10d. It imports the original-equation stopping audit and
native regression tests from owning PR223 commit698d8ac98563675a07d72d60b43162ea4748c563
with full cherry-pick provenance. The CMake conflict retains both existing
binary-split coverage and imported iterator tests.

The Coulomb-only selector is copied from owning PR214 commit
55f183f3b0d3d06f79a5df03f384c16dafde1246, COMMON/collision_model_mod.f90.
Only native/pitch_angle/test_particle/full_fp/monoenergetic_lorentz mode mapping,
its namelist default and config readback are retained. Krook reactions, geometry,
electromagnetic-work changes and field/curl/Phi1 changes are excluded. Existing
native operator matrices and their normalization are unchanged.

`pitch_angle` retains the energy expansion, with isw_lorentz=0, isw_energy=0,
isw_integral=0. `monoenergetic_lorentz` selects the separate single-speed path,
which does not support the intended multispecies experiment. `native` preserves
all legacy switches. A binary becomes a candidate only after native CMake/CTest,
original-equation manufactured oracles and controlled default/model parity.

With an integral part, native JSON covers18 coupled identities for three species,
three sources and two equation systems. Without an integral part it covers six
simultaneous rank-local source solves, each with all three species residuals.
These six vectors represent18 independent block solves, not omitted couplings.
Both paths require every recorded species residual below the same deck tolerance,
all convergence statuses true, exact source/system coverage and no iteration cap.
The physical field contract and transport convergence remain separate admissions.
