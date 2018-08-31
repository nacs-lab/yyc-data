# For asymmetric trap:

## Parameter used:

* Cs traping frequencies 20, 130, 140
* Na traping frequencies 84% that of Cs

Shifts are calculated based on the first order shift on the ground state. (δ₀).
Quantum number are given in Cs axial, Cs radial X, Cs radial Y,
Na axial, Na radial X, Na radial Y, which is in the same order as the trapping above.

## Parity
The `*_000.csv` files contains the results for even parity along each axis,
which includes the ground state. The `*_010.csv` files contains the results
for odd parity only along the X axis.

## Energies
First line is "State" followed by the δ₀ for each column.
Each following lines contains a state number (space separated Cs and Na motional states),
followed by the energies of the state for the corresponding δ₀.

* [Energy for 000 parity](energy_000.csv)

* [Energy for 010 parity](energy_010.csv)

## Overlaps
Same as the energies, but the energy result is replaced with the wavefunction overlap
with the non-interacting ground state, or in case of the `010` parity, the overlap
with the Na n=1 state in X direction.

* [Overlap for 000 parity](overlap_000.csv)

* [Overlap for 010 parity](overlap_010.csv)
