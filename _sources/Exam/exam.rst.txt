Exam
====

Morse potential
---------------

Morse potential sometimes used to describe covalent bonds.
The potential for two connected atoms is:

    .. math::

        V^{Morse}_{ij}=D_e^{ij}\left(1-e^{-a_{ij}(r_{ij}-r_{ij}^0)}\right)^2

Compute the atomic forces due to this potential.

FENE potential
--------------

In some coarse-grained model it important to allow for extensibility of a covalent bond.
These models sometimes use Finitely Extensible Non-linear Elastic potential.
Its mathematical form is:

    .. math::

        V^{FENE}_{ij}=-\frac{k_{ij}}{2}R_{0}^2\log\left[1-(r_{ij}-r_{ij}^0)^2\right]

Compute the atomic forces due to this potential.

Units
-----

Consider a system of particles linked by a set of harmonic springs with spring coefficient :math:`K_s`.

1. Masses of particles are in g/mol, distances are in Angstroms, time is in femtoseconds.

2. Masses of particles are in g/mol, distances are in nm, time is in picoseconds.

Fix the topology file
---------------------

You are given a topology/coordinates pair to run simulations.
However, ``gmx grompp`` exits with error.
Fix the topology file:

1. Use these :download:`topology <FixTopology/1/topol.top>` and :download:`coordinates <FixTopology/1/conf.gro>` files.

1. Use these :download:`topology <FixTopology/2/topol.top>` and :download:`coordinates <FixTopology/2/conf.gro>` files.

End-to-end distance
-------------------

Compute the distribution of the end-to-end distance for hydrocarbon (i.e. distribution of distances between first and last carbon atom throughout the MD simulation run).
Hint: use VMD label functionality.

Surface tension
---------------

Consider a system, in which there are two phase-separation surfaces, parallel to X-Y plane.
In this case, the surface tension can be computed using the following relation [AlejandreJCP95]_:

    .. math::

        \gamma=\frac{L_{X}}{2}\langle P_{ZZ}-\frac{1}{2}\left(P_{XX}+P_{YY}\right)\rangle

Given the energy output of GROMACS, compute the surface tension.

.. [AlejandreJCP95] Alejandre, J., Tildesley, D. J., & Chapela, G. A. (1995). `Molecular dynamics simulation of the orthobaric densities and surface tension of water. <https://aip.scitation.org/doi/pdf/10.1063/1.469505>`_ J. Chem. Phys., 102(11), 4574--4583.