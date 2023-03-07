Bonded interactions
===================

    .. math::

        \begin{split}
            V_{bonded} = &\sum_{(i,j)\in bonds}\frac{k^s_{ij}}{2}(r_{ij}-r_{ij}^0)^2 + \\
                         &\sum_{(i,j,k)\in angles}\frac{k^{\theta}_{ijk}}{2}(\theta_{ijk}-\theta_{ijk}^0)^2 + \\
                         &\sum_{(i,j,k,l)\in dihedrals}k^{\phi}_{ijkl}(1-\cos(n\phi_{ijkl} - \phi_{ijkl}^0)) + \\
                         &\sum_{(i,j,k,l)\in impropers}k^{\psi}_{ijkl}(\psi_{ijkl}-\psi_{ijkl}^0)^2
        \end{split}

Harmonic potential
------------------

Harmonic potential usually used to describe covalent bonds.
Its mathematical formula is:

    .. math::

        V_{bond}=\sum_{(i,j)\in bonds}\frac{k^s_{ij}}{2}(r_{ij}-r_{ij}^0)^2=\sum_{(i,j)\in bonds}V_{bond}^{ij}

Here, the sum goes over all bonded pairs :math:`(i,j)`.
The potential parameters :math:`k^s_{ij}` and :math:`r_{ij}^0` are normally based on atom types for atoms :math:`i` and :math:`j`: for each possible combination of types of bonded atoms, there is a set of the two parameters in the forcefield.
The current distance between atoms :math:`r_{ij}` depends on the coordinates of the two connected atoms.
Let us compute the interatomic force for two atoms, connected by a harmonic bond.
The potential function term that describes this binary interaction is:

    .. math::

        V_{bond}^{ij}=\frac{k^s_{ij}}{2}(r_{ij}-r_{ij}^0)^2

Hence, for each bonded interaction, each atom will be experiencing the force:

    .. math::

        \mathbf{f}_{bond_{ij}}^i=\nabla_iV_{bond}^{ij}=\frac{\partial V_{bond}}{\partial r_{ij}}\nabla r_{ij} = -k^s_{ij}(r_{ij}-r_{ij}^0)\frac{\mathbf{r}_{ij}}{r_{ij}}

One can show that :math:`\mathbf{f}_{bond_{ij}}^j = -\mathbf{f}_{bond_{ij}}^i`, which corresponds to the third Newtons law.