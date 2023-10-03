Non-bonded interactions
=======================

    .. math::

        V_{nb} = \sum_{i,j}\left(\frac{q_{i}q_{j}}{4\pi\varepsilon_{0}\varepsilon r_{ij}} + 
                 \varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right]\right)


Van der Waals forces and Lennard-Jones potential
------------------------------------------------

The van der Waals forces describe a combination of strong repulsion of short distances die to Pauli exclusion principle and weak attractive forces on larger distances due to dipole inducing (London dispersion forces). In computational chemistry it is usually described by the Lennard-Jones potential with the following expression:

    .. math::

        V_{nb} = \sum_{i,j}\varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right]

Neighbors or Verlet lists
^^^^^^^^^^^^^^^^^^^^^^^^^

The computational complexity of evaluating forces when each particle is interacting with each other is :math:`\approx N^2`, where :math:`N` is the total number of particles. Various techniques are used to improve the computational performance. One of the most common and most straightforward approach is neighbors oor Verlet lists.

The potential vanishes as the distance between particles increase.


    .. figure:: Figures/NonBonded/VerletListsAndSwitchFunction.png
        :name: Fig:NonBonded:VerletListsAndSwitchFunction

        Verlet or neighboring lists.