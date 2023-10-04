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

The computational complexity of evaluating forces when each particle is interacting with each other is :math:`\approx N^2`, where :math:`N` is the total number of particles. Various techniques are used to improve the computational performance. One of the most common and most straightforward approach is neighbors or Verlet lists.

The potential vanishes as the distance between particles increase. This allows us to assume that after a certain distance, potential is zero. Hence, the forces for those atoms that are far away can be ignored. We can also assume that during a certain amount of steps atom don't move much. Under this assumption we can construct a list of pairs of particles that are within certain distance :math:`r_{lc}` at the time. Than, assuming that the potential is small enough to be ignored starting from a distance :math:`r_{off}`, we can estimated the frequency of updating the list in number of steps, :math:`s_{lc}`. This number should be small enough so that no particles are traveling father than a difference between :math:`r_{lc}` and :math:`r_{off}`, so that the particles that are outside the list cutoff distance :math:`r_{lc}` at step :math:`s` will not travel into the potential cutoff :math:`r_{off}` at step :math:`s+s_{lc}`. In some MD software the distances traveled by all particles since the last list update is computed. When the maximum of these distances become larger than the difference :math:`r_{lc}-r_{off}`, the list update is triggered. Some software packages require user to define the update frequencies or compute :math:`r_{lc}` and :math:`s_{lc}` from the performance considerations. Note that values of :math:`r_{lc}` and :math:`s_{lc}` do not affect the result, only the performance. Larger the :math:`r_{lc}`, less frequently the list should be updated (larger the :math:`s_{lc}`) and vice-versa.

One important note is that the Lennard-Jones potential asymptotically goes to zero, but never zero. So if we set the potential to zero at a certain length, there will be step (albeit very small step) in the function at this distance. To avoid this step, the function should be smoothed out so it gradually goes to zero before :math:`r_{off}`. Historically two different approaches to do so are used in MD simulations: shift function and switch function.

Shift function
^^^^^^^^^^^^^^

the list is updated after :math:`S_{lc}` number of steps, we can define a constant 


    .. figure:: Figures/NonBonded/VerletListsAndSwitchFunction.png
        :name: Fig:NonBonded:VerletListsAndSwitchFunction

        Verlet or neighboring lists.