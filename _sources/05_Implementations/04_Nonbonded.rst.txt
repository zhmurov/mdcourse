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

Switching function
^^^^^^^^^^^^^^^^^^

Consider an arbitraty binary potential :math:`V({\bf r})` (:math:`{\bf r}=(\vec{r}_{1},\vec{r}_{2},\ldots,\vec{r}_{N})`) that is subjected to a swithing function :math:`sw({\bf r})``:

    .. math::

        V^{sw}({\bf r})=V({\bf r})sw({\bf r})=\sum_{i=1}^{N}\sum_{j=1}^{N}V_{ij}(r_{ij})sw_{ij}(r_{ij})

According to the chain rule, the atomic force acting on :math:`l`-th particle due to this potential will then be:


    .. math::

        \begin{split}
            \vec{f}_{l}({\bf r})&=-\nabla_{l} V^{sw}=-\nabla_{l}\left(\sum_{i=01}^{N}\sum_{j=1}^{N}V_{ij}(r_{ij})sw_{ij}(r_{ij})\right)=\\
            &=-\sum_{i=1}^{N}\sum_{j=1}^{N}\left(\nabla_{l}V_{ij}(r_{ij})\right)sw_{ij}(r_{ij})+V_{ij}(r_{ij})\left(\nabla_{l}sw_{ij}(r_{ij})\right)
        \end{split}


The actual form of the switching function :math:`sw(r)` is given by :cite:`Charmm83`:

    .. math::
    sw_{ij}(r_{ij})=
    \begin{cases}
        1 & \text{when } r_{ij}\leq r_{on}\\
	    \frac{\left(r_{off}^{2}-r_{ij}^{2}\right)^{2}\left(r_{off}^{2}+2r_{ij}^{2}-3r_{on}^{2}\right)}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}} & \text{when } r_{on}<r_{ij}\leq r_{off}\\
	    0 & \text{when } r_{ij}>r_{off}\\
    \end{cases}

where :math:`r_{on}` and :math:`r_{off}` are switching and cut-off distances respectively.

Since cases when :math:`r_{ij}\leq r_{on}` and when :math:`r_{ij}>r_{off}` are trivial, let us consider conditions when :math:`r_{on}<r_{ij}\leq r_{off}`.


According to the chain rule, one need to evaluate :math:`V_{ij}(r_{ij})`, :math:`\nabla_{l}V_{ij}(r_{ij})`, :math:`sw_{ij}(r_{ij})` and :math:`\nabla_{l}sw_{ij}(r_{ij})`. Since all four quantities can be computed independently, here we will only consider evaluation of :math:`sw_{ij}(r_{ij})` and :math:`\nabla_{l}sw_{ij}(r_{ij})`. Rearranging Eq.~\ref{eq:sw} ($r_{on}<r_{ij}\leq r_{off}$), we get:

    .. math::
    sw_{ij}(r_{ij})= \frac{2}{{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}}\left(r_{off}^{2}-r_{ij}^{2}\right)^{2}\left(\frac{r_{off}^{2}-3r_{on}^{2}}{2}+r_{ij}^{2}\right)

Since :math:`r_{on}` and :math:`r_{off}` are constant, one can use pre-computed values for:

    .. math::
    C^{sw}_{1} = \frac{2}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}} \text{,  }
    C^{sw}_{2} = r_{off}^{2} \text{, and  }
    C^{sw}_{3} = \frac{r_{off}^{2}-3r_{on}^{2}}{2}

When $r_{on}<r_{ij}\leq r_{off}$ becomes:

    .. math::
    sw_{ij}(r_{ij})=C^{sw}_{1}\left(C^{sw}_{2}-r_{ij}^{2}\right)^{2}\left(C^{sw}_{3}+r_{ij}^{2}\right)

Using Eq.~\ref{eq:sw} ($r_{on}<r_{ij}\leq r_{off}$), we also can compute $\nabla_{l}sw_{ij}(r_{ij})$:

    .. math::
    \begin{split}
    &\nabla_{l}sw_{ij}(r_{ij}) = \frac{1}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}\left[2\left(r_{off}^{2}-r_{ij}^{2}\right)\left(-2r_{ij}\right)\left(r_{off}^{2}+2r_{ij}^{2}-3r_{on}^{2}\right)+\left(r_{off}^{2}-r_{ij}^2\right)^{2}4r_{ij}\right]\frac{\vec{r}_{ij}}{r_{ij}} = \\
    &= \frac{1}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}\left[-4r_{off}^{4}-8r_{off}^{2}r_{ij}^{2}+12r_{off}^{2}r_{on}^{2}+4r_{ij}^{2}r_{off}^{2}+8r_{ij}^{4}-12r_{ij}^{2}r_{on}^{2}+4r_{off}^{4}-8r_{off}^{2}r_{ij}^{2}+4r_{ij}^{4}\right]\vec{r}_{ij} = \\
    &= \frac{1}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}\left[12r_{ij}^{4}-12r_{ij}^{2}r_{off}^{2}-12r_{ij}^{2}r_{on}^{2}+12r_{off}^{2}r_{on}^{2}\right]\vec{r}_{ij} = \\
    &= \frac{12}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}\left(r_{ij}^{2}-r_{off}^{2}\right)\left(r_{ij}^2-r_{on}^{2}\right)\vec{r}_{ij}\\

Using pre-computed values for:


    .. math::
    C^{dsw} = \frac{12}{\left(r_{off}^{2}-r_{on}^{2}\right)^{3}}\text{,  } r_{off}^{2} \text{, and  } r_{on}^{2}

Equation above becomes:

    .. math::
    \nabla_{l}sw_{ij}(r_{ij})=C^{dsw}\left(r_{ij}^{2}-r_{off}^{2}\right)\left(r_{ij}^2-r_{on}^{2}\right)\vec{r}_{ij}
