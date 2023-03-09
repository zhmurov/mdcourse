From potential energy to atomic forces
======================================

First, let us derive a couple of general expression that should help us later.
Consider a system of :math:`N` particles.
Assume, that the potential energy of the system depends on the coordinates of the particles :math:`V=V(\{\mathbf{r}_i\})=V(\mathbf{r}_1, \mathbf{r}_2, \ldots \mathbf{r}_N)`.
In this case, the force exerted by a particle :math:`i` is given by:

    .. math::

        \mathbf{f}_i = \nabla_i V = \begin{bmatrix}\frac{\partial}{\partial x_i}\\\frac{\partial}{\partial y_i}\\\frac{\partial}{\partial z_i}\end{bmatrix}V

Potential function usually depends on the relative position of the particles, not their absolute coordinates (e.g. the potential energy of a harmonic spring depends on its extension, not the spacial orientation).
Hence, while differentiating and applying the chain rule, we will frequently encounter expressions like :math:`\nabla_i r_{ij}`, where :math:`r_{ij}` is the distance between two particles.
Expanding this expression, we get:

    .. math::

        \begin{split}
        \nabla_i r_{ij} =& \nabla_i\sqrt{(x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2} = \\
                        =& \frac{1}{2r_{ij}}\nabla_i\left((x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2\right) = \\
                        =& \frac{1}{2r_{ij}}
                            \begin{bmatrix}
                                \frac{\partial(x_i-x_j)^2}{\partial x_i}\\
                                \frac{\partial(y_i-y_j)^2}{\partial y_i}\\
                                \frac{\partial(z_i-z_j)^2}{\partial z_i}
                            \end{bmatrix}
                        = \frac{1}{2r_{ij}}\begin{bmatrix}-2(x_i-x_j)\\-2(y_i-y_j)\\-2(z_i-z_j)\end{bmatrix}
                        = \frac{\mathbf{r}_{ij}}{r_{ij}}
        \end{split}

Note that :math:`\nabla_j r_{ij}=-\nabla_i r_{ij}`, as suggested by Newtons third law.

    .. exercise::

        Compute :math:`\nabla_i\theta_{ijk}`, :math:`\nabla_j\theta_{ijk}` and :math:`\nabla_k\theta_{ijk}`, where :math:`\theta_{ijk}` is the angle between vectors, connecting particles :math:`j-i` and :math:`j-k` (i.e. angle between :math:`\mathbf{r_{ji}}` and :math:`\mathbf{r_{jk}}`).

    .. exercise::

    Compute :math:`\nabla_i\phi_{ijkl}`, :math:`\nabla_j\phi_{ijkl}`, :math:`\nabla_k\phi_{ijkl}` and :math:`\nabla_l\phi_{ijkl}`, where :math:`\phi_{ijkl}` is the dihedral (torsion) angle between :math:`i-j-k` and :math:`j-k-l` planes (i.e. torsion angle for the :math:`\mathbf{r_{jk}}` bond).
