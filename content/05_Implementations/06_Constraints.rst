Constraints
===========

Equations of motion for the :math:`k`-th particle:

    .. math::

        m_k\frac{d^2\mathbf{r}_k}{dt^2}=-\frac{\partial}{\partial \mathbf{r}_k}V\left(\{\mathbf{r}\}\right)=-\nabla_kV\left(\{\mathbf{r}\}\right)

Same equations with constraints:

    .. math::

        m_k\frac{d^2\mathbf{r}_k}{dt^2}=-\frac{\partial}{\partial \mathbf{r}_k}\left(V\left(\{\mathbf{r}\}\right)+\sum_{c=1}^C\lambda_c\sigma_c\left(\{\mathbf{r}\}\right)\right)=\mathbf{f}_k+\delta\mathbf{f}_k

Where

    .. math::

        \sigma_c\left(\{\mathbf{r}\}\right)=\mathbf{r}_{ml}^2-d_c^2=\mathbf{r}_{ml}^2-d_{ml}^2

Here, :math:`d_c=d_{ml}` is the target distance between :math:`m`-th and :math:`l`-th particles (constrain :math:`c` connects particles :math:`m` and :math:`l`).

The ''constrain'' force :math:`\delta\mathbf{f}_k` is then:

    .. math::

        \delta\mathbf{f}_k=\frac{\partial}{\partial\mathbf{r}_k}\sum_{c=1}^{C}\lambda_c\sigma_c\left(\{\mathbf{r}\}\right)=\sum_{c=1}^{C}\lambda_c\frac{\partial}{\partial\mathbf{r}_k}\sigma_c\left(\{\mathbf{r}\}\right)=\sum_{c=1}^{C}\lambda_c2\delta_{km}\delta_{kl}\mathbf{r}_{ml}

Numerical integration
---------------------

Since

    .. math::

        m_k\frac{d^2\mathbf{r}_k}{dt^2}=\mathbf{f}_k+\delta\mathbf{f}_k

(add numerical scheme here to show how lambda and mu are related)

The numerical integration can be done in two steps: (1) integration as if there are no constraints and (2) adding the displacement due to ''constrain force'' :math:`\delta\mathbf{f}_k`, so that the constraints are satisfied:

    .. math::

        \mathbf{r}_k^{n+1}=\mathbf{r}_k^* + \delta\mathbf{r}_k

Where :math:`\mathbf{r}_k^*` is the displaced coordinates as if there are no constraints and:

    .. math::
        
        \delta\mathbf{r}_k=\sum_{c=1}^C\frac{\mu_c}{m_k}\left(\delta_{km}-\delta_{kl}\right)\mathbf{r}_{ml}

Here, :math:`\mu_c=2\lambda_c\tau^2` is the normalized Lagrange multiplier.

Note that since :math:`\delta\mathbf{\upsilon}=\frac{\delta\mathbf{r}}{\tau}`, the scaled Lagrange multiplier for the velocities is :math:`\mu_c^\upsilon=\frac{\mu_c}{\tau}`.

Solving for :math:`\mu_c`
-------------------------

For each constrain c, after the numerical integration, the distance between constrained atoms should be equal to the target distance of that constraint, :math:`d_c`. Hence, for each constraint we have :math:`|\mathbf{r}_{ij}^{n+1}|=d_{ij}` or :math:`{\mathbf{r}_{ij}^{n+1}}^2=d_{ij}^2`, where :math:`\mathbf{r}^{n+1}_{ij}=\mathbf{r}^{n+1}_{i}-\mathbf{r}^{n+1}_{j}`.

These conditions should be satisfied by finding appropriate values for :math:`\mu_c` (or equally :math:`\lambda_c`). The general form for these quadratic equations is:

    .. math::

        \begin{split}
        \left(\mathbf{r}_{ij}^*+\delta\mathbf{r}_{ij}\right)^2-d_{ij}^2=0 \\
        \left(\mathbf{r}_{ij}^*+\left(\delta\mathbf{r}_{i}-\delta\mathbf{r}_{j}\right)\right)^2-d_{ij}^2=0 \\
        \left(\mathbf{r}_{ij}^*+\left(\sum_{c=1}^C\frac{\mu_c}{m_i}\left(\delta_{im}-\delta_{il}\right)\mathbf{r}_{ml} - \sum_{c=1}^C\frac{\mu_c}{m_j}\left(\delta_{jm}-\delta_{jl}\right)\mathbf{r}_{ml}\right)\right)^2-d_{ij}^2=0
        \end{split}

Finally:

    .. math::
        :label: mu

        \left(\mathbf{r}_{ij}^*+\sum_{c=1}^C\mu_c\left(\frac{\delta_{im}-\delta_{il}}{m_i} - \frac{\delta_{jm}-\delta_{jl}}{m_j}\right)\mathbf{r}_{ml}\right)^2-d_{ij}^2=0

Note that the sum over all constraints :math:`c=\{m,l\}` should include the pair :math:`\{i,j\}`, in which case only two of the kronecker's deltas will survive. For a coupled constraint, when either :math:`i` or :math:`j` is equal to one index in the :math:`\{m,l\}` pair, only one delta will not be zero. Since :math:`i\ne j` and :math:`m\ne l`, there is never going to be the case when more than two deltas are not zero.

Notations
---------

Before applying the above formula to specific cases, let us simplify the notations. First, let us us the :math:`r_i^0` for the coordinates of the atom before the integration and :math:`r_i` as the coordinates for the same particle after the integration and constraining is done. Second, let us assume that the constrain connecting atoms :math:`m` and :math:`l` has an index :math:`c`, and we will use the constraint index :math:`c` and indices of the connected particles :math:`m` and :math:`l` interchangeably (we already did this for :math:`d_c=d_{ml}` above). We will be assigning these indexes for each case before deriving the formulas. This should allow us to have less noise in the formulas and the generalization for the :math:`n`-th step and for arbitrary indices of particles should be straightforward.

One constraint
--------------

The simplest case by far is when we have a single constraint, which is not coupled to any other. An example should be the constrained C-H bond in the protein backbone in the case when only bonds with hydrogens are constrained (with the only exception of glycine, where there is another hydrogen instead of a side-chain). Indeed, this case is very simple: we have the direction of the bond in which atoms should be moved. And the sum of the displacements for two moved atoms should be so that the final distance is equal to the target. The displacements are reversely proportional to masses of atoms. Let us derive this formally, using the :eq:`mu` above.

We have one constraint :math:`c=1={0,1}` that connects two atoms :math:`i=0` and :math:`j=1` with masses :math:`m_0` and :math:`m_1`. We will also denote all the connecting vectors with constrain index, i.e. :math:`\mathbf{r}_1=\mathbf{r}_{01}=\mathbf{r}_0-\mathbf{r}_1`, etc. so there is less noise in formulas. Note that we could remove the index altogether since there is only one constrained bond, but we will keep the index here for consistency with the rest of the cases. The sole equation derived from :eq:`mu` is:

    .. math::

        \left(\mathbf{r}_{1}^*+\mu_1\left(\frac{1}{m_0}+\frac{1}{m_1}\right)\mathbf{r}_1^0\right)^2-d_1^2=0

With reduced mass :math:`M_1=\frac{m_0m_1}{m_0+m_1}`:

    .. math::

        \left(\mathbf{r}_{1}^*+\frac{\mu_1}{M_1}\mathbf{r}_1^0\right)^2-d_1^2=0

Expanding:

    .. math::

        \mu_1^2\frac{{\mathbf{r}_1^0}^2}{M_1^2}+2\mu_1\frac{\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}+\left({\mathbf{r}_1^*}^2-d_1^2\right)=0

This is a single quadratic equation for :math:`\mu_1`, with the following solution:

    .. math::

        \begin{split}
        D_2=\frac{1}{M_1^2}\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right) \\
        {\mu_1}_{1,2}=\frac{\pm\sqrt{D_2}-\frac{\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}}{\frac{{\mathbf{r}_1^0}^2}{M_1^2}}
        \end{split}

Plugging the :math:`D_2` into the solution gives:

    .. math::

        {\mu_1}_{1,2}=M_1\frac{\pm\sqrt{\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right)}-\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{{\mathbf{r}_1^0}^2}

As expected for quadratic equation, we have two solutions. However only one is valid for our problem: the one that gives smaller absolute value of :math:`mu`. This is because we assume that the bond length after the numerical integration without constraints is close to the target length, i.e. the displacements of the atoms is small. Indeed, we can move the atoms so they satisfy the constraint two ways, one of which employs flipping the bond. The solution we are interested in is the one with the plus sign before the square root. So the final solution for :math:`\mu_1` is:

    .. math::

        {\mu_1}_{1,2}=M_1\frac{\sqrt{\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right)}-\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{{\mathbf{r}_1^0}^2}

Two coupled constraints
-----------------------

