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

        \delta\mathbf{f}_k=-\frac{\partial}{\partial\mathbf{r}_k}\sum_{c=1}^{C}\lambda_c\sigma_c\left(\{\mathbf{r}\}\right)=-\sum_{c=1}^{C}\lambda_c\frac{\partial}{\partial\mathbf{r}_k}\sigma_c\left(\{\mathbf{r}\}\right)=-2\sum_{c=1}^{C}\lambda_c\left(\delta_{km}-\delta_{kl}\right)\mathbf{r}_{ml}

Numerical integration
---------------------

Since

    .. math::

        m_k\frac{d^2\mathbf{r}_k}{dt^2}=\mathbf{f}_k+\delta\mathbf{f}_k

For example, with timestep :math:`\tau`, the leapfrog scheme is

    .. math::

        \begin{split}
        \mathbf{v}_k^{n+1/2}&=\mathbf{v}_k^{n-1/2}+\frac{\tau}{m_k}\left(\mathbf{f}_k^n+\delta\mathbf{f}_k^n\right), \\
        \mathbf{r}_k^{n+1}&=\mathbf{r}_k^n+\tau\mathbf{v}_k^{n+1/2}.
        \end{split}

The constraint force is evaluated along the bond directions at the old positions :math:`\mathbf{r}^n`; its multipliers are chosen so that the new positions satisfy the constraints.

The numerical integration can be done in two steps: (1) integration as if there are no constraints and (2) adding the displacement due to ''constrain force'' :math:`\delta\mathbf{f}_k`, so that the constraints are satisfied:

    .. math::

        \mathbf{r}_k^{n+1}=\mathbf{r}_k^* + \delta\mathbf{r}_k

Where :math:`\mathbf{r}_k^*` is the displaced coordinates as if there are no constraints and:

    .. math::
        
        \delta\mathbf{r}_k=\frac{\tau^2}{m_k}\delta\mathbf{f}_k^n=\sum_{c=1}^C\frac{\mu_c}{m_k}\left(\delta_{km}-\delta_{kl}\right)\mathbf{r}_{ml}^n

Here, :math:`\mu_c=-2\lambda_c^n\tau^2` is the displacement multiplier for this leapfrog scheme and the sign convention used above. The relation between force and displacement multipliers depends on the integration scheme.

For the leapfrog half-step velocities, :math:`\delta\mathbf{v}_k^{n+1/2}=\delta\mathbf{r}_k/\tau`, so the corresponding velocity multiplier is :math:`\mu_c^v=\mu_c/\tau`. This update does not by itself impose the end-of-step velocity constraints in a velocity-Verlet scheme; those require a separate velocity projection, as in RATTLE.

Solving for :math:`\mu_c`
-------------------------

For each constrain c, after the numerical integration, the distance between constrained atoms should be equal to the target distance of that constraint, :math:`d_c`. Hence, for each constraint we have :math:`|\mathbf{r}_{ij}^{n+1}|=d_{ij}` or :math:`{\mathbf{r}_{ij}^{n+1}}^2=d_{ij}^2`, where :math:`\mathbf{r}^{n+1}_{ij}=\mathbf{r}^{n+1}_{i}-\mathbf{r}^{n+1}_{j}`.

These conditions should be satisfied by finding appropriate values for :math:`\mu_c` (or equally :math:`\lambda_c`). The general form for these quadratic equations is:

    .. math::

        \begin{split}
        \left(\mathbf{r}_{ij}^*+\delta\mathbf{r}_{ij}\right)^2-d_{ij}^2=0 \\
        \left(\mathbf{r}_{ij}^*+\left(\delta\mathbf{r}_{i}-\delta\mathbf{r}_{j}\right)\right)^2-d_{ij}^2=0 \\
        \left(\mathbf{r}_{ij}^*+\left(\sum_{c=1}^C\frac{\mu_c}{m_i}\left(\delta_{im}-\delta_{il}\right)\mathbf{r}_{ml}^n - \sum_{c=1}^C\frac{\mu_c}{m_j}\left(\delta_{jm}-\delta_{jl}\right)\mathbf{r}_{ml}^n\right)\right)^2-d_{ij}^2=0
        \end{split}

Finally:

    .. math::
        :label: mu

        \left(\mathbf{r}_{ij}^*+\sum_{c=1}^C\mu_c\left(\frac{\delta_{im}-\delta_{il}}{m_i} - \frac{\delta_{jm}-\delta_{jl}}{m_j}\right)\mathbf{r}_{ml}^n\right)^2-d_{ij}^2=0

Note that the sum over all constraints :math:`c=\{m,l\}` should include the pair :math:`\{i,j\}`, in which case only two of the kronecker's deltas will survive. For a coupled constraint, when either :math:`i` or :math:`j` is equal to one index in the :math:`\{m,l\}` pair, only one delta will not be zero. Since :math:`i\ne j` and :math:`m\ne l`, there is never going to be the case when more than two deltas are not zero.

Notations
---------

Before applying the above formula to specific cases, let us simplify the notation. We write :math:`\mathbf{x}_i^0`, :math:`\mathbf{x}_i^*`, and :math:`\mathbf{x}_i` for the old, trial, and corrected positions of atom :math:`i`. For an oriented constraint :math:`c=(i,j)`, we use :math:`\mathbf{r}_c^0=\mathbf{x}_i^0-\mathbf{x}_j^0`, :math:`\mathbf{r}_c^*=\mathbf{x}_i^*-\mathbf{x}_j^*`, and :math:`\mathbf{r}_c=\mathbf{x}_i-\mathbf{x}_j` for its bond vectors. Thus atom indices on :math:`\mathbf{x}` and constraint indices on :math:`\mathbf{r}` are distinct. The superscript :math:`0` denotes the old integration step, not a Newton iteration. Constraint orientations and indices are specified separately for each case. All bond vectors must use consistent periodic images.

One constraint
--------------

The simplest case by far is when we have a single constraint, which is not coupled to any other. An example should be the constrained C-H bond in the protein backbone in the case when only bonds with hydrogens are constrained (with the only exception of glycine, where there is another hydrogen instead of a side-chain). Indeed, this case is very simple: we have the direction of the bond in which atoms should be moved. And the sum of the displacements for two moved atoms should be so that the final distance is equal to the target. The displacements are reversely proportional to masses of atoms. Let us derive this formally, using the :eq:`mu` above.

    .. figure:: Figures/Constraints/CH1.png
        :name: Fig:Constraints:CH1
        :width: 200

        A schematic representation of a constrained bond. Index on top of the bond represents the constraint index :math:`c`. Arrow indicates the direction of the vector :math:`\mathbf{r}_1`. Masses subscripts correspond to aton indices.

We have one constraint :math:`c=1=\{0,1\}` that connects two atoms :math:`i=0` and :math:`j=1` with masses :math:`m_0` and :math:`m_1`. We will also denote all the connecting vectors with constrain index, i.e. :math:`\mathbf{r}_1=\mathbf{x}_0-\mathbf{x}_1`, etc. so there is less noise in formulas. Note that we could remove the index altogether since there is only one constrained bond, but we will keep the index here for consistency with the rest of the cases. The sole equation derived from :eq:`mu` is:

    .. math::

        \left(\mathbf{r}_{1}^*+\mu_1\left(\frac{1}{m_0}+\frac{1}{m_1}\right)\mathbf{r}_1^0\right)^2-d_1^2=0

With reduced mass :math:`M_1=\frac{m_0m_1}{m_0+m_1}`:

    .. math::

        \left(\mathbf{r}_{1}^*+\frac{\mu_1}{M_1}\mathbf{r}_1^0\right)^2-d_1^2=0

Expanding:

    .. math::
        :label: EquationOneUncoupled

        \mu_1^2\frac{{\mathbf{r}_1^0}^2}{M_1^2}+2\mu_1\frac{\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}+\left({\mathbf{r}_1^*}^2-d_1^2\right)=0

This is a single quadratic equation for :math:`\mu_1`, with the following solution:

    .. math::

        \begin{split}
        D_2=\frac{1}{M_1^2}\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right) \\
        {\mu_1}^{1,2}=\frac{\pm\sqrt{D_2}-\frac{\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}}{\frac{{\mathbf{r}_1^0}^2}{M_1^2}}
        \end{split}

Plugging the :math:`D_2` into the solution gives:

    .. math::

        {\mu_1}^{1,2}=M_1\frac{\pm\sqrt{\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right)}-\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{{\mathbf{r}_1^0}^2}

For a positive discriminant, there are two real solutions. We assume that the trial bond vector is close to the old one, not merely that its length is close to the target. In this regime, :math:`\mathbf{r}_1^*\cdot\mathbf{r}_1^0>0`, and the solution of smaller absolute magnitude gives the small correction connected to the previous configuration. The other root reverses the component of the bond along its old direction. A zero discriminant gives a repeated root; a negative discriminant means that no real correction along the old bond direction can satisfy this isolated constraint.

    .. figure:: Figures/Constraints/CH1_mu12.png
        :name: Fig:Constraints:CH1_mu12
        :width: 400

        Two roots give two different directions of the bond after constraining. Given that the timestep is small, the displacement should be small as well. Hence the correct solution for :math:`\mu_1` is the smallest by absolute value.

The solution with the plus sign before the square root is the small root under the positive-dot-product assumption above. The resulting expression for :math:`\mu_1` is:

    .. math::
        :label: UncoupledMu

        \mu_1=M_1\frac{\sqrt{\left(\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)^2-{\mathbf{r}_1^0}^2\left({\mathbf{r}_1^*}^2-d_1^2\right)\right)}-\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{{\mathbf{r}_1^0}^2}

For numerical evaluation near an already satisfied constraint, the subtraction in :eq:`UncoupledMu` can lose significant digits. Define :math:`a={\mathbf{r}_1^0}^2`, :math:`b=\mathbf{r}_1^*\cdot\mathbf{r}_1^0`, and :math:`c={\mathbf{r}_1^*}^2-d_1^2`. An algebraically equivalent, cancellation-resistant expression for :math:`b>0` is

    .. math::
        :label: UncoupledMuStable

        \mu_1=-M_1\frac{c}{b+\sqrt{b^2-ac}}.

The discriminant must be checked before taking the square root. If :math:`b\leq0`, the small-step assumption has failed and the plus-root choice is no longer justified by this argument.

Though this equation is quite simple, we need to prepare ourselves to the following cases. So let us introduce the notations for all the coefficients in the :eq:`EquationOneUncoupled`, re-writing it as:

    .. math::

        k_1^{11}\mu_1^2+k_1^1\mu_1+k_1^0=0

Here, the indices for parameters :math:`k_x^y` are chosen the following way. The bottom index is the number of equation (we have one so far). The top indexes indicate what combination of variables :math:`mu` this coefficient is multiplied by: :math:`k_1^{11}` is for :math:`\mu_1\mu_1=\mu^2`, :math:`k_1^1` is for :math:`\mu_1` and :math:`0` is for free coefficient. Later we will encounter other variables (:math:`\mu_2`, :math:`\mu_3`, etc.), for which we will have :math:`k_1^{22}`, :math:`k_1^{33}` as well as some mixed products, e.g. :math:`\mu_1\mu_2` with :math:`k_1^{12}`, etc. As follows from :eq:`EquationOneUncoupled`, in the case of one uncoupled constraint, there are three coefficient:

    .. math::

        k_1^{11}=\frac{{\mathbf{r}_1^0}^2}{M_1^2}\mathrm{,~~}
        k_1^1=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}\mathrm{,~~}
        k_1^0={\mathbf{r}_1^*}^2-d_1^2

In these notations, the solution for :math:`\mu_1` is:

    .. math::
        :label: UncoupledMuThroughK

        \mu_1=\frac{\sqrt{{k_1^1}^2-4k_1^{11}k_1^0}-k_1^1}{2k_1^{11}}

Equations :eq:`UncoupledMuThroughK` and :eq:`UncoupledMu` are algebraically equivalent; :eq:`UncoupledMuStable` avoids cancellation when the constraint residual is small.

Two coupled constraints
-----------------------

    .. figure:: Figures/Constraints/CH2.png
        :name: Fig:Constraints:CH2
        :width: 300

        A schematic representation of two coupled constrained bonds.

Two coupled constraints will give us a much more complicated case of two coupled quadratic equations. Assuming that the central atom has index :math:`0`, with atoms :math:`1` and :math:`2` connected to it with two constrained bonds :math:`c=1=\{0,1\}` and :math:`c=2=\{0,2\}`, and using the reduced masses :math:`M_1=\frac{m_0m_1}{m_0+m_1}` and :math:`M_2=\frac{m_0m_2}{m_0+m_2}`, from :eq:`mu` we get:

    .. math::

        \begin{cases}
        \left(\mathbf{r}_{1}^*+\frac{\mu_1}{M_1}\mathbf{r}_1^0+\frac{\mu_2}{m_0}\mathbf{r}_2^0\right)^2-d_1^2=0 \\
        \left(\mathbf{r}_{2}^*+\frac{\mu_1}{m_0}\mathbf{r}_1^0+\frac{\mu_2}{M_2}\mathbf{r}_2^0\right)^2-d_2^2=0
        \end{cases}

    .. math::
        :label: SystemTwoCoupled

        \begin{cases}
        k_1^{11}\mu_1^2+k_1^{22}\mu_2^2+k_1^{12}\mu_1\mu_2+k_1^1\mu_1+k_1^2\mu_2+k_1^0=0 \\
        k_2^{11}\mu_1^2+k_2^{22}\mu_2^2+k_2^{12}\mu_1\mu_2+k_2^1\mu_1+k_2^2\mu_2+k_2^0=0
        \end{cases}

Where

    .. math::
        :label: SystemTwoCoupledKs

        \begin{split}
            k_1^{11}=\frac{{\mathbf{r}_1^0}^2}{M_1^2}\mathrm{,~~}
            k_1^{22}=\frac{{\mathbf{r}_2^0}^2}{m_0^2}\mathrm{,~~}
            k_1^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_1m_0}\mathrm{,} \\
            k_1^1=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}\mathrm{,~~}
            k_1^2=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_2^0\right)}{m_0}\mathrm{,~~}
            k_1^0={\mathbf{r}_1^*}^2-d_1^2\mathrm{,~~} \\ \\
            k_2^{11}=\frac{{\mathbf{r}_1^0}^2}{m_0^2}\mathrm{,~~}
            k_2^{22}=\frac{{\mathbf{r}_2^0}^2}{M_2^2}\mathrm{,~~}
            k_2^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_2m_0}\mathrm{,} \\
            k_2^1=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_1^0\right)}{m_0}\mathrm{,~~}
            k_2^2=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_2^0\right)}{M_2}\mathrm{,~~}
            k_2^0={\mathbf{r}_2^*}^2-d_2^2\mathrm{.}
        \end{split}

:eq:`SystemTwoCoupled` are two coupled quadratic equations of a general form. Unfortunately, there is no simple analytical solution to it. Hence the numerical method should be used. Before describing the methods used currently in most MD software packages, let us try to introduce a method for the specific cases we are considering here.

If we introduce two dimensional vector-function :math:`F(\mu_1, \mu_2)` in :eq:`SystemTwoCoupled` it becomes:

    .. math::

        \mathbf{F}(\mu_1, \mu_2) = \mathbf{0}

This is a system of non-linear equations, which are notoriously hard to solve if we (1) need to find all the solutions and (2) don't have good initial guess for these solutions. However, this is not the case here. Firstly, we assume that atoms did satisfy the constraints on the previous step (so we seek the solution near the zero multiplier vector). Secondly, we only need to find this solution (see how we dropped one of the roots when solving for uncoupled constraint). Hence we can apply the Newtons method (also known as Newton-Raphson method). The general form of in in multi-dimensional case is:

    .. math::
        :label: SystemTwoCoupledNewton

        \mathbf{\mu}^{n+1}=\mathbf{\mu}^n-J_\mathbf{F}\left(\mathbf{\mu}^n\right)^{-1}\mathbf{F}\left(\mathbf{\mu}^n\right)

Here, :math:`\mathbf{\mu}^n=(\mu_1^n,\mu_2^n)^T` is a vector with the solutions on :math:`n`-th iteration, :math:`J_\mathbf{F}\left(\mathbf{\mu}\right)=J_\mathbf{F}\left(\mu_1, \mu_2\right)` is the Jacobian matrix for the system. Using :eq:`SystemTwoCoupled`, we can compute :math:`J_\mathbf{F}\left(\mu_1, \mu_2\right)`:

    .. math::

        J_\mathbf{F}=
        \begin{pmatrix}
        2k_1^{11}\mu_1+k_1^{12}\mu_2+k_1^1 & k_1^{12}\mu_1+2k_1^{22}\mu_2+k_1^2 \\
        2k_2^{11}\mu_1+k_2^{12}\mu_2+k_2^1 & k_2^{12}\mu_1+2k_2^{22}\mu_2+k_2^2
        \end{pmatrix}

With everything in the :eq:`SystemTwoCoupledNewton` derived, we need good initial approximation for :math:`\mu_1` and :math:`mu_2`. One can use zeroes, assuming that atoms did not moved much during the integration time step atoms did not moved far away from their previous positions (which did satisfy the constraints). An alternative initial approximation is obtained by assuming that the constraints are not coupled, i.e. by using :eq:`UncoupledMu` or :eq:`UncoupledMuThroughK` for two constraints we have in our case:

    .. math::
        :label: CoupledMyAsUncoupledThroughK

        \begin{split}
        \mu_1^0=\frac{\sqrt{{k_1^1}^2-4k_1^{11}k_1^0}-k_1^1}{2k_1^{11}} \\
        \mu_2^0=\frac{\sqrt{{k_2^2}^2-4k_2^{22}k_2^0}-k_2^2}{2k_2^{22}}
        \end{split}

The numerical procedure is as follows. (1) Compute the coefficients in :eq:`SystemTwoCoupledKs`. (2) Use the independent-constraint roots in :eq:`CoupledMyAsUncoupledThroughK` as initial estimates, evaluating them with the stable form :eq:`UncoupledMuStable` when applicable. (3) Evaluate :math:`\mathbf{F}` and :math:`J_\mathbf{F}`. (4) Solve :math:`J_\mathbf{F}\Delta\boldsymbol{\mu}=-\mathbf{F}` and update :math:`\boldsymbol{\mu}\leftarrow\boldsymbol{\mu}+\Delta\boldsymbol{\mu}`. Solving this small linear system avoids explicitly forming the inverse. (5) Repeat until :math:`\max_c |F_c|/d_c^2` is below the chosen tolerance, then apply the atomic displacements and verify the resulting bond lengths. This is a relative squared-distance tolerance; for small errors it is approximately twice the relative length error. Newton convergence is local, so invalid initial roots, singular Jacobians, and failure to converge within an iteration limit must be handled explicitly.


Three constraints coupled through the central atom
--------------------------------------------------

    .. figure:: Figures/Constraints/CH3.png
        :name: Fig:Constraints:CH3
        :width: 300

        A schematic representation of three atoms constrained to one central atom.

    .. math::

        \begin{cases}
        \left(\mathbf{r}_{1}^*+\frac{\mu_1}{M_1}\mathbf{r}_1^0+\frac{\mu_2}{m_0}\mathbf{r}_2^0+\frac{\mu_3}{m_0}\mathbf{r}_3^0\right)^2-d_1^2=0 \\
        \left(\mathbf{r}_{2}^*+\frac{\mu_1}{m_0}\mathbf{r}_1^0+\frac{\mu_2}{M_2}\mathbf{r}_2^0+\frac{\mu_3}{m_0}\mathbf{r}_3^0\right)^2-d_2^2=0 \\
        \left(\mathbf{r}_{3}^*+\frac{\mu_1}{m_0}\mathbf{r}_1^0+\frac{\mu_2}{m_0}\mathbf{r}_2^0+\frac{\mu_3}{M_3}\mathbf{r}_3^0\right)^2-d_3^2=0 \\
        \end{cases}

    .. math::
        :label: SystemThreeCoupled

        \begin{cases}
        k_1^{11}\mu_1^2+k_1^{22}\mu_2^2+k_1^{33}\mu_3^2+k_1^{12}\mu_1\mu_2+k_1^{13}\mu_1\mu_3+k_1^{23}\mu_2\mu_3+k_1^1\mu_1+k_1^2\mu_2+k_1^3\mu_3+k_1^0=0 \\
        k_2^{11}\mu_1^2+k_2^{22}\mu_2^2+k_2^{33}\mu_3^2+k_2^{12}\mu_1\mu_2+k_2^{13}\mu_1\mu_3+k_2^{23}\mu_2\mu_3+k_2^1\mu_1+k_2^2\mu_2+k_2^3\mu_3+k_2^0=0 \\
        k_3^{11}\mu_1^2+k_3^{22}\mu_2^2+k_3^{33}\mu_3^2+k_3^{12}\mu_1\mu_2+k_3^{13}\mu_1\mu_3+k_3^{23}\mu_2\mu_3+k_3^1\mu_1+k_3^2\mu_2+k_3^3\mu_3+k_3^0=0
        \end{cases}

Where

    .. math::
        :label: SystemThreeCoupledKs

        \begin{split}
            k_1^{11}=\frac{{\mathbf{r}_1^0}^2}{M_1^2}\mathrm{,~~}
            k_1^{22}=\frac{{\mathbf{r}_2^0}^2}{m_0^2}\mathrm{,~~}
            k_1^{33}=\frac{{\mathbf{r}_3^0}^2}{m_0^2}\mathrm{,} \\
            k_1^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_1m_0}\mathrm{,~~}
            k_1^{13}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{M_1m_0}\mathrm{,~~}
            k_1^{23}=\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{m_0^2}\mathrm{,} \\
            k_1^1=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}\mathrm{,~~}
            k_1^2=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_2^0\right)}{m_0}\mathrm{,~~}
            k_1^3=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_3^0\right)}{m_0}\mathrm{,~~}
            k_1^0={\mathbf{r}_1^*}^2-d_1^2\mathrm{,~~} \\ \\
            k_2^{11}=\frac{{\mathbf{r}_1^0}^2}{m_0^2}\mathrm{,~~}
            k_2^{22}=\frac{{\mathbf{r}_2^0}^2}{M_2^2}\mathrm{,~~}
            k_2^{33}=\frac{{\mathbf{r}_3^0}^2}{m_0^2}\mathrm{,} \\
            k_2^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_2m_0}\mathrm{,~~}
            k_2^{13}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{m_0^2}\mathrm{,~~}
            k_2^{23}=\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{M_2m_0}\mathrm{,} \\
            k_2^1=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_1^0\right)}{m_0}\mathrm{,~~}
            k_2^2=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_2^0\right)}{M_2}\mathrm{,~~}
            k_2^3=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_3^0\right)}{m_0}\mathrm{,~~}
            k_2^0={\mathbf{r}_2^*}^2-d_2^2\mathrm{,~~} \\ \\
            k_3^{11}=\frac{{\mathbf{r}_1^0}^2}{m_0^2}\mathrm{,~~}
            k_3^{22}=\frac{{\mathbf{r}_2^0}^2}{m_0^2}\mathrm{,~~}
            k_3^{33}=\frac{{\mathbf{r}_3^0}^2}{M_3^2}\mathrm{,} \\
            k_3^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{m_0^2}\mathrm{,~~}
            k_3^{13}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{M_3m_0}\mathrm{,~~}
            k_3^{23}=\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{M_3m_0}\mathrm{,} \\
            k_3^1=\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_1^0\right)}{m_0}\mathrm{,~~}
            k_3^2=\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_2^0\right)}{m_0}\mathrm{,~~}
            k_3^3=\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_3^0\right)}{M_3}\mathrm{,~~}
            k_3^0={\mathbf{r}_3^*}^2-d_3^2\mathrm{,~~} \\ \\
        \end{split}


    .. math::

        \begin{split}
        &J_\mathbf{F}=\\
        &\begin{pmatrix}
        2k_1^{11}\mu_1+k_1^{12}\mu_2+k_1^{13}\mu_3+k_1^1 & k_1^{12}\mu_1+2k_1^{22}\mu_2+k_1^{23}\mu_3+k_1^2 & k_1^{13}\mu_1+k_1^{23}\mu_2+2k_1^{33}\mu_3+k_1^3 \\
        2k_2^{11}\mu_1+k_2^{12}\mu_2+k_2^{13}\mu_3+k_2^1 & k_2^{12}\mu_1+2k_2^{22}\mu_2+k_2^{23}\mu_3+k_2^2 & k_2^{13}\mu_1+k_2^{23}\mu_2+2k_2^{33}\mu_3+k_2^3 \\
        2k_3^{11}\mu_1+k_3^{12}\mu_2+k_3^{13}\mu_3+k_3^1 & k_3^{12}\mu_1+2k_3^{22}\mu_2+k_3^{23}\mu_3+k_3^2 & k_3^{13}\mu_1+k_3^{23}\mu_2+2k_3^{33}\mu_3+k_3^3 \\
        \end{pmatrix}
        \end{split}

Triangle of constraints
-----------------------

    .. figure:: Figures/Constraints/H2O.png
        :name: Fig:Constraints:H2O
        :width: 300

        Three constraints can form a rigid triangle, e.g. in water molecule.

Label the atoms :math:`1,2,3`, as in the figure, and orient the bond vectors along its arrows:

    .. math::

        \mathbf{r}_1=\mathbf{x}_2-\mathbf{x}_1,\qquad
        \mathbf{r}_2=\mathbf{x}_3-\mathbf{x}_2,\qquad
        \mathbf{r}_3=\mathbf{x}_1-\mathbf{x}_3.

Use the same orientations for the old and trial vectors. The reduced masses are :math:`M_1=m_1m_2/(m_1+m_2)`, :math:`M_2=m_2m_3/(m_2+m_3)`, and :math:`M_3=m_3m_1/(m_3+m_1)`. Then :eq:`mu` gives

    .. math::

        \begin{cases}
        \left(\mathbf{r}_{1}^*+\frac{\mu_1}{M_1}\mathbf{r}_1^0-\frac{\mu_2}{m_2}\mathbf{r}_2^0-\frac{\mu_3}{m_1}\mathbf{r}_3^0\right)^2-d_1^2=0 \\
        \left(\mathbf{r}_{2}^*-\frac{\mu_1}{m_2}\mathbf{r}_1^0+\frac{\mu_2}{M_2}\mathbf{r}_2^0-\frac{\mu_3}{m_3}\mathbf{r}_3^0\right)^2-d_2^2=0 \\
        \left(\mathbf{r}_{3}^*-\frac{\mu_1}{m_1}\mathbf{r}_1^0-\frac{\mu_2}{m_3}\mathbf{r}_2^0+\frac{\mu_3}{M_3}\mathbf{r}_3^0\right)^2-d_3^2=0 \\
        \end{cases}


    .. math::
        :label: SystemTriangleKs

        \begin{split}
            k_1^{11}=\frac{{\mathbf{r}_1^0}^2}{M_1^2}\mathrm{,~~}
            k_1^{22}=\frac{{\mathbf{r}_2^0}^2}{m_2^2}\mathrm{,~~}
            k_1^{33}=\frac{{\mathbf{r}_3^0}^2}{m_1^2}\mathrm{,} \\
            k_1^{12}=-\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_1m_2}\mathrm{,~~}
            k_1^{13}=-\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{M_1m_1}\mathrm{,~~}
            k_1^{23}=\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{m_1m_2}\mathrm{,} \\
            k_1^1=\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_1^0\right)}{M_1}\mathrm{,~~}
            k_1^2=-\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_2^0\right)}{m_2}\mathrm{,~~}
            k_1^3=-\frac{2\left(\mathbf{r}_1^*\cdot\mathbf{r}_3^0\right)}{m_1}\mathrm{,~~}
            k_1^0={\mathbf{r}_1^*}^2-d_1^2\mathrm{,~~} \\ \\
            k_2^{11}=\frac{{\mathbf{r}_1^0}^2}{m_2^2}\mathrm{,~~}
            k_2^{22}=\frac{{\mathbf{r}_2^0}^2}{M_2^2}\mathrm{,~~}
            k_2^{33}=\frac{{\mathbf{r}_3^0}^2}{m_3^2}\mathrm{,} \\
            k_2^{12}=-\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{M_2m_2}\mathrm{,~~}
            k_2^{13}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{m_2m_3}\mathrm{,~~}
            k_2^{23}=-\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{M_2m_3}\mathrm{,} \\
            k_2^1=-\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_1^0\right)}{m_2}\mathrm{,~~}
            k_2^2=\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_2^0\right)}{M_2}\mathrm{,~~}
            k_2^3=-\frac{2\left(\mathbf{r}_2^*\cdot\mathbf{r}_3^0\right)}{m_3}\mathrm{,~~}
            k_2^0={\mathbf{r}_2^*}^2-d_2^2\mathrm{,~~} \\ \\
            k_3^{11}=\frac{{\mathbf{r}_1^0}^2}{m_1^2}\mathrm{,~~}
            k_3^{22}=\frac{{\mathbf{r}_2^0}^2}{m_3^2}\mathrm{,~~}
            k_3^{33}=\frac{{\mathbf{r}_3^0}^2}{M_3^2}\mathrm{,} \\
            k_3^{12}=\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_2^0\right)}{m_1m_3}\mathrm{,~~}
            k_3^{13}=-\frac{2\left(\mathbf{r}_1^0\cdot\mathbf{r}_3^0\right)}{M_3m_1}\mathrm{,~~}
            k_3^{23}=-\frac{2\left(\mathbf{r}_2^0\cdot\mathbf{r}_3^0\right)}{M_3m_3}\mathrm{,} \\
            k_3^1=-\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_1^0\right)}{m_1}\mathrm{,~~}
            k_3^2=-\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_2^0\right)}{m_3}\mathrm{,~~}
            k_3^3=\frac{2\left(\mathbf{r}_3^*\cdot\mathbf{r}_3^0\right)}{M_3}\mathrm{,~~}
            k_3^0={\mathbf{r}_3^*}^2-d_3^2\mathrm{,~~} \\ \\
        \end{split}

The cyclic bond vectors satisfy :math:`\mathbf{r}_1^0+\mathbf{r}_2^0+\mathbf{r}_3^0=0`. Consequently,

    .. math::

        \begin{split}
        \mathbf{r}_1^0\cdot\mathbf{r}_2^0&=\frac{{\mathbf{r}_3^0}^2-{\mathbf{r}_1^0}^2-{\mathbf{r}_2^0}^2}{2}, \\
        \mathbf{r}_2^0\cdot\mathbf{r}_3^0&=\frac{{\mathbf{r}_1^0}^2-{\mathbf{r}_2^0}^2-{\mathbf{r}_3^0}^2}{2}, \\
        \mathbf{r}_3^0\cdot\mathbf{r}_1^0&=\frac{{\mathbf{r}_2^0}^2-{\mathbf{r}_3^0}^2-{\mathbf{r}_1^0}^2}{2}.
        \end{split}

These dot products are generally different. The triangle area is instead :math:`S_0=\tfrac12\lVert\mathbf{r}_1^0\times\mathbf{r}_2^0\rVert`.

SETTLE
------

The triangle of constraints case is very common for the biomolecular simulations since most of water models constrain both bonds between hydrogens and oxygen and H-O-H angle. And there is an analytical solution for this specific case, called SETTLE :cite:`Constraints-MiyamotoJCC92`. This algorithm is in fact used in most of the MD software. In brief, the idea of SETTLE is first to reduce the problem into two dimensions, and than compute the rotation angle of the constraints triangle using two-dimensional Euler angles.

More than three coupled constraints
-----------------------------------

We can keep expanding our formulas to the cases with more coupled constraints, and they will get more and more complex. However, there are only a few cases that is useful to solve this way: four constraints coupled to the central (e.g. methane molecule) atom and various rings that should be kept planar (e.g. benzene ring in phenylalanine or double ring in the tryptophane side-chain). These are not very common, the rings are usually kept planar with improper dihedral potentials, so we will not derive them here. Much more important is the case where all covalent bonds are constraints, where universal formalism is needed. The two most common algorithms for this case are SHAKE :cite:`Constraints-RyckaertJCP77` (and its "velocity version" RATTLE :cite:`Constraints-AndersenJCP83`) and LINCS :cite:`Constraints-HessJCC97,Constraints-HessJCTC07`. Note that the latter is the default option in GROMACS.

.. bibliography::
    :keyprefix: Constraints-
