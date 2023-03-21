Numerical integration of equations of motion 
============================================

Basic integration scheme
------------------------

We start with Newtonian equations of motion, which state that at any given time :math:`t`:

    .. math::

        m_i\frac{d\mathbf{v}_{i}}{dt}=\mathbf{f}_{i}
        
        \frac{d\mathbf{r}_i}{dt}=\mathbf{v}_{i}

Here, :math:`i` is the particle index, :math:`m_i` is its mass, :math:`\mathbf{r}_i` its position, :math:`\mathbf{v}_i` its velocity and the force acting on a particle due to interactions with other particles is :math:`\mathbf{f}_i`.


To integrate these equations numerically, we need to introduce a small time increment :math:`\Delta t`, called timestep.
If :math:`\Delta t` is infinitely small, we can estimate the derivatives of :math:`\mathbf{v}_i` and :math:`\mathbf{r}_i` with respect to time as:

    .. math::

        \frac{d\mathbf{r}_{i}}{dt} = \lim_{\Delta t \rightarrow\infty} \frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t)}{\Delta t}

        \frac{d\mathbf{v}_{i}}{dt} = \lim_{\Delta t \rightarrow\infty} \frac{\mathbf{v}_{i}(t+\Delta t)-\mathbf{v}_{i}(t)}{\Delta t}

Obviously, numerical integration assumes that :math:`\Delta t` is finitely small.
Otherwise we will have to compute infinite number of steps, which is impossible.
So, we have to agree on approximation:

    .. math::

        \frac{d\mathbf{r}_{i}}{dt} \approx \frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t)}{\Delta t}

        \frac{d\mathbf{v}_{i}}{dt} \approx \frac{\mathbf{v}_{i}(t+\Delta t)-\mathbf{v}_{i}(t)}{\Delta t}

Plugging these into Newtonian equations gives:

    .. math::

        m_i\frac{\mathbf{v}_{i}(t+\Delta t)-\mathbf{v}_{i}(t)}{\Delta t} = \mathbf{f}_{i}(t)

        \frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t)}{\Delta t} = \mathbf{v}_{i}(t)

Or:

    .. math::

        \mathbf{v}_{i}(t+\Delta t) = \mathbf{v}_{i}(t) + \frac{\mathbf{f}_{i}(t)}{m_i}\Delta t

        \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t)\Delta t

Since functions are represented by numbers on a computer and to make these more readable, let us introduce the following notations: :math:`\mathbf{x}_{i}(t)=\mathbf{x}_{i}^n`, :math:`\mathbf{v}_{i}(t)=\mathbf{v}_{i}^n` and :math:`\mathbf{f}_{i}(t)=\mathbf{f}_{i}^n`.
Here, :math:`n` is the number of the timestep, which we will notate as :math:`tau = \Delta t` (i.e. current time is :math:`t_n=n\times\tau`):

    .. math::

        \mathbf{v}_{i}^{n+1} = \mathbf{v}_{i}^{n} + \frac{\mathbf{f}_{i}^{n}}{m_i}\tau

        \mathbf{r}_{i}^{n+1} = \mathbf{r}_{i}^{n} + \mathbf{v}_{i}^{n}\tau

Which is the simplest explicit numerical integration scheme.

Velocity-verlet algorithm
-------------------------

Using Taylor expansion:

    .. math::

        \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\frac{d\mathbf{r}_{i}(t)}{dt}\Delta t+\frac{1}{2}\frac{d^2\mathbf{r}_{i}(t)}{dt}\Delta t^{2}+\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})

        \mathbf{r}_{i}(t-\Delta t)=\mathbf{r}_{i}(t)-\frac{d\mathbf{r}_{i}(t)}{dt}\Delta t+\frac{1}{2}\frac{d^2\mathbf{r}_{i}(t)}{dt}\Delta t^{2}-\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})
        
If we plug Newtonian equations into relation above, we get:

    .. math::

        \mathbf{r}_{i}(t+\Delta t)=\mathbf{r}_{i}(t)+\mathbf{v}_{i}(t)\Delta t+\frac{\mathbf{f}_{i}(t)}{2m}\Delta t^{2}+\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})
        
        \mathbf{r}_{i}(t-\Delta t)=\mathbf{r}_{i}(t)-\mathbf{v}_{i}(t)\Delta t+\frac{\mathbf{f}_{i}(t)}{2m}\Delta t^{2}-\frac{1}{6}\frac{d^{3}\mathbf{r}_{i}(t)}{dt^{3}}\Delta t^{3}+\mathbf{o}(\Delta t^{4})

Adding these up:

    .. math::

        \mathbf{r}_{i}(t+\Delta t)=2\mathbf{r}_{i}(t)-\mathbf{r}_{i}(t-\Delta t)+\frac{\mathbf{f}_{i}(t)}{m}\Delta t^{2}+\mathbf{o}(\Delta t^{4}).

To compute velocities, one can use:

    .. math::

        \mathbf{v}_{i}(t)=\frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t-\Delta t)}{2\Delta t}+\mathbf{o}(\Delta t^{2}).

Or, in numerical notations:

    .. math::

        \mathbf{r}_{i}^{n+1}=2\mathbf{r}_{i}^{n}-\mathbf{r}_{i}^{n-1}+\frac{\mathbf{f}_{i}^{n}}{m}\tau^{2},

        \mathbf{v}_{i}^{n}=\frac{\mathbf{r}_{i}^{n+1}-\mathbf{r}_{i}^{n-1}}{2\tau}.
        

This is called Velocity-Verlet algorithm.

Leap-Frog algorithm
-------------------

Newtonian equation of motion are second order differential equation for positions.
By introducing velocities, we can re-write them as first-order equation.
Notice that in this case, the velocities and positions are never on both right- and left-hand side.
If we compute position on full steps and velocities on half steps, we can increase the approximation to the second order.
Note that this is only possible because (and when) the forces :math:`\mathbf{f}_{i}` are functions of coordinates.

We start with re-writing the basic scheme but using the central point to where the derivative is estimated:

    .. math::

        m_i\frac{\mathbf{v}_{i}(t+\frac{\Delta t}{2})-\mathbf{v}_{i}(t-\frac{\Delta t}{2})}{\Delta t} = \mathbf{f}_{i}(t)

        \frac{\mathbf{r}_{i}(t+\Delta t)-\mathbf{r}_{i}(t)}{\Delta t} = \mathbf{v}_{i}(t+\frac{\Delta t}{2})

Re-arranging the formulas:

    .. math::

        \mathbf{v}_{i}(t+\frac{\Delta t}{2}) = \mathbf{v}_{i}(t-\frac{\Delta t}{2}) + \frac{\mathbf{f}_{i}(t)}{m_i}\Delta t

        \mathbf{r}_{i}(t+\Delta t) = \mathbf{r}_{i}(t) + \mathbf{v}_{i}(t+\frac{\Delta t}{2})\Delta t

Or, in numerical notations:

    .. math::

        \mathbf{v}_{i}^{n+\frac{1}{2}} = \mathbf{v}_{i}^{n-\frac{1}{2}} + \frac{\mathbf{f}_{i}^{n}}{m_i}\tau

        \mathbf{r}_{i}^{n+1} = \mathbf{r}_{i}^{n} + \mathbf{v}_{i}^{n+\frac{1}{2}}\tau

Langevin equation
-----------------

Langevin equations of motion often used to describe coarse-grained systems, i.e. when more than one atom is considered a particle.
In addition to the deterministic force due to inter-particle interactions, there are two terms in these equations: drag force and random force.
Combined together, these imitate the effect on the medium on the system.
The functional form of Langevin equations are:

    .. math::

            m_i\frac{d^2\mathbf{r}_i}{dt^2}=-\nabla_iV-\lambda\frac{d\mathbf{r}_i}{dt}+\eta_i(t),

where :math:`m_i` is the mass of the particle, :math:`V({\mathbf{r}_i})` --- potential energy, :math:`\eta_i(t)` --- random force that is normally distributed:

    .. math::

            \langle\eta_i(t)\eta_j(t)\rangle = 2\lambda k_BT\delta_{ij}\delta(t-t')

Note that the Langevin equations are stochastic and can not be integrated in normal sense.
The numerical integration is possible though, with the help of random variables.

The coefficient :math:`\lambda` is called damping and describe how strong the interactions with the media are.
If :math:`\lambda` is big, the system becomes Brownian, when :math:`lambda` is small the system is Newtonian.
Take a look at the equation and its limits for small and big :math:`lambda`.
Typical value of the damping coefficient is :math:`\lambda=50`.

To integrate Langevin equations numerically one can use the following integration:

    .. math::

        \begin{split}
         m_i\frac{\mathbf{v}_{i}^{n+\frac{1}{2}}-\mathbf{v}_{i}^{n-\frac{1}{2}}}{\tau} &= \mathbf{f}_{i}^{n} - \lambda m_i\frac{\mathbf{v}_{i}^{n+\frac{1}{2}}+\mathbf{v}_{i}^{n-\frac{1}{2}}}{2}+\sqrt{\frac{2k_BT\lambda m_i}{\tau}}\mathbf{r}_i^f\\
        \frac{\mathbf{r}_{i}^{n+1}-\mathbf{r}_{i}^{n}}{\tau} &= \mathbf{v}_{i}^{n+\frac{1}{2}}

        \end{split}

Here, :math:`\mathbf{r}_i^f` --- a vector of three normally distributed random variables with dispersion :math:`1` and expectancy of :math:`0`. :math:`k_B=8.31\times10^{-3}` kJ/molK is Boltzman constant, :math:`T` is temperature.
The numerical scheme above is explicit.