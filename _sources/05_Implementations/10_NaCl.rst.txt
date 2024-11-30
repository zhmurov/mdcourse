Implementation of a simple molecular dynamics code: NaCl in vacuum
==================================================================

Let us consider a basic example to show the underpinnings of any molecular dynamics simulations code.
We will create a simulation box of size 4nm x 4nm x 4nm and fill it with Na and Cl atoms.
The motion of atoms will be governed by Newton equations of motion and periodic boundary conditions will be applied to keep the atoms inside the box.
The atoms will interact with each other via Lennard-Jones and Coulomb potentials.

Parameters and units
--------------------

Before, we can create and fill an array containing all the particles in the system, we need to get and define all the parameters for our simulations (temperature, masses, charges, etc).
This marks an important decision to make: we need to choose the units.
It is tempting to use familiar SI system, but this will lead to a very small numbers --- atoms are very small.
We don't need to invent anything here though, since there are two system that are widely used in Molecular Dynamics.
Here, we are going to use the same units as MD software GROMACS: nm for length, ps for time, g/mol for mass, K for temperature.
These units conveniently converge to kJ/mol for energy.
We will also use the 

    .. exercise::

        Suppose that you are using nm for length, ps for time, g/mol for mass, K for temperature and elementary charge for charge.
        Using Newton equation of motion, show that the units of energy are kJ/mol.
        Using Maxwell-Boltzmann distribution for particle velocities, find a numerical value and units for Boltzmann constant.
        Using equation for Coulomb energy for two interacting particles, define the units and find numerical value for Coulomb constant.

With the units chosen, it is quite easy to define masses and charges for the particles.
The charges are :math:`+1` for sodium ions (Na) and :math:`-1` for chloride (Cl)
Masses can be taken from periodic table of elements: :math:`23.0` g/mol for sodium and :math:`35.5` g/mol for chloride.

We also need to define Van-der-Waals parameters (:math:`\sigma` and :math:`\varepsilon`) for all three possible interactions: Na-Na, Cl-Cl and Na-Cl.
In force-field, these parameters are usually listed on atom base: one set per atom type.
When the forces and energies are computed, these are averaged using arithmetic average for :math:`sigma` and geometric average for :math:`\varepsilon`.
The main reason for defining parameters this way is to minimize the number of parameters needed for an atomic system: if we have :math:`N_{at}` atom types, we will have :math:`N_{par}=N_{at}(N_{at}-1)/2+N_{at}` different pairs.
Though this is not a problem in our case (:math:`N_{par}=2\times1/2+2=3`), we are going to stick with this strategy to able to use the force-field files as a source of parameters.
If you did not do so already, `download <http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-jul2022.ff.tgz>`_ and unpack the force-field files.
The file that we need is called ``ffnonbonded.itp``.
There we should look for sodium and chloride, the respective lines are:

    .. code::

        ...
        ; type atnum         mass   charge ptype           sigma  epsilon
        ...
        SOD    11    22.989770    0.000     A  0.251367073323  0.1962296  ; Sodium Ion
        ...
        CLA    17    35.450000    0.000     A  0.404468018036  0.6276000  ; Chloride Ion
        ...

Two parameters we are looking for are in the last two columns, which list :math:`sigma` and :math:`epsilon` for the respective atoms in GROMACS units.
Note that the masses are also listed there (with more precision), feel free to take those from there as well.
The charges on the other hand are not in this list, since they are more volatile with respect to local environment --- same atom type can have different charges depending on the atom position in the molecule.
The charges thus are defined when the molecule is constructed as we saw when we were making topology of a small molecule.

Let us summarize all the parameters we found in the code, by adding this on top of the file:

    .. code::

        #define T  300.0 // K

        #define KB 8.314462e-3 // kJ/mol/K
        #define KC 138.9118    // nm*kJ/mol/e^2 

        #define Q_NA  1.0 // Elementary charge
        #define Q_CL  -1.0

        #define M_NA  23.0 // g/mol
        #define M_CL  35.5

        #define SIGMA_NA 0.25
        #define SIGMA_CL 0.40

        #define EPSILON_NA 0.196
        #define EPSILON_CL 0.628

Creating the system and saving output for visualization
-------------------------------------------------------

Initially, we will place particle in the box randomly.
This may cause a problem if the atoms will be placed too close to one another, so we should expect that.
We will use basic XYZ file format to save output so that we will be able to visualize it in third-party software (e.g. VMD).
But to start with let use introduce the type ``Atom`` which we will use to store all the atom-related data (positions, velocities, forces, mass, etc.):

    .. code::

        struct Atom
        {
            float x, y, z;
            float vx, vy, vz;
            float fx, fy, fz;
            float q;
            float m;
            float sigma;
            float epsilon;
            char name;
        };

Note that from the performance standpoint this is not the most efficient way to store positions and other data that we will be frequently interacting with.
It would be better if these were more localized, for more efficient memory access pattern.
But this should be fine in our example.

Before we place the atom, we need to define the dimensions of the system.
Let us use quadratic box of size L, filled with N randomly placed atoms.
First, we define these values.

    .. code::

        #define N  100
        #define L  5.0  // nm

Now let us create the ``main`` function in which we define parameters for each atom and fill the box with randomly placed atoms.
Half of the atoms will be sodium ions, half --- chloride.

    .. code::

        int main(int argc, char* argv[])
        {
            std::vector<Atom> atoms(N);

            std::random_device randomDevice;
            std::mt19937 randomGenerator(randomDevice());
            std::uniform_real_distribution<> distributionX(0, L);

            for (int i = 0; i < N; i++)
            {
                atoms[i].x = distributionX(randomGenerator);
                atoms[i].y = distributionX(randomGenerator);
                atoms[i].z = distributionX(randomGenerator);

                if (i < N/2)
                {
                    atoms[i].name = 'N';
                    atoms[i].q = Q_NA;
                    atoms[i].m = M_NA;
                    atoms[i].sigma = SIGMA_NA;
                    atoms[i].epsilon = EPSILON_NA;
                }
                else
                {
                    atoms[i].name = 'C';
                    atoms[i].q = Q_CL;
                    atoms[i].m = M_CL;
                    atoms[i].sigma = SIGMA_CL;
                    atoms[i].epsilon = EPSILON_CL;
                }

            }
        }

In this function, we first define random numbers generator using Mersene-Twister PRNG.
We need uniform distribution from 0 to L to draw coordinates of atoms from.
In a loop, we also add a condition on which atom the ``i``-th atom is going to be: sodium or chloride, and assign the parameters accordingly.

We also need to define initial velocities and zero the initial forces.
Though, we are going to have an excess in the potential energy in our initial state, the temperature control will also be needed.
Otherwise, the temperature will increase from its initial value.
We will deal with this on the later stages.

To instantiate the velocities, we need normally distributed random numbers.
The mean of this distribution should be zero, and the dispersion should be :math:`\sqrt{k_BT/m}`.
Since the masses of Na and Cl ions are different, we are going to define the normal distribution with dispersion of :math:`\sqrt{k_BT}` and then divide the number over :math:`\sqrt{m}` for each particle.
This way we can use one distribution for all particles:

    .. code::

        int main(int argc, char* argv[])
        {
            ...
            std::uniform_real_distribution<> distributionX(0, L);
            std::normal_distribution<> distributionV(0.0, sqrtf(KB*T));
            ...
            for (int i = 0; i < N; i++)
            {
                ...

                float mult = 1.0/sqrtf(atoms[i].m);
                atoms[i].vx = mult*distributionV(randomGenerator);
                atoms[i].vy = mult*distributionV(randomGenerator);
                atoms[i].vz = mult*distributionV(randomGenerator);

                atoms[i].fx = 0.0;
                atoms[i].fy = 0.0;
                atoms[i].fz = 0.0;
            }
        }

Now, all the positions are defined, we want to be able to save and visualize them.
To do so, let us define a function that will save the positions of the atoms in XYZ format:

    .. code::

        void saveFrame(const char* filename, const char* modifier, std::vector<Atom> atoms)
        {
            FILE* out = fopen(filename, modifier);
            fprintf(out, "%ld\nNa+Cl\n", atoms.size());
            for (int i = 0; i < atoms.size(); i++)
            {
                fprintf(out, "%c\t%f\t%f\t%f\n",
                    atoms[i].name, 
                    atoms[i].x*10.0,
                    atoms[i].y*10.0,
                    atoms[i].z*10.0);
            }
            fclose(out);
        }

The function takes three arguments: the name of the file to save data into, modifier (``w`` to overwrite or ``a`` to append data to the end of the file) and the vector of atoms containing the data that we are saving.
Note that XYZ format allows one to save multiple frames by appending the data to the same file.
Also note, that XYZ uses angstroms as length units, hence we are multiplying the positions by 10.
The first line of the file should contain the number of atoms, which is the length of the atoms vector.
The second line should have the description of the system.
These two lines are repeated for each frame if we are appending them.
XYZ uses only one character for atom name, which we defined earlier, while filling the vector of atoms.

Now we can add the call to this function to write the initial frame:

    .. code::

        saveFrame("nacl.xyz", "w", atoms);

If you compile and execute the code, you should have ``nacl.xyz`` file, which can be visualized with e.g. VMD software.
It is better to look at in VdW representation, that shows Van-der-Waals spheres for particles.

Integrating equations of motion and applying periodic boundary conditions
-------------------------------------------------------------------------

We have the static picture of the system.
Let us now make the atoms move by integrating their equations of motion numerically.
Before implementing the numerical integration, we need to define its parameters: the timestep and the number of steps:

    .. code::

        #define NSTEPS 1000000
        #define tau 0.001  // ps

We used the timestep of 1 fs, which should be plenty small: we already showed that Leap-Frog method is stable with this timestep even if we have fast harmonic bond oscillations.
Now, let us add a main loop in which we are going to do these integration steps.
On each step, for each atom, we need to re-compute its velocities and positions according to the integration scheme.
We will use the Leap-Frog algorithm:

    .. math::

        \mathbf{v}_{i}^{n+\frac{1}{2}} = \mathbf{v}_{i}^{n-\frac{1}{2}} + \frac{\mathbf{f}_{i}^{n}}{m_i}\tau

        \mathbf{r}_{i}^{n+1} = \mathbf{r}_{i}^{n} + \mathbf{v}_{i}^{n+\frac{1}{2}}\tau

In the code, this renders to:

    .. code::

        for (int n = 0; n < NSTEPS; n++)
        {
            for (int i = 0; i < N; i++)
            {
                float mult = tau/atoms[i].m;

                atoms[i].vx = atoms[i].vx + mult*atoms[i].fx;
                atoms[i].vy = atoms[i].vy + mult*atoms[i].fy;
                atoms[i].vz = atoms[i].vz + mult*atoms[i].fz;

                atoms[i].x = atoms[i].x + tau*atoms[i].vx;
                atoms[i].y = atoms[i].y + tau*atoms[i].vy;
                atoms[i].z = atoms[i].z + tau*atoms[i].vz;
            }
        }

We also want to save the positions of atoms every now and then (say, every 100 steps).
To do so, we can use the function we implemented earlier by calling it with ``a`` modifier:

    .. code::

        if (n % 100 == 0)
        {
            saveFrame("nacl.xyz", "a", atoms);
        }

Note that this should be placed inside the loop over the time steps, but not inside the loop over the particles.

If one to compile and run the code, they will see that particles are floating away from the initial box they we put in.
This is because we are not using any boundary conditions yet.
In molecular dynamics, the most common boundary conditions are Periodic Boundary Conditions or PBC.
These affect system in two ways.
First, when particle floats from the box, it is moved as if it enters from the other side.
Second, when the distances are computed, one should select the distance between the nearest images of the particles.
We will deal with the second part later, let us do the first.

What we need to do is transfer the particle to the other side if the box if it crosses the border.
Since our box is quadratic, we can do this one component at a time.
There are two options to go beyond the border then: when the component of the position of the particle becomes negative, or when it is larger than the size of the box L:

    .. code::

        float transferPBC(float x)
        {
            if (x < 0)
            {
                return x + L;
            }
            else if (x > L)
            {
                return x - L;
            }
            return x;
        }

This method only works when particles moved no further than one unit cell (i.e. crossed the border only once).
There are more robust methods, but in the interest of clarity we stick with this simple method for now.
Now we need to call this function during the simulation so that the particles are transferred to the initial box.
Because of the way we implemented the function, is should be called for every component of the position for every atom, which can be done right after moving atoms in the integrator:

    .. code::

        atoms[i].x = transferPBC(atoms[i].x);
        atoms[i].y = transferPBC(atoms[i].y);
        atoms[i].z = transferPBC(atoms[i].z);

Compile and run the code to see that everything worked.
The atoms should now be in the same box throughout the simulation.

Adding interactions between atoms
---------------------------------

The atoms should be interacting via Van-der-Waals and electrostatic interactions.
The functional form of the interaction potential is:

    .. math::

        V_{nb} = \sum_{i,j}\left(\varepsilon_{ij}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{12}-2\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\right] + k_c\frac{q_{i}q_{j}}{r_{ij}}\right)

To get the forces out of this, we need to compute the gradient of the potential with respect to the coordinates of the particle in question.
The formula for the interatomic force acting on the ``i``-th particle is:

    .. math::

        \mathbf{f}_i = \sum_{j}\left(-12\varepsilon_{ij}\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}\left[\left(\frac{\sigma_{ij}}{r_{ij}}\right)^{6}-1\right]\frac{1}{r_{ij}^2}\frac{\mathbf{r}_{ij}}{r_{ij}}-k_c\frac{q_{i}q_{j}}{r_{ij}^2}\frac{\mathbf{r}_{ij}}{r_{ij}} \right)

    .. exercise::

        Derive the expression above.

Here, :math:`\mathbf{r}_{ij}` is the vector, connecting ``i``-th and ``j``-th particles and ::math:`r_{ij}` is its length.

    .. code::

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                float dx = atoms[i].x - atoms[j].x;
                float dy = atoms[i].y - atoms[j].y;
                float dz = atoms[i].z - atoms[j].z;

                float dr2 = dx*dx + dy*dy + dz*dz;
                float sigma = 0.5*(atoms[i].sigma + atoms[j].sigma);
                float epsilon = sqrtf(atoms[i].epsilon*atoms[j].epsilon);
                float sigma2 = sigma*sigma;
                float sor2 = sigma2/dr2;
                float sor6 = sor2*sor2*sor2;
                float df = 12.0*epsilon*sor6*(sor6 - 1.0)/dr2;

                float dr = sqrtf(dr2);
                df += KC*atoms[i].q*atoms[j].q/(dr2*dr);

                atoms[i].fx += df*dx;
                atoms[i].fy += df*dy;
                atoms[i].fz += df*dz;

                atoms[j].fx -= df*dx;
                atoms[j].fy -= df*dy;
                atoms[j].fz -= df*dz;
            }
        }

Note that in the code above we took advantage of the Newton third law and compute force between two particle once and add the increment to both particles force but with opposite signs.
This requires us to use separate loop over atoms: all the forces should be computed before we integrate the equations of motion.

Before we test the code, there are two thing to fix.
First, the vector, connecting particles ``i`` and ``j`` and its length should be computed with periodic boundary conditions in mind: the closest distance between any periodic image of particle ``i`` and any periodic image of particle ``j`` should be taken.
Likely, this is easier than going through all the images: from every component of the vector :math:`\mathbf{r}_{ij}`, we should subtract the integer number of the box length in it.
To fix this, edit the code:

    .. code::

        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < i; j++)
            {
                float dx = atoms[i].x - atoms[j].x;
                float dy = atoms[i].y - atoms[j].y;
                float dz = atoms[i].z - atoms[j].z;

                dx -= rint(dx/L)*L;
                dy -= rint(dy/L)*L;
                dz -= rint(dz/L)*L;

                ...
            }
        }

Secondly, we are accumulating forces for each particle.
To make sure that we don't keep these increments from the previous steps, we need to set forces to zero after the numerical integration timestep is done:

    .. code::

        atoms[i].fx = 0.0;
        atoms[i].fy = 0.0;
        atoms[i].fz = 0.0;


Monitoring and controlling temperature of the system
----------------------------------------------------

We would expect that after running the simulation for a while, we will have a partly of fully formed NaCl crystal.
However, this is not what is happening.
The reason is following.
When we assigned the initial velocities to the atoms, we used Maxwell-Boltzmann distribution.
These velocities correspond to the temperature of 300K and should give us the kinetic energy that also corresponds to 300K.
However, we did not do the same with positions and potential energy --- we just placed the particles randomly in the selected compartment.
Because of that, the initial configuration of the system (i.e. randomly placed particles) has an excess in potential energy.
During the simulations, the extra energy makes its way into the kinetic energy, particles accelerate and never form crystal.
The effective temperature is just too high.

Let us check if this is the case.
We can use Maxwell-Boltzmann distribution again, only now we are going to get the actual temperature of the system out of particles velocities.
The formula that we are going to use is:

    .. math::

        T=\langle \frac{m_iv_i^2}{3k_B}\rangle

Where :math:`m_i` and :math:`v_i` are the mass and velocity of the particle, :math:`k_B` is Boltzmann constant.
The angular brackets represent average over all particles.
Since our molecular system is small, we also are going to take the average over time.
This will allow us to get better average with less noise.

To do the averaging, we are going to need two numbers: a real number to accumulate the sum of all the temperatures for all particles across the system and an integer that will keep the quantity of the numbers added.
The integer is not required since we can always compute it as a product of the number of particles and number of steps that we were accumulating.
But to avoid mistakes it is better to add one to this number every time we accumulate.
We are also adding the variable to keep last computed temperature, which is going to be useful later.
Initially, we set it to ``T`` since this is the temperature that corresponds to the initial velocities of the particles.

    .. code::

        double temperature = 0.0;
        int nTemperature = 0;
        float currentTemperature = T;

Note that we used ``double`` here.
We are going to add small number to a relatively large sum.
It is better to use double precision in such cases to avoid round-up errors.

Now, on each timestep and for each particle, we are going to add :math:`m_iv_i^2` to ``temperature`` and add one to ``nTemperature``:

    .. code::

        float v2 = atoms[i].vx*atoms[i].vx + 
        atoms[i].vy*atoms[i].vy +
        atoms[i].vz*atoms[i].vz;

        temperature += atoms[i].m*v2;
        nTemperature ++;

To compute and print the temperature:

    .. code::

        currentTemperature = (temperature/(3.0*KB))/nTemperature;
        printf("%d\t%f\n", n, currentTemperature);
        temperature = 0.0;
        nTemperature = 0;

Note that we divide by :math:`3k_B` on the later stage in order to do the division only once, which is better from the performance standpoint.

If we ran the code, we will see that the temperature rises from its initial value of 300K, which is expected: we have too much extra potential energy in our initial conformation.
In order to control the current temperature, we are going to employ basic thermostat, that will scale the velocities of particles on each timestep.
Since the computed temperature is not exact and will oscillate, we are going to smooth the scaling by employing the relaxation time constant, :math:`\tau_r`.
With this constant, the scaling factor is:

    .. math::

        s_v=\sqrt{1 - \frac{(T_c-T)*\tau}{T*\tau_r}}

Here, :math:`T_c` and :math:`T` are our current (computed) and target temperatures, :math:`\tau` is the timestep and :math:`\tau_r` is the relaxation time.
Note that the relaxation time should be larger than the timestep and larger than periods with which we update the current temperature.
In the code, the temperature control renders to the definition of relaxation time and scaling of the velocities when we integrate:

    .. code::

        ...
        #define relax 10.0
        ...
        float scale = sqrtf(1.0 - ((currentTemperature-T)*tau)/(T*relax));
        atoms[i].vx *= scale;
        atoms[i].vy *= scale;
        atoms[i].vz *= scale;

Run the code to see that the temperature is leveling off on the target of 300K.

The final code
--------------

If one puts all the above together, they end up with something similar to the following code snippet:

    .. code::

        #include <stdio.h>
        #include <math.h>
        #include <limits.h>
        #include <vector>
        #include <random>

        #define N  100
        #define L  5.0  // nm
    
        #define T  300.0 // K

        #define KB 8.314462e-3 // kJ/mol/K
        #define KC 138.9118    // nm*kJ/mol/e^2 

        #define NSTEPS 1000000
        #define tau 0.001  // ps

        #define Q_NA  1.0 // Elementary charge
        #define Q_CL  -1.0

        #define M_NA  23.0 // g/mol
        #define M_CL  35.5

        #define SIGMA_NA 0.25
        #define SIGMA_CL 0.40

        #define EPSILON_NA 0.196
        #define EPSILON_CL 0.628

        #define relax 10.0

        struct Atom
        {
            float x, y, z;
            float vx, vy, vz;
            float fx, fy, fz;
            float q;
            float m;
            float sigma;
            float epsilon;
            char name;
        };

        float transferPBC(float x)
        {
            if (x < 0)
            {
                return x + L;
            }
            else if (x > L)
            {
                return x - L;
            }
            return x;
        }

        void saveFrame(const char* filename, const char* modifier, std::vector<Atom> atoms)
        {
            FILE* out = fopen(filename, modifier);
            fprintf(out, "%ld\nNa+Cl\n", atoms.size());
            for (int i = 0; i < atoms.size(); i++)
            {
                fprintf(out, "%c\t%f\t%f\t%f\n",
                    atoms[i].name, 
                    atoms[i].x*10.0,
                    atoms[i].y*10.0,
                    atoms[i].z*10.0);
            }
            fclose(out);
        }

        int main(int argc, char* argv[])
        {
            std::vector<Atom> atoms(N);

            std::random_device randomDevice;
            std::mt19937 randomGenerator(randomDevice());
            std::uniform_real_distribution<> distributionX(0, L);
            std::normal_distribution<> distributionV(0.0, sqrtf(KB*T));

            for (int i = 0; i < N; i++)
            {
                atoms[i].x = distributionX(randomGenerator);
                atoms[i].y = distributionX(randomGenerator);
                atoms[i].z = distributionX(randomGenerator);

                if (i < N/2)
                {
                    atoms[i].name = 'N';
                    atoms[i].q = Q_NA;
                    atoms[i].m = M_NA;
                    atoms[i].sigma = SIGMA_NA;
                    atoms[i].epsilon = EPSILON_NA;
                }
                else
                {
                    atoms[i].name = 'C';
                    atoms[i].q = Q_CL;
                    atoms[i].m = M_CL;
                    atoms[i].sigma = SIGMA_CL;
                    atoms[i].epsilon = EPSILON_CL;
                }

                float mult = 1.0/sqrtf(atoms[i].m);
                atoms[i].vx = mult*distributionV(randomGenerator);
                atoms[i].vy = mult*distributionV(randomGenerator);
                atoms[i].vz = mult*distributionV(randomGenerator);

                atoms[i].fx = 0.0;
                atoms[i].fy = 0.0;
                atoms[i].fz = 0.0;
            }

            saveFrame("nacl.xyz", "w", atoms);

            double temperature = 0.0;
            int nTemperature = 0;
            float currentTemperature = T;

            for (int n = 0; n < NSTEPS; n++)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        float dx = atoms[i].x - atoms[j].x;
                        float dy = atoms[i].y - atoms[j].y;
                        float dz = atoms[i].z - atoms[j].z;

                        dx -= rint(dx/L)*L;
                        dy -= rint(dy/L)*L;
                        dz -= rint(dz/L)*L;

                        float dr2 = dx*dx + dy*dy + dz*dz;
                        float sigma = 0.5*(atoms[i].sigma + atoms[j].sigma);
                        float epsilon = sqrtf(atoms[i].epsilon*atoms[j].epsilon);
                        float sigma2 = sigma*sigma;
                        float sor2 = sigma2/dr2;
                        float sor6 = sor2*sor2*sor2;
                        float df = 12.0*epsilon*sor6*(sor6 - 1.0)/dr2;

                        float dr = sqrtf(dr2);
                        df += KC*atoms[i].q*atoms[j].q/(dr2*dr);

                        atoms[i].fx += df*dx;
                        atoms[i].fy += df*dy;
                        atoms[i].fz += df*dz;

                        atoms[j].fx -= df*dx;
                        atoms[j].fy -= df*dy;
                        atoms[j].fz -= df*dz;
                    }
                }

                for (int i = 0; i < N; i++)
                {
                    float scale = sqrtf(1.0 - ((currentTemperature-T)*tau)/(T*relax));
                    atoms[i].vx *= scale;
                    atoms[i].vy *= scale;
                    atoms[i].vz *= scale;

                    float mult = tau/atoms[i].m;

                    atoms[i].vx = atoms[i].vx + mult*atoms[i].fx;
                    atoms[i].vy = atoms[i].vy + mult*atoms[i].fy;
                    atoms[i].vz = atoms[i].vz + mult*atoms[i].fz;

                    atoms[i].x = atoms[i].x + tau*atoms[i].vx;
                    atoms[i].y = atoms[i].y + tau*atoms[i].vy;
                    atoms[i].z = atoms[i].z + tau*atoms[i].vz;

                    atoms[i].x = transferPBC(atoms[i].x);
                    atoms[i].y = transferPBC(atoms[i].y);
                    atoms[i].z = transferPBC(atoms[i].z);

                    atoms[i].fx = 0.0;
                    atoms[i].fy = 0.0;
                    atoms[i].fz = 0.0;

                    float v2 = atoms[i].vx*atoms[i].vx + 
                            atoms[i].vy*atoms[i].vy +
                            atoms[i].vz*atoms[i].vz;

                    temperature += atoms[i].m*v2;
                    nTemperature ++;

                }

                if (n % 100 == 0)
                {
                    currentTemperature = (temperature/(3.0*KB))/nTemperature;
                    printf("%d\t%f\n", n, currentTemperature);
                    temperature = 0.0;
                    nTemperature = 0;
                    saveFrame("nacl.xyz", "a", atoms);
                }
            }
        }

Big problem with the code above
-------------------------------

There is one big problem with the code above, which can not be easily addressed.
When we compute the distance between two particles, we take into account the periodic boundary conditions.
This means, that any component of the :math:`r_{ij}` vector that connect two particles can not be larger than half of the simulation box.
What happens is when the distance becomes larger, we switch to another periodic image of the particle.
This creates an artificial mimima in the interaction potential on the half-box length.
Throughout the simulations, some of the particles can be dragged to this minima