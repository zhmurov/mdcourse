Implementation of a simple molecular dynamics code: NaCl in vacuum
==================================================================

Let us consider a basic example to show the underpinnings of any molecular dynamics simulations code.
We will create a simulation box of size 4nm x 4nm x 4nm and fill it with Na and Cl atoms.
The motion of atoms will be governed by Newton equations of motion and periodic boundary conditions will be applied to keep the atoms inside the box.
The atoms will interact with each other via Lennard-Jones and Coulomb potentials.

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

Before, we can to create and fill an array containing all the particles in the system, we need to get and define all the parameters for our simulations (temperature, masses, charges, etc).
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


    .. code::

        #include <stdio.h>
        #include <math.h>
        #include <limits.h>
        #include <vector>
        #include <random>

        #define N  100
        #define NSTEPS 1000000
        #define T  300.0 // K

        #define KB 8.314462e-3 // kJ/mol/K
        #define KC 138.9118    // nm*kJ/mol/e^2 

        #define L  5.0  // nm
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
            float currentTemperature = T;
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

            for (int n = 0; n < NSTEPS; n++)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        float dx = atoms[i].x - atoms[j].x;
                        float dy = atoms[i].y - atoms[j].y;
                        float dz = atoms[i].z - atoms[j].z;

                        // Check these
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