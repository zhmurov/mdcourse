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