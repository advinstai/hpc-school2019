/* 
 * simple lennard-jones potential MD code with velocity verlet.
 * units: Length=Angstrom, Mass=amu; Energy=kcal
 *
 * optimized c version with threads.
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

#if defined(_OPENMP)
#include <omp.h>
#endif

/* generic file- or pathname buffer length */
#define BLEN 200

/* a few physical constants */
const double kboltz=0.0019872067;     /* boltzman constant in kcal/mol/K */
const double mvsq2e=2390.05736153349; /* m*v^2 in kcal/mol */

/* structure to hold the complete information 
 * about the MD system */
struct _mdsys {
    double dt, mass, epsilon, sigma, box, rcut;
    double ekin, epot, temp, __pad1;
    double *rx, *ry, *rz;
    double *vx, *vy, *vz;
    double *fx, *fy, *fz;
    int natoms,nfi,nsteps,nthreads;
};
typedef struct _mdsys mdsys_t;

/* helper function: read a line and then return
   the first string with whitespace stripped off */
static int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}
 
/* helper function: apply minimum image convention */
__attribute__((pure))
static double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= boxby2 + boxby2;
    while (x < -boxby2) x += boxby2 + boxby2;
    return x;
}

/* compute kinetic energy */
static void ekin(mdsys_t *sys)
{   
    double ekin=0.0;
    const int natoms = sys->natoms;
    const double *vx = sys->vx;
    const double *vy = sys->vy;
    const double *vz = sys->vz;

#if defined(_OPENMP)
#pragma omp parallel for reduction(+:ekin)
#endif
    for (int i=0; i < natoms; ++i) {
        ekin += vx[i]*vx[i];
        ekin += vy[i]*vy[i];
        ekin += vz[i]*vz[i];
    }
    sys->ekin = ekin*0.5*mvsq2e*sys->mass;
    sys->temp = 2.0*sys->ekin/(3.0*natoms-3.0)/kboltz;
}

/* compute forces */
static void force(mdsys_t *sys) 
{
    double epot;
    int i;
    const int natoms = sys->natoms;

    /* zero energy and forces */
    epot=0.0;
    memset(sys->fx,0,sizeof(double)*natoms*sys->nthreads);
    memset(sys->fy,0,sizeof(double)*natoms*sys->nthreads);
    memset(sys->fz,0,sizeof(double)*natoms*sys->nthreads);

    {
        double c12,c6,boxby2,rsq,rcsq;
        double *fx,*fy,*fz;
        const double *rx,*ry,*rz;

        c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
        c6  = 4.0*sys->epsilon*pow(sys->sigma, 6.0);
        rcsq= sys->rcut*sys->rcut;
        boxby2 = 0.5*sys->box;

        rx = sys->rx; ry = sys->ry; rz = sys->rz;

#if defined(_OPENMP)
#pragma omp parallel for private(i) reduction(+:epot)
#endif
        for(i=0; i < natoms-1; ++i) {
            double dx,dy,dz;
#if defined(_OPENMP)
            int tid = omp_get_thread_num();
#else
            int tid = 0;
#endif
            fx = sys->fx + tid*natoms;
            fy = sys->fy + tid*natoms;
            fz = sys->fz + tid*natoms;
            for(int j=i+1; j < natoms; ++j) {

                /* get distance between particle i and j */
                dx=pbc(rx[i] - rx[j], boxby2);
                dy=pbc(ry[i] - ry[j], boxby2);
                dz=pbc(rz[i] - rz[j], boxby2);
                rsq = dx*dx + dy*dy + dz*dz;
      
                /* compute force and energy if within cutoff */
                if (rsq < rcsq) {
                    double r6,rinv,ffac;

                    rinv  = 1.0/rsq;
                    r6    = rinv*rinv*rinv;
                    ffac  = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                    epot += r6*(c12*r6 - c6);

                    fx[i] += dx*ffac;
                    fy[i] += dy*ffac;
                    fz[i] += dz*ffac;
                    fx[j] -= dx*ffac;
                    fy[j] -= dy*ffac;
                    fz[j] -= dz*ffac;
                }
            }
        }
    }
    sys->epot = epot;
#if defined(_OPENMP)
#pragma omp parallel for private(i)
#endif
    for (i=0; i < natoms; ++i) {
        for (int j = 1; j < sys->nthreads; ++j) {
            sys->fx[i] += sys->fx[i+j*natoms];
            sys->fy[i] += sys->fy[i+j*natoms];
            sys->fz[i] += sys->fz[i+j*natoms];
        }
    }
}

/* velocity verlet */
static void velverlet(mdsys_t *sys)
{
    int i=0;
    double dt = sys->dt;
    double dtfm = 0.5*dt / mvsq2e / sys->mass;
    const int natoms = sys->natoms;

    /* first part: propagate velocities by half and positions by full step */
#if defined(_OPENMP)
#pragma omp parallel for private(i)
#endif
    for (i=0; i < natoms; ++i) {
        sys->vx[i] += dtfm * sys->fx[i];
        sys->vy[i] += dtfm * sys->fy[i];
        sys->vz[i] += dtfm * sys->fz[i];
        sys->rx[i] += dt*sys->vx[i];
        sys->ry[i] += dt*sys->vy[i];
        sys->rz[i] += dt*sys->vz[i];
    }

    /* compute forces and potential energy */
    force(sys);

    /* second part: propagate velocities by another half step */
#if defined(_OPENMP)
#pragma omp parallel for private(i)
#endif
    for (i=0; i < natoms; ++i) {
        sys->vx[i] += dtfm * sys->fx[i];
        sys->vy[i] += dtfm * sys->fy[i];
        sys->vz[i] += dtfm * sys->fz[i];
    }
}

/* append data to output. */
static void output(mdsys_t *sys, FILE *erg, FILE *traj)
{
    int i;
    
    printf("% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(erg,"% 8d % 20.8f % 20.8f % 20.8f % 20.8f\n", sys->nfi, sys->temp, sys->ekin, sys->epot, sys->ekin+sys->epot);
    fprintf(traj,"%d\n nfi=%d etot=%20.8f\n", sys->natoms, sys->nfi, sys->ekin+sys->epot);
    for (i=0; i<sys->natoms; ++i) {
        fprintf(traj, "Ar  %20.8f %20.8f %20.8f\n", sys->rx[i], sys->ry[i], sys->rz[i]);
    }
}


/* main */
int main(int argc, char **argv) 
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;

#if defined(_OPENMP)
    sys.nthreads = omp_get_max_threads();
#else
    sys.nthreads = 1;
#endif

    /* read input file */
    if(get_a_line(stdin,line)) return 1;
    sys.natoms=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.mass=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.epsilon=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.sigma=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.rcut=atof(line);
    if(get_a_line(stdin,line)) return 1;
    sys.box=atof(line);
    if(get_a_line(stdin,restfile)) return 1;
    if(get_a_line(stdin,trajfile)) return 1;
    if(get_a_line(stdin,ergfile)) return 1;
    if(get_a_line(stdin,line)) return 1;
    sys.nsteps=atoi(line);
    if(get_a_line(stdin,line)) return 1;
    sys.dt=atof(line);
    if(get_a_line(stdin,line)) return 1;
    nprint=atoi(line);

    /* allocate memory */
    sys.rx=(double *)malloc(sys.natoms*sizeof(double));
    sys.ry=(double *)malloc(sys.natoms*sizeof(double));
    sys.rz=(double *)malloc(sys.natoms*sizeof(double));
    sys.vx=(double *)malloc(sys.natoms*sizeof(double));
    sys.vy=(double *)malloc(sys.natoms*sizeof(double));
    sys.vz=(double *)malloc(sys.natoms*sizeof(double));
    sys.fx=(double *)malloc(sys.natoms*sizeof(double)*sys.nthreads);
    sys.fy=(double *)malloc(sys.natoms*sizeof(double)*sys.nthreads);
    sys.fz=(double *)malloc(sys.natoms*sizeof(double)*sys.nthreads);

    /* read restart */
    fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
        }
        for (i=0; i<sys.natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
        }
        fclose(fp);
        memset(sys.fx, 0, sizeof(double)*sys.natoms);
        memset(sys.fy, 0, sizeof(double)*sys.natoms);
        memset(sys.fz, 0, sizeof(double)*sys.natoms);
    } else {
        perror("cannot read restart file");
        return 3;
    }

    /* initialize forces and energies.*/
    sys.nfi=0;
    force(&sys);
    ekin(&sys);
    
    erg=fopen(ergfile,"w");
    traj=fopen(trajfile,"w");

    printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
#if defined(_OPENMP)
    printf("Using %d threads\n",sys.nthreads);
#endif
    printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
    output(&sys, erg, traj);

    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {

        /* write output, if requested */
        if ((sys.nfi % nprint) == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
        velverlet(&sys);
        ekin(&sys);
    }
    /**************************************************/

    /* clean up: close files, free memory */
    printf("Simulation Done.\n");
    fclose(erg);
    fclose(traj);

    free(sys.rx);
    free(sys.ry);
    free(sys.rz);
    free(sys.vx);
    free(sys.vy);
    free(sys.vz);
    free(sys.fx);
    free(sys.fy);
    free(sys.fz);

    return 0;
}
