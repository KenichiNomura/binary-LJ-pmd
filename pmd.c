/*----------------------------------------------------------------------
Program pmd.c performs parallel molecular-dynamics for Lennard-Jones
systems using the Message Passing Interface (MPI) standard.
----------------------------------------------------------------------*/
#include "pmd.h"

/*--------------------------------------------------------------------*/
int main(int argc, char **argv) {
    /*--------------------------------------------------------------------*/
    double cpu1;

    MPI_Init(&argc,&argv); /* Initialize the MPI environment */
    MPI_Comm_rank(MPI_COMM_WORLD, &sid);  /* My processor ID */
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  /* # of processors */

    init_params();
    make_tables();
    set_topology();
    init_conf();
    atom_copy();
    compute_accel(); /* Computes initial accelerations */

    analysis_manager(0);
    analysis_manager(1);

    cpu1 = MPI_Wtime();
    for (stepCount=currCount; stepCount<=currCount+StepLimit; stepCount++) {
        single_step();
        if (stepCount%StepAvg == 0) {
            eval_props();
            write_config(stepCount);

            analysis_manager(1);
        }
    }
    cpu = MPI_Wtime() - cpu1;
    if (sid == 0) printf("CPU & COMT = %le %le\n",cpu,comt);

    /* save last config */
    write_config(-1);

    analysis_manager(2);

    MPI_Finalize(); /* Clean up the MPI environment */
}

/*--------------------------------------------------------------------*/
void init_params() {
    /*----------------------------------------------------------------------
    Initializes parameters.
    ----------------------------------------------------------------------*/
    int a;
    FILE *fp;

    /* Read control parameters */
    fp = fopen("pmd.in","ro");
    fscanf(fp,"%d",&mdmode);
    fscanf(fp,"%d%d%d",&vproc[0],&vproc[1],&vproc[2]);
    fscanf(fp,"%d%d%d",&InitUcell[0],&InitUcell[1],&InitUcell[2]);
    fscanf(fp,"%le",&Density);
    fscanf(fp,"%le",&InitTemp);
    fscanf(fp,"%le",&DeltaT);
    fscanf(fp,"%d",&StepLimit);
    fscanf(fp,"%d",&StepAvg);
    fclose(fp);

    nproc=vproc[0]*vproc[1]*vproc[2];

    /* Vector index of this processor */
    vid[0] = sid/(vproc[1]*vproc[2]);
    vid[1] = (sid/vproc[2])%vproc[1];
    vid[2] = sid%vproc[2];

    /* Compute basic parameters */
    DeltaTH = 0.5*DeltaT;
    for (a=0; a<3; a++) al[a] = InitUcell[a]/pow(Density/4.0,1.0/3.0);
    for (a=0; a<3; a++) {
        gl[a]=al[a]*vproc[a];
        ol[a]=al[a]*vid[a];
    }
    if (sid == 0) {
        printf("mdmode: %d\n",mdmode);
        printf("Global MD box (gl): %9.2f %9.2f %9.2f\n",gl[0],gl[1],gl[2]);
        printf(" Local MD box (al): %9.2f %9.2f %9.2f\n",al[0],al[1],al[2]);
    }
    //printf("MD box origin (ol): %9.2f %9.2f %9.2f\n",ol[0],ol[1],ol[2]);

    /* Compute the # of cells for linked cell lists */
    for (a=0; a<3; a++) {
        lc[a] = al[a]/RCUT;
        rc[a] = al[a]/lc[a];
    }
    if (sid == 0) {
        printf("# of LL per domain (lc):  %d %d %d\n",lc[0],lc[1],lc[2]);
        printf("LL: rc = %e %e %e\n",rc[0],rc[1],rc[2]);
    }

}

/*--------------------------------------------------------------------*/
void set_topology() {
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/
    /* Integer vectors to specify the six neighbor nodes */
    int iv[6][3] = {
        {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
    };
    int ku,a,k1[3];

    /* Set up neighbor tables, nn & sv */
    for (ku=0; ku<6; ku++) {

        /* Vector index of neighbor ku */
        for (a=0; a<3; a++)
            k1[a] = (vid[a]+iv[ku][a]+vproc[a])%vproc[a];
        nns[ku] = k1[0]*vproc[1]*vproc[2]+k1[1]*vproc[2]+k1[2];

        /* Scalar neighbor ID, nn */
        nn[ku] = k1[0]*vproc[1]*vproc[2]+k1[1]*vproc[2]+k1[2];

        /* Vector index of neighbor ku */
        for (a=0; a<3; a++)
            k1[a] = (vid[a]-iv[ku][a]+vproc[a])%vproc[a];
        nnr[ku] = k1[0]*vproc[1]*vproc[2]+k1[1]*vproc[2]+k1[2];

        /* Shift vector, sv */
        for (a=0; a<3; a++) sv[ku][a] = al[a]*iv[ku][a];
    }

    /* Set up the node parity table, myparity */
    for (a=0; a<3; a++) {
        if (vproc[a] == 1)
            myparity[a] = 2;
        else if (vid[a]%2 == 0)
            myparity[a] = 0;
        else
            myparity[a] = 1;
    }
}

/*--------------------------------------------------------------------*/
void init_conf() {
    /*----------------------------------------------------------------------
    r are initialized to face-centered cubic (fcc) lattice positions.
    rv are initialized with a random velocity corresponding to Temperature.
    ----------------------------------------------------------------------*/
    double c[3],gap[3],e[3],vSum[3],gvSum[3],vMag;
    int j,a,nex,nX,nY,nZ;
    double seed;
    FILE *fp = NULL;
    long int *nary;
    int *tary;
    double *rary, *vary;
    /* FCC atoms in the original unit cell */
    double origAtom[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
        {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}
    };
    /* Atom type for multicomponent system */
    int origType[4] = {1, 1, 0, 0};

    if(mdmode!=0) {

        // only one MPIrank reads data, then broadcast the data to all other MPIranks
        if(sid==0) {
            fp = fopen("pmd.d","ro");
            fscanf(fp,"%d %d %lf %lf %lf\n",&nglob,&currCount,&gl[0],&gl[1],&gl[2]);
            currCount*=1000;
        }
        MPI_Bcast(&nglob, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&currCount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // allocate buffer to receive the broadcast
        nary=malloc(sizeof(long int)*nprocs);
        tary=malloc(sizeof(int)*nglob);
        rary=malloc(sizeof(double)*3*nglob);
        vary=malloc(sizeof(double)*3*nglob);

        // read data by MPIrank==0
        if(sid==0) {
            for(a=0; a<nprocs; a++) fscanf(fp,"%ld",&nary[a]);
            for(a=0; a<nglob; a++)
                fscanf(fp,"%d %lf %lf %lf %lf %lf %lf\n",
                       &tary[a], &rary[3*a],&rary[3*a+1],&rary[3*a+2],
                       &vary[3*a],&vary[3*a+1],&vary[3*a+2]);
        }

        // then, broadcast data from MPIrank==0
        MPI_Bcast(nary, nprocs, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(tary, nglob, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(rary, 3*nglob, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(vary, 3*nglob, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // get exscan to identify array indices to store
        n=nary[sid];
        nex=0;
        MPI_Exscan(&n,&nex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

        //printf("sid=%d, n=%d, nex=%d\n", sid, n, nex);
        if(sid==nprocs-1) {
            for(a=0; a<nprocs; a++) printf("n[%d]=%ld  ",a,nary[a]);
            printf("\n");
        }

        // resident atom information is stored between nex and nex+n
        for(j=0, a=nex; a<nex+n; a++, j++) {
            atype[j]=tary[a];
            r[j][0]=rary[3*a];
            r[j][1]=rary[3*a+1];
            r[j][2]=rary[3*a+2];
            rv[j][0]=vary[3*a];
            rv[j][1]=vary[3*a+1];
            rv[j][2]=vary[3*a+2];
        }

        // close file, free memory
        if(sid==0) fclose(fp);
        free(nary);
        free(tary);
        free(rary);
        free(vary);

    } else {

        if(sid==0) printf("Start from FCC config.\n");
        /* Set up a face-centered cubic (fcc) lattice */
        for (a=0; a<3; a++) gap[a] = al[a]/InitUcell[a];
        n = 0;
        for (nZ=0; nZ<InitUcell[2]; nZ++) {
            c[2] = nZ*gap[2];
            for (nY=0; nY<InitUcell[1]; nY++) {
                c[1] = nY*gap[1];
                for (nX=0; nX<InitUcell[0]; nX++) {
                    c[0] = nX*gap[0];
                    for (j=0; j<4; j++) {
                        for (a=0; a<3; a++)
                            r[n][a] = c[a] + gap[a]*origAtom[j][a];
                        atype[n]=origType[j];
                        n++;
                    }
                }
            }
        }

        /* Generate random velocities */
        seed = 13597.0+sid;
        vMag = sqrt(3*InitTemp);
        for(j=0; j<n; j++) {
            RandVec3(e,&seed);
            for (a=0; a<3; a++) {
                rv[j][a] = vMag*e[a];
            }
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);

    /* Total # of atoms summed over processors */
    MPI_Allreduce(&n,&nglob,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (sid == 0) printf("nglob = %d\n",nglob);

    for(a=0; a<3; a++) vSum[a] = 0.0;
    for(j=0; j<n; j++) {
        for (a=0; a<3; a++) {
            vSum[a] = vSum[a] + rv[j][a];
        }
    }

    MPI_Allreduce(vSum,gvSum,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    /* Make the total momentum zero */
    for (a=0; a<3; a++) gvSum[a] /= nglob;
    for (j=0; j<n; j++)
        for(a=0; a<3; a++) rv[j][a] -= gvSum[a];

    /* Initialize the initial coordinate */
    for (j=0; j<n; j++)
        for(a=0; a<3; a++)
        {
            r0[j][a] = r[j][a];
            v0[j][a] = rv[j][a];
        }
}

/*--------------------------------------------------------------------*/
void single_step() {
    /*----------------------------------------------------------------------
    r & rv are propagated by DeltaT using the velocity-Verlet scheme.
    ----------------------------------------------------------------------*/
    int i,a;

    half_kick(); /* First half kick to obtain v(t+Dt/2) */
    for (i=0; i<n; i++) /* Update atomic coordinates to r(t+Dt) */
        for (a=0; a<3; a++) r[i][a] = r[i][a] + DeltaT*rv[i][a];
    atom_move();
    atom_copy();
    compute_accel(); /* Computes new accelerations, a(t+Dt) */
    half_kick(); /* Second half kick to obtain v(t+Dt) */
}

/*--------------------------------------------------------------------*/
void half_kick() {
    /*----------------------------------------------------------------------
    Accelerates atomic velocities, rv, by half the time step.
    ----------------------------------------------------------------------*/
    int i,a;
    for (i=0; i<n; i++)
        for (a=0; a<3; a++) rv[i][a] = rv[i][a]+DeltaTH*ra[i][a];
}

/*--------------------------------------------------------------------*/
void compute_accel() {
    /*----------------------------------------------------------------------
    Given atomic coordinates, r[0:n+nb-1][], for the extended (i.e.,
    resident & copied) system, computes the acceleration, ra[0:n-1][], for
    the residents.
    ----------------------------------------------------------------------*/
    int i,j,a,lc2[3],lcyz2,lcxyz2,mc[3],c,mc1[3],c1;
    int bintra;
    double dr[3],rr,rrCut,lpe,dr2,dr2i;
    int idiv,ia,ib;
    double rdiv,vr,fr;

    /* Reset the potential & forces */
    lpe = 0.0;
    for (i=0; i<n; i++) for (a=0; a<3; a++) ra[i][a] = 0.0;

    /* Make a linked-cell list, lscl------------------------------------*/

    for (a=0; a<3; a++) lc2[a] = lc[a]+2;
    lcyz2 = lc2[1]*lc2[2];
    lcxyz2 = lc2[0]*lcyz2;

    /* Reset the headers, head */
    for (c=0; c<lcxyz2; c++) head[c] = EMPTY;

    /* Scan atoms to construct headers, head, & linked lists, lscl */

    for (i=0; i<n+nb; i++) {
        for (a=0; a<3; a++) mc[a] = (r[i][a]+rc[a])/rc[a];

        /* Translate the vector cell index, mc, to a scalar cell index */
        c = mc[0]*lcyz2+mc[1]*lc2[2]+mc[2];

        /* Link to the previous occupant (or EMPTY if you're the 1st) */
        lscl[i] = head[c];

        /* The last one goes to the header */
        head[c] = i;
    } /* Endfor atom i */

    /* Calculate pair interaction---------------------------------------*/

    rrCut = RCUT*RCUT;
    dr2 = rrCut/NTMAX;
    dr2i = 1.0/dr2;

    /* Scan inner cells */
    for (mc[0]=1; mc[0]<=lc[0]; (mc[0])++)
        for (mc[1]=1; mc[1]<=lc[1]; (mc[1])++)
            for (mc[2]=1; mc[2]<=lc[2]; (mc[2])++) {

                /* Calculate a scalar cell index */
                c = mc[0]*lcyz2+mc[1]*lc2[2]+mc[2];
                /* Skip this cell if empty */
                if (head[c] == EMPTY) continue;

                /* Scan the neighbor cells (including itself) of cell c */
                for (mc1[0]=mc[0]-1; mc1[0]<=mc[0]+1; (mc1[0])++)
                    for (mc1[1]=mc[1]-1; mc1[1]<=mc[1]+1; (mc1[1])++)
                        for (mc1[2]=mc[2]-1; mc1[2]<=mc[2]+1; (mc1[2])++) {

                            /* Calculate the scalar cell index of the neighbor cell */
                            c1 = mc1[0]*lcyz2+mc1[1]*lc2[2]+mc1[2];
                            /* Skip this neighbor cell if empty */
                            if (head[c1] == EMPTY) continue;

                            /* Scan atom i in cell c */
                            i = head[c];
                            while (i != EMPTY) {
                                ia = atype[i];

                                /* Scan atom j in cell c1 */
                                j = head[c1];
                                while (j != EMPTY) {
                                    ib = atype[j];

                                    /* No calculation with itself */
                                    if (j != i) {
                                        /* Logical flag: intra(true)- or inter(false)-pair atom */
                                        bintra = (j < n);

                                        /* Pair vector dr = r[i] - r[j] */
                                        for (rr=0.0, a=0; a<3; a++) {
                                            dr[a] = r[i][a]-r[j][a];
                                            rr += dr[a]*dr[a];
                                        }

                                        /* Calculate potential & forces for intranode pairs (i < j)
                                           & all the internode pairs if rij < RCUT; note that for
                                           any copied atom, i < j */
                                        if (i<j && rr<rrCut) {
                                            rdiv = rr * dr2i;
                                            idiv = (int) rdiv;
                                            rdiv = rdiv - idiv;
                                            vr = vtable[ia][ib][idiv] - rdiv*(vtable[ia][ib][idiv] - vtable[ia][ib][idiv+1]);
                                            fr = ftable[ia][ib][idiv] - rdiv*(ftable[ia][ib][idiv] - ftable[ia][ib][idiv+1]);

                                            lpe += 0.5*vr;
                                            ra[i][0] += fr*dr[0];
                                            ra[i][1] += fr*dr[1];
                                            ra[i][2] += fr*dr[2];
                                            if(bintra) {
                                                lpe += 0.5*vr;
                                                ra[j][0] -= fr*dr[0];
                                                ra[j][1] -= fr*dr[1];
                                                ra[j][2] -= fr*dr[2];
                                            }
                                        }

                                    } /* Endif not self */

                                    j = lscl[j];
                                } /* Endwhile j not empty */

                                i = lscl[i];
                            } /* Endwhile i not empty */

                        } /* Endfor neighbor cells, c1 */

            } /* Endfor central cell, c */

    /* save local potential energy */
    potEnergy0=lpe;
}

/*--------------------------------------------------------------------*/
void eval_props() {
    /*----------------------------------------------------------------------
    Evaluates physical properties: kinetic, potential & total energies.
    ----------------------------------------------------------------------*/
    double vv,lke,msd,vac;
    int i,a;
    double fs;

    /* Global potential energy */
    MPI_Allreduce(&potEnergy0,&potEnergy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    /* Total kinetic energy */
    for (lke=0.0, i=0; i<n; i++) {
        for (vv=0.0, a=0; a<3; a++) vv += rv[i][a]*rv[i][a];
        lke += vv;
    }
    lke *= 0.5;
    MPI_Allreduce(&lke,&kinEnergy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    /* Energy per atom */
    kinEnergy /= nglob;
    potEnergy /= nglob;
    totEnergy = kinEnergy + potEnergy;
    temperature = kinEnergy*2.0/3.0;

    /* Temperature control */
    if(mdmode==4) {
        fs = sqrt(InitTemp/temperature);
        for (lke=0.0, i=0; i<n; i++)
            for (vv=0.0, a=0; a<3; a++) rv[i][a]=fs*rv[i][a];
    }

    msd=compute_msd();
    vac=compute_vac();

    /* Print the computed properties */
    if (sid == 0) printf("%9.2f %9.6f %9.6f %9.6f  %9.6lf %9.6lf\n",
                             stepCount*DeltaT,temperature,potEnergy,totEnergy,msd,vac);
}

/*----------------------------------------------------------------------
Bit condition functions:

1. bbd(ri,ku) is .true. if coordinate ri[3] is in the boundary to
     neighbor ku.
2. bmv(ri,ku) is .true. if an atom with coordinate ri[3] has moved out
     to neighbor ku.
----------------------------------------------------------------------*/
inline int bbd(double* ri, int ku) {
    int kd,kdd;
    kd = ku/2; /* x(0)|y(1)|z(2) direction */
    kdd = ku%2; /* Lower(0)|higher(1) direction */
    if (kdd == 0)
        return ri[kd] < RCUT;
    else
        return al[kd]-RCUT < ri[kd];
}

inline int bmv(double* ri, int ku) {
    int kd,kdd;
    kd = ku/2; /* x(0)|y(1)|z(2) direction */
    kdd = ku%2; /* Lower(0)|higher(1) direction */
    if (kdd == 0)
        return ri[kd] < 0.0;
    else
        return al[kd] < ri[kd];
}

/*--------------------------------------------------------------------*/
void atom_copy() {
    /*----------------------------------------------------------------------
    Exchanges boundary-atom coordinates among neighbor nodes:  Makes
    boundary-atom list, LSB, then sends & receives boundary atoms.
    ----------------------------------------------------------------------*/
    int kd,kdd,i,ku,nsd,nrc,a;
    int inode_s,inode_r;
    int nbnew = 0; /* # of "received" boundary atoms */
    double com1;

    /* Main loop over x, y & z directions starts--------------------------*/

    for (kd=0; kd<3; kd++) {

        /* Make a boundary-atom list, LSB---------------------------------*/

        /* Reset the # of to-be-copied atoms for lower&higher directions */
        for (kdd=0; kdd<2; kdd++) lsb[2*kd+kdd][0] = 0;

        /* Scan all the residents & copies to identify boundary atoms */
        for (i=0; i<n+nbnew; i++) {
            for (kdd=0; kdd<2; kdd++) {
                ku = 2*kd+kdd; /* Neighbor ID */
                /* Add an atom to the boundary-atom list, LSB, for neighbor ku
                   according to bit-condition function, bbd */
                if (bbd(r[i],ku)) lsb[ku][++(lsb[ku][0])] = i;
            }
        }

        /* Message passing------------------------------------------------*/

        com1=MPI_Wtime(); /* To calculate the communication time */

        /* Loop over the lower & higher directions */
        for (kdd=0; kdd<2; kdd++) {

            inode_s = nns[ku=2*kd+kdd]; /* Neighbor node ID */
            inode_r = nnr[ku=2*kd+kdd]; /* Neighbor node ID */

            /* Send & receive the # of boundary atoms-----------------------*/

            nsd = lsb[ku][0]; /* # of atoms to be sent */

            /* Even node: send & recv */
            if (myparity[kd] == 0) {
                MPI_Send(&nsd,1,MPI_INT,inode_s,10,MPI_COMM_WORLD);
                MPI_Recv(&nrc,1,MPI_INT,inode_r,20,MPI_COMM_WORLD,&status);
            }
            /* Odd node: recv & send */
            else if (myparity[kd] == 1) {
                MPI_Recv(&nrc,1,MPI_INT,inode_r,10,MPI_COMM_WORLD,&status);
                MPI_Send(&nsd,1,MPI_INT,inode_s,20,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else
                nrc = nsd;
            /* Now nrc is the # of atoms to be received */

            /* Send & receive information on boundary atoms-----------------*/

            /* Message buffering */
            for (i=1; i<=nsd; i++) {
                dbuf[4*(i-1)] = (double) atype[lsb[ku][i]];
                for (a=0; a<3; a++) /* Shift the coordinate origin */
                    dbuf[4*(i-1)+a+1] = r[lsb[ku][i]][a]-sv[ku][a];
            }

            /* Even node: send & recv */
            if (myparity[kd] == 0) {
                MPI_Send(dbuf,4*nsd,MPI_DOUBLE,inode_s,30,MPI_COMM_WORLD);
                MPI_Recv(dbufr,4*nrc,MPI_DOUBLE,inode_r,40,MPI_COMM_WORLD,&status);
            }
            /* Odd node: recv & send */
            else if (myparity[kd] == 1) {
                MPI_Recv(dbufr,4*nrc,MPI_DOUBLE,inode_r,30,MPI_COMM_WORLD,&status);
                MPI_Send(dbuf,4*nsd,MPI_DOUBLE,inode_s,40,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else {
                for (i=0; i<4*nrc; i++) dbufr[i] = dbuf[i];
            }

            /* Message storing */
            for (i=0; i<nrc; i++) {
                atype[n+nbnew+i] = (int) dbufr[4*i];
                for (a=0; a<3; a++) r[n+nbnew+i][a] = dbufr[4*i+a+1];
            }

            /* Increment the # of received boundary atoms */
            nbnew = nbnew+nrc;
        } /* Endfor lower & higher directions, kdd */

        comt += MPI_Wtime()-com1; /* Update communication time, COMT */

    } /* Endfor x, y & z directions, kd */

    /* Main loop over x, y & z directions ends--------------------------*/

    /* Update the # of received boundary atoms */
    nb = nbnew;
}

/*--------------------------------------------------------------------*/
void atom_move() {
    /*----------------------------------------------------------------------
    Sends moved-out atoms to neighbor nodes and receives moved-in atoms
    from neighbor nodes.  Called with n, r[0:n-1] & rv[0:n-1], atom_move
    returns a new n' together with r[0:n'-1] & rv[0:n'-1].
    ----------------------------------------------------------------------*/

    /* Local variables------------------------------------------------------

    mvque[6][NBMAX]: mvque[ku][0] is the # of to-be-moved atoms to neighbor
      ku; MVQUE[ku][k>0] is the atom ID, used in r, of the k-th atom to be
      moved.
    ----------------------------------------------------------------------*/
    int mvque[6][NBMAX];
    int newim = 0; /* # of new immigrants */
    int ku,kd,i,kdd,kul,kuh,ipt,a,nsd,nrc;
    int inode_s,inode_r;
    double com1;

    /* Reset the # of to-be-moved atoms, MVQUE[][0] */
    for (ku=0; ku<6; ku++) mvque[ku][0] = 0;

    /* Main loop over x, y & z directions starts------------------------*/

    for (kd=0; kd<3; kd++) {

        /* Make a moved-atom list, mvque----------------------------------*/

        /* Scan all the residents & immigrants to list moved-out atoms */
        for (i=0; i<n+newim; i++) {
            kul = 2*kd  ; /* Neighbor ID */
            kuh = 2*kd+1;
            /* Register a to-be-copied atom in mvque[kul|kuh][] */
            if (r[i][0] > MOVED_OUT) { /* Don't scan moved-out atoms */
                /* Move to the lower direction */
                if (bmv(r[i],kul)) mvque[kul][++(mvque[kul][0])] = i;
                /* Move to the higher direction */
                else if (bmv(r[i],kuh)) mvque[kuh][++(mvque[kuh][0])] = i;
            }
        }

        /* Message passing with neighbor nodes----------------------------*/

        com1 = MPI_Wtime();

        /* Loop over the lower & higher directions------------------------*/

        for (kdd=0; kdd<2; kdd++) {

            inode_s = nns[ku=2*kd+kdd]; /* Neighbor node ID */
            inode_r = nnr[ku=2*kd+kdd]; /* Neighbor node ID */

            /* Send atom-number information---------------------------------*/

            nsd = mvque[ku][0]; /* # of atoms to-be-sent */

            /* Even node: send & recv */
            if (myparity[kd] == 0) {
                MPI_Send(&nsd,1,MPI_INT,inode_s,110,MPI_COMM_WORLD);
                MPI_Recv(&nrc,1,MPI_INT,inode_r,120,MPI_COMM_WORLD,&status);
            }
            /* Odd node: recv & send */
            else if (myparity[kd] == 1) {
                MPI_Recv(&nrc,1,MPI_INT,inode_r,110,MPI_COMM_WORLD,&status);
                MPI_Send(&nsd,1,MPI_INT,inode_s,120,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else
                nrc = nsd;
            /* Now nrc is the # of atoms to be received */

            /* Send & receive information on boundary atoms-----------------*/

            /* Message buffering */
            for (i=1; i<=nsd; i++) {
                dbuf[13*(i-1)] = (double)atype[mvque[ku][i]];
                for (a=0; a<3; a++) {
                    /* Shift the coordinate origin */
                    dbuf[13*(i-1)+1+a] = r [mvque[ku][i]][a]-sv[ku][a];
                    dbuf[13*(i-1)+4+a] = rv[mvque[ku][i]][a];
                    dbuf[13*(i-1)+7+a] = r0[mvque[ku][i]][a]-sv[ku][a];
                    dbuf[13*(i-1)+10+a] = v0[mvque[ku][i]][a];
                    r[mvque[ku][i]][0] = MOVED_OUT; /* Mark the moved-out atom */
                }
            }

            /* Even node: send & recv, if not empty */
            if (myparity[kd] == 0) {
                MPI_Send(dbuf,13*nsd,MPI_DOUBLE,inode_s,130,MPI_COMM_WORLD);
                MPI_Recv(dbufr,13*nrc,MPI_DOUBLE,inode_r,140,MPI_COMM_WORLD,&status);
            }
            /* Odd node: recv & send, if not empty */
            else if (myparity[kd] == 1) {
                MPI_Recv(dbufr,13*nrc,MPI_DOUBLE,inode_r,130,MPI_COMM_WORLD,&status);
                MPI_Send(dbuf,13*nsd,MPI_DOUBLE,inode_s,140,MPI_COMM_WORLD);
            }
            /* Single layer: Exchange information with myself */
            else
                for (i=0; i<13*nrc; i++) dbufr[i] = dbuf[i];

            /* Message storing */
            for (i=0; i<nrc; i++) {
                atype [n+newim+i] = (int) dbufr[13*i];
                for (a=0; a<3; a++) {
                    r [n+newim+i][a] = dbufr[13*i+1+a];
                    rv[n+newim+i][a] = dbufr[13*i+4+a];
                    r0[n+newim+i][a] = dbufr[13*i+7+a];
                    v0[n+newim+i][a] = dbufr[13*i+10+a];
                }
            }

            /* Increment the # of new immigrants */
            newim = newim+nrc;

        } /* Endfor lower & higher directions, kdd */

        comt=comt+MPI_Wtime()-com1;

    } /* Endfor x, y & z directions, kd */

    /* Main loop over x, y & z directions ends--------------------------*/

    /* Compress resident arrays including new immigrants */

    ipt = 0;
    for (i=0; i<n+newim; i++) {
        if (r[i][0] > MOVED_OUT) {
            atype [ipt] = atype [i];
            for (a=0; a<3; a++) {
                r [ipt][a] = r [i][a];
                rv[ipt][a] = rv[i][a];
                r0[ipt][a] = r0[i][a];
                v0[ipt][a] = v0[i][a];
            }
            ++ipt;
        }
    }

    /* Update the compressed # of resident atoms */
    n = ipt;
}

/*---------------------------------------------------------------------*/
void make_tables() {
    /*----------------------------------------------------------------------
     * make potential & force tables. the tables have values upto RCUT/NTMAX,
     * the last point in the two tables to make sure r=RCUT being treated propery.
     * ----------------------------------------------------------------------*/
    int i,ia,ib;
    double ri, ri2, ri6, ri12;
    double currentr, currentr2, voffset, foffset;
    double dr2,rcab;
    double aa[NC][NC],ss[NC][NC];

// V(r)[i][j] = aa[i][j]*( (ss[i][j]/r)^12 - (ss[i][j]/r)^6 )

    for (ia=0; ia<NC; ia++) {
        for (ib=0; ib<NC; ib++) {
            aa[ia][ib]=0.5*(Ac[ia]+Ac[ib]);   // Energy coefficient
            ss[ia][ib]=0.5*(Sg[ia]+Sg[ib]);   // Sigma
        }
    }

//aa&ss values taken from Srikanth et al., Nature, p554, vol. 398 (1998)
    aa[0][0]=1.0;
    aa[0][1]=1.5;
    aa[1][0]=1.5;
    aa[1][1]=0.5;

    ss[0][0]=1.0;
    ss[0][1]=0.8;
    ss[1][0]=0.8;
    ss[1][1]=0.88;

// For each sample point, calculate the value of potential table elements
// First creat the tables without the energy cofficient, then multiply aa[][].
//dr2=RCUT*RCUT/NTMAX;
    for (ia=0; ia<NC; ia++) {
        for (ib=0; ib<NC; ib++) {

            rcab=RCUT*ss[ia][ib];
            dr2=rcab*rcab/NTMAX;
            // Compute potential energy at RCUT
            ri = ss[ia][ib] / rcab;
            ri2 = ri * ri;
            ri6 = ri2 * ri2 * ri2;
            ri12 = ri6 * ri6;
            voffset = 4.0 * (ri12 - ri6);
            foffset = (48.0 * ri12 - 24.0 * ri6) * ri;

            for(i = 0; i<NTMAX; i++) {
                currentr2 = (double)i * dr2;
                currentr = sqrt(currentr2);
                ri = ss[ia][ib] / currentr;
                ri2 = ri * ri;
                ri6 = ri2 * ri2 * ri2;
                ri12 = ri6 * ri6;
                vtable[ia][ib][i] = 4.0 * (ri12 - ri6) - voffset + foffset * (currentr - rcab);
                ftable[ia][ib][i] = (48.0 * ri12 - 24.0 * ri6) * ri2 - foffset * ri;
                vtable[ia][ib][i] *= aa[ia][ib];
                ftable[ia][ib][i] *= aa[ia][ib];
                //if(i%20==0 && sid==0) printf("%d %d %f %f %f %f\n", ia, ib, currentr, currentr2, vtable[ia][ib][i],ftable[ia][ib][i]);
            }

            //Take care of the last point
            vtable[ia][ib][NTMAX-1] = 0.0;
            ftable[ia][ib][NTMAX-1] = 0.0;
            vtable[ia][ib][NTMAX] = 0.0;
            ftable[ia][ib][NTMAX] = 0.0;
        }
    }

}

/*---------------------------------------------------------------------*/
void gather_coordinates() {
    /*---------------------------------------------------------------------
    ----------------------------------------------------------------------*/
    int i,a,nex;
    double *sbuf,*rbuf;

    nex=0;
    MPI_Exscan(&n,&nex,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    //printf("sid=%d n=%d nex=%d\n",sid,n,nex);

    // allocate and clean up arrays
    sbuf = (double *)malloc(sizeof(double)*NMAX*3);
    for(i=0; i<NMAX; i++)
        for(a=0; a<3; a++)
            sbuf[i*3+a]=0.0;

    // initialize dbuf with global coord.
    for(i=0; i<n; i++)
        for(a=0; a<3; a++)
            sbuf[(i+nex)*3+a]=r[i][a]+ol[a];

    rbuf = (double *)malloc(sizeof(double)*NMAX*3);
    MPI_Allreduce(sbuf,rbuf,NMAX*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    // rbuf has all atom global coordinates by now.
    for(i=0; i<nglob; i++)
        for(a=0; a<3; a++)
            rg[i][a]=rbuf[i*3+a];

    free(sbuf);
    free(rbuf);
}
/*---------------------------------------------------------------------*/
void compute_gr(double* hist, int nbin, double grrcut, double grdr) {
    /*---------------------------------------------------------------------
    ----------------------------------------------------------------------*/
    int i,j,a;
    int idx;
    double *rbuf, *sbuf;
    double rdif,rr;
    FILE *fp;

    //nbin = grrcut/grdr;
    //printf("sid = %d, grdr = %lf, nbin=%d\n", sid, grdr, nbin);
    //hist = (int *)malloc(sizeof(int )*nbin);

    sbuf = (double *)malloc(sizeof(double)*nbin);
    for(i=0; i<nbin; i++) sbuf[i]=0.0;

    for(i=sid; i<nglob-1; i+=nprocs) {
        for(j=i+1; j<nglob; j++) {

            if(i==j) continue;

            for(rr=0.0, a=0; a<3; a++)
            {
                rdif = rg[j][a]-rg[i][a];
                if(rdif>=0.5*gl[a]) rdif-=gl[a];
                if(rdif<-0.5*gl[a]) rdif+=gl[a];
                rr += rdif*rdif;
            }
            rr = sqrt(rr);
            if(rr<grrcut)
            {
                idx = rr/grdr;
                sbuf[idx]+=2.0;
            }
        }
    }

    rbuf = (double *)malloc(sizeof(double)*nbin);
    MPI_Allreduce(sbuf,rbuf,nbin,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    for(int i=0; i<nbin; i++)
        hist[i]+=rbuf[i];

    free(sbuf);
    free(rbuf);
}

/*---------------------------------------------------------------------*/
double compute_vac()
/*---------------------------------------------------------------------
 * Compute velocity auto-correlation
----------------------------------------------------------------------*/
{
    int i, a;
    double  vv,v0nrm,v1nrm,vacl,vac=0.0;

    vacl=0.0;
    for (i=0; i<n; i++) {
        vv=0.0;
        v0nrm=0.0;
        v1nrm=0.0;
        for (a=0; a<3; a++) {
            v0nrm += v0[i][a]*v0[i][a];
            v1nrm += rv[i][a]*rv[i][a];
            vv += v0[i][a]*rv[i][a];
        }
        vacl += vv/(sqrt(v0nrm)*sqrt(v1nrm));
    }

    MPI_Allreduce(&vacl,&vac,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    vac /= nglob;

    return vac;
}

/*---------------------------------------------------------------------*/
double compute_msd()
/*---------------------------------------------------------------------
 * Compute mean square displacement
----------------------------------------------------------------------*/
{
    int i, a;
    double rr,rr1,msd=0.0;

    rr=0.0;
    for (i=0; i<n; i++) {
        for (a=0; a<3; a++) {
            rr1 = r[i][a]-r0[i][a];
            rr += rr1*rr1;
        }
    }

    MPI_Allreduce(&rr,&msd,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    msd /= nglob;

    return msd;
}

/*---------------------------------------------------------------------*/
void analysis_manager(int phase) {
    /*---------------------------------------------------------------------*/
    /*---------------------------------------------------------------------*/

    static int idx = 0, NSAMPLES = 50;
    static int *COUNT;
    static double *MSD, *VAC, *GR;
    static double grdr=0.1, grrcut=5;
    static int nbin, grcount;
    int i,j,a;
    double gr,rr,cr,denom,pi=3.14159265359;
    FILE *fp;

    enum {init,record,final};

    switch(phase)
    {
    case (init):
        COUNT = (int *) malloc(sizeof(int)*NSAMPLES);
        MSD = (double *) malloc(sizeof(double)*NSAMPLES);
        VAC = (double *) malloc(sizeof(double)*NSAMPLES);
        for(int i=0; i<NSAMPLES; i++)
        {
            COUNT[i]=0;
            MSD[i]=0.0;
            VAC[i]=0.0;
        }

        /* get number of bins and reset gr*/
        nbin=grrcut/grdr;
        GR = (double *) malloc(sizeof(double)*nbin);
        for(int i=0; i<nbin; i++) GR[i]=0.0;
        grcount = 0;

        /* Initialize the initial coordinate */
        for (j=0; j<n; j++)
            for(a=0; a<3; a++)
            {
                r0[j][a] = r[j][a];
                v0[j][a] = rv[j][a];
            }

        break;

    case (record):

        if(idx==NSAMPLES)
        {
            idx=0;

            /* Initialize the initial coordinate */
            for (j=0; j<n; j++)
            {
                for(a=0; a<3; a++)
                {
                    r0[j][a] = r[j][a];
                    v0[j][a] = rv[j][a];
                }
            }
        }

        grcount++;
        gather_coordinates();
        compute_gr(GR, nbin, grrcut, grdr);

        MSD[idx] += compute_msd();
        VAC[idx] += compute_vac();
        COUNT[idx]++;

        idx++;

        break;

    case (final):

        if(sid==0)
        {
            fp = fopen("msd.d","w");
            for(int i=0; i<NSAMPLES; i++)
            {
                if(COUNT[i]>0)
                    fprintf(fp,"%d %lf %d\n", i, MSD[i]/COUNT[i], COUNT[i]);
            }
            fclose(fp);
            free(MSD);

            fp = fopen("vac.d","w");
            for(int i=0; i<NSAMPLES; i++)
            {
                if(COUNT[i]>0)
                    fprintf(fp,"%d %lf %d\n", i, VAC[i]/COUNT[i], COUNT[i]);
            }
            fclose(fp);
            free(VAC);

            fp=fopen("gr.d","w");
            cr=0.0;
            for(int i=1; i<nbin; i++)
            {
                rr = i*grdr;
                denom = 4.0*pi*(rr*rr*grdr*Density)*nglob;
                gr = (double)GR[i]/denom/grcount;
                cr += gr;
                fprintf(fp,"%lf  %lf %lf\n", rr, gr, cr);
            }
            fclose(fp);
        }

        break;

    default:
        break;
    }
}

/*---------------------------------------------------------------------*/
void write_config(int nstep) {
    /*---------------------------------------------------------------------
     * This funcition serializes writing file process in order to
     * save all data in one file. A write token will be passed from
     * MPIrank==0 to nprocs-1. whoever holds the token can access file.
    ----------------------------------------------------------------------*/
    FILE *fp;
    char fname[19];
    uint32_t i,a;
    long int *iary, *nary;
    double rs[3];

    // allocate temporary arrays to store local number of atoms
    iary=malloc(sizeof(long int)*nprocs);
    nary=malloc(sizeof(long int)*nprocs);
    for(i=0; i<nprocs; i++) iary[i]=0;
    iary[sid]=n;

    // copy all local numbers of atoms into one array
    MPI_Allreduce(iary,nary,nprocs,MPI_LONG,MPI_SUM,MPI_COMM_WORLD);

    // if MPIrank != sid, receive write token from previous MPI task
    if(0!=sid) MPI_Recv(&i,1,MPI_INT,sid-1,sid-1,MPI_COMM_WORLD,&status);

    // Open File
    if(nstep==-1) {
        // nstep==1 save final configulation for next run
        if(sid==0) {
            fp = fopen("pmd.d","w");
        }
        else {
            fp = fopen("pmd.d","a");
        }
        // set shift vector zero, save local atom position
        for(a=0; a<3; a++) rs[a]=0;
    } else {
        // save check point data for visualization & analysis
        sprintf(fname,"data/pmd%06d.xyz",nstep/1000);
        if(sid==0) {
            fp = fopen(fname,"w");
        }
        else {
            fp = fopen(fname,"a");
        }
        // set shift vector domain dims, save global atom position
        for(a=0; a<3; a++) rs[a]=ol[a];
    }

    //--- the serialization part begins ---//

    if(sid==0) { // if MPIrank = 0, write header part
        fprintf(fp, "%d %d %f %f %f\n", nglob, (currCount+StepLimit)/1000, gl[0],gl[1],gl[2]);
        for(i=0; i<nprocs; i++) fprintf(fp, "%ld ",nary[i]);
        fprintf(fp, "\n");
    }

    // write atom type, position, velocity & close file
    for(a=0; a<n; a++) {
        fprintf(fp,"%d %lf %lf %lf %lf %lf %lf\n",
                atype[a], r[a][0]+rs[0],r[a][1]+rs[1],r[a][2]+rs[2],
                rv[a][0],rv[a][1],rv[a][2]);
    }
    fclose(fp);

    // if MPIrank != nprocs-1, send write token to next MPI task
    if(sid!=nprocs-1) MPI_Send(&i,1,MPI_INT,sid+1,sid,MPI_COMM_WORLD);

    //--- the serialization part over ---//

    // free temporary arrays
    free(iary);
    free(nary);
}
