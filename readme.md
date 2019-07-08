Code Author : Adil Ansari

![Image of Visualization](https://github.com/aansari2/Distributed-Computing-on-MPI/blob/master/visualization.png)

The command 'make' should compile both the mpi and coarray
programs but intel mpi must be loaded and coarray must be 
set to core count as as follows:
export FOR_COARRAY_NUM_IMAGES=16

The Single Core Implementation gives:

    ./jacobi 4000 1e-4; #Run Single Core Implement
    d = 9.993779e-05        c = 1.195785e-01

The MPI results are:

    mpirun -np 16 ./jacobimpi 4000 1e-4; #Run MPI Program
    d = 9.993779e-05        c = 1.195785e-01

The Coarray Results are:

    ./jacobicoarray 4000 1e-4; #Run Coarray Program
     d =   9.993779345127130E-005   c =   0.119578496530036

As we can see all results match up.

The coarray is significantly slower than mpi which is slower
than naive implementation because multiple cores aren't called
and even then we can't be certain if comm overheads will be 
small.
