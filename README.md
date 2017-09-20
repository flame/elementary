# Elementary

Elementary is a pedagogical sandbox used for teaching collective communication and distributed memory parallel implementations of matrix multiplication, covering the examples in [a systematic journey](http://epubs.siam.org/doi/abs/10.1137/140993478). It is based on the early version of [Elemental](http://libelemental.org), a high-level high-performance framework for distributed-memory dense linear algebra that is inspired by [PLAPACK](http://www.cs.utexas.edu/users/plapack/new/using.html) and [FLAME](https://www.cs.utexas.edu/~flame/web/libFLAME.html).

The key idea behind the project is that an appropriate level of abstraction 
allows one to focus on *algorithms* rather than getting tied down in 
implementation details.

# Build

There are two main build modes, *debug* and *release*. The former maintains a call stack and performs judicious error-checking. Thus debug-mode should be used when testing a new algorithm. You can build it with
```
      make debug
```

If you would also like to build the debug test drivers for the distributed 
BLAS routines, run
```
      make test-debug
```

Similarly, you can build the baremetal version of the library by running 
```
      make release
```
The command
```
      make test-release
```
perform the equivalent action as in the debug case.

If neither debug nor release modes are specified, both are built. The relevant
commands are
```
      make
      make test
```

# Example Driver

Check test folder.

