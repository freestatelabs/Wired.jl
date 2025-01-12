# C Kernel

`Wired.jl` includes an optional computational kernel written in C. It is meant to 
be compiled natively to the user's machine, thereby allowing additional optimizations to 
occur for further efficiency. 

The `Wired.jl` C kernel often allows for speedups of 2-4x compared to the pure Julia
implementation, but requires additional installation steps.

## Installation 

Installation is performed via the `installkernel()` function. This requires use of 
the GCC compiler, callable with `gcc` from the command line. Therefore, the program
must be run on Linux, MacOS, or WSL on Windows. 

Once the kernel is compiled (should take just a few seconds), the kernel can be 
switched via setting `Wired.kernel = "c"`. Future calls to the `bfield()` function
will now be directly to the C kernel. 