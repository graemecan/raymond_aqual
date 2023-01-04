RAyMOND is a patch to modify the N-body/hydrodynamics code RAMSES in order to be able to run simulations in MOND gravity. It includes both the fully non-linear AQUAL formulation of MOND, and the quasi-linear QUMOND formulation (in separate repositories as two separate patches). For details and examples of using the code, see these papers:

* [RAyMOND: an N-body and hydrodynamics code for MOND (Candlish, Smith & Fellhauer 2015)](10.1093/mnras/stu2158)
* [The velocity field in MOND cosmology (Candlish 2016)](10.1093/mnras/stw1130)
* [Consequences of the external field effect for MOND disc galaxies in galaxy clusters (Candlish et al. 2018)](10.1093/mnras/sty2228)

Including either patch automatically switches the gravitational solver to that of QUMOND/AQUAL. For both patches an additional MOND\_PARAMS section must be added to the namelist with the following parameters:

* imond\_a0: This is the MOND acceleration scale in m/s^2. This is converted to code units automatically.
* mond\_n: This is the exponent used in the MOND interpolation function which is hard-coded to be of the form x/(1+x^n)^(1/n) for AQUAL, or its inverse in the case of QUMOND (see Famaey & McGaugh 2012, eqs. 49 and 50). It is an easy matter to use a different mu (or nu) function: just change the appropriate line in the function "mu\_function" (or "nu\_function") at the end of poisson\_commons.f90. Note that it's probably numerically more efficient to remove mond\_n completely from the interpolation function if you are using the "simple" version where mond_n = 1.

There is one additional patch-specific parameter for AQUAL runs:

* maxstoredcells: This is a technical parameter used to set the size of an array that stores extra grid cell values required due to the enlarged numerical stencil of the non-linear solver. Setting this to about 20% of the value of ngridmax is usually more than enough.

This patch has been checked to **compile** with the [stable 19_10](https://bitbucket.org/rteyssie/ramses/branch/stable_19_10) branch of the RAMSES code. No tests have yet been undertaken to check runtime compatibility.
