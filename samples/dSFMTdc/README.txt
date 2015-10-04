double precision SIMD oriented Fast Mersenne Twister sample
===========================================================

dSFMTdc
usage:
./dSFMTdc [-s seed] [-v] [-c count] [-f outputfile] [-x [fixedSL]] mexp

--verbose, -v        Verbose mode. Output parameters, calculation time, etc.
--file, -f filename  Parameters are outputted to this file. without this
                     option, parameters are outputted to standard output.
--count, -c count    Output count. The number of parameters to be outputted.
--seed, -s seed      seed of randomness.
--fixed, -x fixedSL  fix the parameter sl1 to given value.
mexp                 mersenne exponent.

calc_equidist
usage:
./calc_equidist [-v] mexp,pos1,sl1,sr1,msk1,msk2,fixed1,fixed2,parity1,parity2

--verbose, -v        Verbose mode. Output detailed information.
other arguments are outputs of dSFMTdc

test_linearity
./test_linearity
mexp,pos1,sl1,sr1,msk1,msk2,fixed1,fixed2,parity1,parity2
arguments are outputs of dSFMTdc
