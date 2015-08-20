SIMD oriented Fast Mersenne Twister sample
==========================================

sfmtdc
usage:
./sfmtdc [-s seed] [-v] [-c count] [-f outputfile] mexp

--verbose, -v        Verbose mode. Output parameters, calculation time, etc.
--file, -f filename  Parameters are outputted to this file. without this
                     option, parameters are outputted to standard output.
--count, -c count    Output count. The number of parameters to be outputted.
--seed, -s seed      seed of randomness.
mexp                 mersenne exponent.

calc_equidist
usage:
./calc_equidist [-v] mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,parity1,parity2,parity3,parity4

--verbose, -v        Verbose mode. Output detailed information.
other arguments are outputs of sfmtdc

test_linearity
./test_linearity
 mexp,pos1,sl1,sl2,sr1,sr2,msk1,msk2,msk3,msk4,parity1,parity2,parity3,parity4
arguments are outputs of sfmtdc
