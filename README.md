This repository contains supplementary files to the paper ["Maximal distance minimizers for a rectangle"](https://arxiv.org/abs/2106.00809).

[verifier folder](https://github.com/c0pymaster/minimizer/tree/master/verifier) contains the source code of the computer search program described in the paper. It can be compiled as follows:
```
g++ -o main -std=c++11 -pthread -DCXSC_FORCE_TLS -DCXSC_USE_TLS_PREC -O2 -Wall -Werror -Wextra -mfpmath=sse -msse2 -I<path to C-XSC>/cxsc/include -L<path to C-XSC>/cxsc/lib main.cpp -lcxsc
```
Note that C-XSC library should be preinstalled on your machine. You can find the library source code and installation instructions on http://www.xsc.de/. Please make sure to list flags `-O2 -DCXSC_FORCE_TLS -DCXSC_USE_TLS_PREC` for the installation script to use during the compilation of the library. Without the last two flags the library will not work properly in multi-threaded environment.

After compilation, the program can be run as follows:
```
./main <t>
```
where t is the number of threads you would like to use. If t=0, the program will try to guess the number of concurrent threads supported by the machine.

[Case4](https://github.com/c0pymaster/minimizer/blob/master/Case4.ipynb) contains the computations used in case 4 in Section 6 of the paper.
