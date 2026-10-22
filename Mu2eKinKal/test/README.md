# Single-hit background ANN regression

`BkgANNSHU` passes exactly six floats to `TrainBkg::Session::infer` and consumes
one score. `TrainBkgSingleHit_test.cc` places those six floats immediately before
an inaccessible memory page, detecting oversized BLAS input reads even when BLAS
is not instrumented with sanitizers.

The test writes temporary synthetic coefficients, loads them through the normal
session constructor, and removes the file before inference. All six features,
all four bias tensors, the three ReLUs, and the sigmoid are exercised over 20
calls on one session, including positive/negative impulses, activation boundaries,
saturated scores, and a return to zero input. Each score must agree within 1e-6
with the closed-form synthetic model; inputs must remain unchanged.

In a configured Linux environment with ROOT SOFIE and OpenBLAS, from the Offline
checkout root (whose parent contains the `Offline` directory):

```sh
c++ -std=c++17 -O2 -g $(root-config --cflags) -I.. \
  Mu2eKinKal/test/TrainBkgSingleHit_test.cc \
  $(root-config --libs) -lROOTTMVASofie -lopenblas \
  -o /tmp/TrainBkgSingleHit_test
(ulimit -c 0; OPENBLAS_NUM_THREADS=1 /tmp/TrainBkgSingleHit_test)
```

Compile this same test with the original header to reproduce the invalid read;
with the corrected header it must return one correct score per call and exit
zero. The test uses no trained weights or detector data. It establishes the
single-hit memory contract and bounded inference arithmetic, not reconstruction
performance or calibration. It is standalone and not registered in the general
Offline test suite.
