# Background ANN single-hit regression

**Personal archive notice:** This branch is an unofficial personal archive of
version-specific development code. It is unmaintained as a distribution, is
not an endorsed Mu2e tutorial or result, and may have known limitations.

The background updater supplies one six-feature hit, while the former
generated evaluator requested a batch of 32 hits. This branch changes the
batch-dependent buffers and operations to one hit, preserving coefficient
loading, feature order, dense-layer dimensions, and activation arithmetic.

`TrainBkgSingleHit_test.cc` checks the input contract used by `BkgANNSHU`:
six readable floats must produce one score. It places those floats immediately
before an inaccessible memory page, so an incorrect batch dimension also fails
when the BLAS library was built without address sanitizers.

The test creates and removes its own temporary synthetic weight file. Its
coefficients exercise all six features, the three ReLU layers, and the sigmoid
against a known affine logit over repeated calls. It uses no trained weights,
event files, or reconstruction jobs.

In an existing Offline environment with ROOT SOFIE and OpenBLAS configured,
run from the `Offline` checkout root:

```sh
c++ -std=c++17 -O2 $(root-config --cflags) -I.. \
  Mu2eKinKal/test/TrainBkgSingleHit_test.cc \
  $(root-config --libs) -lROOTTMVASofie -lopenblas \
  -o /tmp/TrainBkgSingleHit_test
(ulimit -c 0; /tmp/TrainBkgSingleHit_test)
```

This standalone Linux test uses `mmap`/`mprotect`. With the former batch-32
header it should fail by an invalid read or by returning the wrong number of
scores; with the single-hit header it must exit zero. It checks memory bounds
and synthetic inference arithmetic, not reconstruction performance.

On 2026-10-15, both versions compiled with GCC 13.3.0, ROOT 6.32.06 and
OpenBLAS 0.3.28 in an existing p103 environment. The single-hit header passed;
the former batch-32 header failed with a segmentation violation in BLAS input
copying. These standalone checks preceded the instruction to stop further tests.
The newly selected upstream p108 environment has not been tested here.
