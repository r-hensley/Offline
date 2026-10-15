# Experimental straw-hit local-reference subtraction

When a straw hit's closest approach moves to another material piece after fitting,
its cached information matrix can still describe the previous measurement
location. Subtracting that information from the newly selected local fitted state
need not produce a positive-definite remainder. Clipping a projected variance
does not correct this state/information mismatch.

This draft ports a previously developed local-reference diagnostic into the
ordinary straw-hit interface. It is a new source integration, not a validated
production correction. It is experimental, version-specific, and not an endorsed
reconstruction prescription or tutorial. No physics results or performance claims
are included here.

## Optional behavior

`KinematicLineFit` adds the module setting
`ExperimentalUnbiasedLocalReference`, with default `false`. Enabling it requires
the companion KinKal `AcceptedFitReference` hook and disables neither existing
selection cuts nor covariance inversion errors. Magnetic-domain transport in
either the initial or extension configuration is rejected for this option.
The option is exposed only by the KinematicLine module; other trajectory models
are outside this draft's intended scope.

The diagnostic follows this sequence:

1. After `updateWeight`, save the actual information matrix/vector, particle and
   sensor measurement times, and a weight-generation identifier. Invalidate any
   older posterior snapshot.
2. During accepted post-fit reference attachment, copy the new trajectory's local
   piece at that saved measurement time, before the ordinary biased closest-approach
   update can select a different piece. Both the trajectory range and the selected
   piece must contain the old time. This is a copy of the accepted post-fit state,
   not a reused pre-fit covariance.
3. For an eligible active hit, subtract the exact cached information from that
   copy and recompute both parameters and covariance. Start the new closest
   approach from the saved particle and sensor times.
4. For initialization, inactive hits, clones without a new accepted association,
   changed weights or unsupported old times, use the original calculation.

In information form the operation is

```text
W_without_hit = inverse(C_postfit_at_old_time) - W_cached_hit
w_without_hit = inverse(C_postfit_at_old_time) * a_postfit_at_old_time - w_cached_hit
C_without_hit = inverse(W_without_hit)
a_without_hit = C_without_hit * w_without_hit
```

The vector and matrix come from the same selected state. No absolute value,
variance floor, eigenvalue clipping, or acceptance bypass is introduced.

## Lifecycle and diagnostic counts

An external thread-local cache leaves the `KKStrawHit` object layout and virtual
interface unchanged. Each module invocation owns a fresh cache scope; nested scopes
restore their predecessor, including on exceptions. Ordinary and copy construction
discard any association at the new object's address. Copying a hit never copies
the diagnostic association, even in versions where its base weight is copied.
Capture and use require the same weight generation, exact full vector/matrix
equality, an active captured state, and supported old measurement time.

When enabled, the module prints completed local-calculation, legacy-fallback,
out-of-piece, out-of-sensor and unusable-local counts. These counts describe
calculations, not independent tracks or rates of improvement. They do not check
positive definiteness or uncertainty calibration. An inversion that throws can
terminate a fit before a completed local calculation is counted.

The cache exists only inside that module invocation. Subsequent modules and generic
`Hit::unbiasedParameters`, residual, pull and chi-squared consumers retain their
existing behavior. The option therefore does not claim framework-wide consistency.

## Scientific and implementation limits

The recomputed closest approach can leave the copied trajectory piece. In that
case this code extrapolates that one piece; it does not propagate the hit-removed
state through intervening material or recompute material crossing locations. It
counts the condition and returns the diagnostic calculation, matching the original
local experiment. Positive projected variance is insufficient to establish a
positive-definite covariance, and positive covariance is insufficient to establish
correct physical uncertainty.

The accepted state and remaining linearizations still depend on the prior all-hit
fit and its selection history. This is not a complete nonlinear fit with the hit
always absent. Field-dependent parameter-frame transport, boundary derivatives,
unsupported times, full material propagation, and calibrated uncertainty remain
outside this local experiment.

The first-piece prior correction is a separate change. This draft neither includes
nor enables it. There is no code dependency on that correction, but the preceding
local experiment used it together with the updated closest-approach hint. A future
comparison must declare and match those baseline choices. The companion KinKal
context hook is a required compile-time dependency.

## Source review and validation status

The source mechanism remains present in KinKal main
`71cc7d928d27b04d3475ac47bea29709481f32f4` and Offline main
`8fdf0eadfc1758810cacd7da6cdedc6895ac4862`: cached hit information is still removed
from the current closest-approach reference. The newer Offline clone constructor
copies the base residual weight and its branch check uses the newly calculated
PCA. Those changes do not implement this old-location posterior association.
Its background ANN's clamp of the scalar distance variance also leaves the full
matrix subtraction unchanged. This source observation does not prove that any
particular event still fails under those versions.

The draft branch is based on the fork's Offline main
`4b764401cfd1707531a4f8654bb591b32f0977c8`; unrelated changes from later upstream main
are not bundled into this proposal. Porting onto a newer base must preserve its
clone, reference-selection and module changes.

This port has not been built or run. No tests or new reconstruction jobs were
executed for this source-only draft. Previous checks of the earlier isolated
implementation do not validate this integration. Before enabling it for a study,
a consistent rebuild, disabled-option parity, cache lifecycle checks, full matrix
checks, support checks and a matched comparison against a transported hit-removed
reference remain necessary. Production use and general reconstruction improvement
are not claimed.
