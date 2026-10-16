# Experimental local-reference subtraction

A hit's closest approach can move to another material piece after fitting while
its cached weight still describes the old measurement location. This option
subtracts that weight from the accepted fitted piece at the old time, recomputing
both the mean and covariance before finding closest approach.

`KinematicLineFit.ModuleSettings.ExperimentalUnbiasedLocalReference` defaults to
`false`. It requires [KinKal PR #2](https://github.com/r-hensley/KinKal/pull/2).
Magnetic-domain correction in either fit phase is rejected when enabled.

Each hit owns its saved weight, particle/sensor times and copied posterior.
Updating the weight or reference invalidates the posterior; only an accepted
fit can replace it. Use requires an active hit, supported old time, exact weight
matrix/vector equality, and the same module-call scope. Clones start without an
association. Missing associations retain the existing calculation. The module
scope prevents use by later consumers; there is no external cache or counter log.

The recomputed closest approach may leave the copied piece. That extrapolates
one piece without transporting the removed information through material. This
is an experimental local calculation, not a nonlinear hit-free refit or calibrated
uncertainty. It changes neither generic hit-removed parameters/residuals nor the
first-piece prior, and does not clip inversion failures or bypass selection.

This branch targets the personal fork's main; it is not an endorsed reconstruction
prescription. The changed hit layout requires a consistent rebuild of consumers.
Validation of this revision is recorded in the PR description; compile and
software checks do not establish reconstruction performance.
