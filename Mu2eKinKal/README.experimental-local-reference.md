# Experimental local-reference subtraction

A hit's closest approach can move to another material piece after fitting while
its cached weight still describes the old measurement location. Subtract that
weight from the accepted fitted piece at the old time, recomputing both mean and
covariance before finding closest approach.

This is automatic for `KinematicLine` hits whose accepted fit has magnetic-domain
correction disabled and no retained magnetic domains. It requires
[KinKal PR #2](https://github.com/r-hensley/KinKal/pull/2). Other trajectory types,
magnetic-domain fits and hits without a valid association use the legacy calculation.
There is no module switch, scope cache or diagnostic counter log.

Each hit owns its saved weight, particle/sensor times and copied posterior.
Every weight/reference update invalidates the posterior; only a supported accepted
fit can replace it. Use requires an active hit, supported old time and exact weight
matrix/vector equality. Clones start without an association. Unlike the former
module-scoped option, the association remains available to later callers until
invalidated. Consumers must be rebuilt together because the hit layout changed.

The new closest approach may leave the copied piece, extrapolating it without
transporting removed information through material. This remains an experimental
local calculation, not a nonlinear hit-free refit or calibrated uncertainty.
Generic hit-removed interfaces and the first-piece prior are unchanged; inversion
failures propagate. Compiler/fixture checks do not establish reconstruction
performance or replace a full linked build and matched event validation.
