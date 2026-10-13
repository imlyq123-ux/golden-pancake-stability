# Multiprecision spectral transitions

Each row: one transition refined by several pipelines (secant on the transition function; every
evaluation is a fresh Newton-corrected orbit plus a separate monodromy integration). Confirmed digits
count only Taylor-vs-GBS agreement (independent implementations); Taylor-vs-Taylor agreement is a
precision-consistency check. Brackets are empirically verified sign changes, not interval enclosures.

| transition | kind | pipelines | confirmed estimate (truncated) | confirmed digits (pair) | Taylor-Taylor agreement | bracket half-width | endpoint values | slope |
|---|---|---|---|---|---|---|---|---|
| T1 ((ABab)^3 segment, s = +2 crossing past mu_min) | tangent | taylor50, gbs80 | 0.475901444550864118752069986423460087896... | 39 (taylor50~gbs80) | - | 1.0e-30 | -7.825e-30, 7.825e-30 | 7.8253722 |
| K0 (upper edge of lower stable stretch) | krein | taylor50, gbs80 | 0.580094414288468996548897180964879640274... | 39 (taylor50~gbs80) | - | 1.0e-30 | 2.225e-29, -2.225e-29 | -22.254192 |
| K1 = mu_L (lower edge of main window) | krein | taylor50, taylor70, taylor90, gbs80 | 0.599043784208899851220728860462068419275413725661780425744347... | 60 (taylor90~gbs80) | 53 | 1.0e-30 | -7.487e-29, 7.487e-29 | 74.874061 |
| K2 = mu_R (upper edge of main window) | krein | taylor50, taylor70, taylor90, gbs80 | 0.617844268209063178727366710600417996702072689168837458676169... | 60 (taylor90~gbs80) | 58 | 1.0e-30 | 3.059e-28, -3.059e-28 | -305.88257 |
| K3 (lower edge of near-fold stable stretch) | krein | taylor50, gbs80 | 0.634197066930554588552605424098150004643... | 39 (taylor50~gbs80) | - | 1.0e-30 | -5.84e-27, 5.84e-27 | 5839.8126 |
