# Notebooks

Two runnable notebooks accompany this book, under `notebooks/`:

- **[quickstart.ipynb](notebooks/quickstart.ipynb)**: implements and validates a small
  finite-difference solver for the governing equation in
  [Solution methodology](methodology.md), checks it against the closed-form `sinf_head_perf` and
  `sinf_head_leak` solutions, then demonstrates a recharge-driven finite aquifer
  (`fin_rech_leak`).
- **[applied-examples.ipynb](notebooks/applied-examples.ipynb)**: runs three scenarios from
  [Applied examples](examples.md) end to end with synthetic data and plots: construction
  dewatering vs. a permit limit, snowmelt timing for a trout stream's baseflow, and a stormwater
  basin's 72-hour drawdown requirement.

**Why these don't call `flowsim.run_model()` directly:** the `flowsim` package currently on PyPI
fails to import in this environment (see [Installation](installation.md)); it installs as
version `0.0.0`, its own metadata describes an unrelated Edcrop script, and it's missing a
required module file. Until that's fixed upstream, both notebooks solve the same governing
equation and boundary conditions numerically (implicit finite differences) instead, and validate
the result against Flowsim's own closed-form analytical solutions before trusting it further. If
you have a working local build of the real package, the `flowsim.yaml` snippets throughout
[Inputs](inputs.md) and [Applied examples](examples.md) are what you'd feed to
`flowsim.run_model()` instead; the notebooks' `simulate_1d()` function is a drop-in stand-in for
that call.

The rest of the scenarios in [Applied examples](examples.md) (mine dewatering, levee
underseepage, riverbank filtration, etc.) follow the same pattern used in
`applied-examples.ipynb`: build a boundary-head or recharge time series, call `simulate_1d`, and
read the answer off the output. Contributions extending the notebook to more of them are welcome
at <https://github.com/SteenChr/flowsim>.
