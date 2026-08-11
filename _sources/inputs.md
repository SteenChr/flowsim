# Inputs

## Overview

To use Flowsim, all input must be in consistent units. Correspondingly, all output from Flowsim
(see [Outputs](outputs.md)) is presented in these same consistent units. The simulation is
carried out at discrete time steps of constant duration, defined by the time step of the
boundary-condition input.

All input must be provided in a required text file named `flowsim.yaml`. This file follows the
[YAML](https://yaml.org/) format, a human-readable data serialization standard. It must be
located in the working directory when Flowsim is executed.

---

## Structure of `flowsim.yaml`

The `flowsim.yaml` file may include:

- Blank lines
- Comments (start with `#`)
- Keys and key-value pairs
- Value continuations across multiple lines

### Value types

A value in the file can be:

- Text string (e.g., `S10_WatLvL`)
- Integer (e.g., `23`)
- Float (e.g., `7.831` or `3.52e-7`)
- Date (`%Y-%m-%d`, e.g., `1998-08-31`)
- Boolean (`True` / `False`)
- List (e.g., `[3.2, 7.32, 7.5e-3]`)
- Dictionary (e.g., `{key_1: 1.3, key_2: 2.7}`)

### Syntax rules

- Lists begin with `[` and end with `]`, entries separated by commas
- Dictionaries begin with `{` and end with `}`, key-value pairs separated by commas
- Use **spaces** for indentation (not tabs!) — using tabs causes an error and program failure
- Indentation defines blocks and sub-blocks of information; use identical indentation for
  sibling blocks (e.g. all aquifer sub-blocks), or file loading will fail

---

## Input to make a simulation

Only one key-value pair is mandatory:

```yaml
simulation_periods:
  - {begin: "2010-01-01", end: "2015-12-31"}
  - {begin: "2015-09-03", end: "2016-11-20", early_begin: "2005-09-03"}
```

Each simulation period requires `begin` and `end` dates. A period may optionally include
`early_begin`, the start of a warm-up period that precedes `begin`; results are only written for
the `begin`–`end` sub-period, but the warm-up lets the simulation "spin up" before results are
recorded.

### Optional setup keys

```yaml
dtformat: "%Y-%m-%d"   # default date format for the keys above
response: "head"        # default is "flux"
x: [0.0, 58.03, 122.89, 247.76, 375.01]   # locations to report head/flux at (default [0.0])
```

### Aquifer information (mandatory)

```yaml
aquifers:
  Upp:
    func:
      fin_rech_leak: "recharge"
      sinf_head_leak: "stream_stage"
    T: 259.2
    S: 0.20
    C: 22.06413
    L: 1000.0
    bcfac:
      recharge: 1.0
```

- Key names (e.g., `Upp`) must be unique across the `aquifers` block.
- `func` lists one or more unit response functions to superimpose for this aquifer — see Table 1
  in [Solution methodology](methodology.md). Each function name maps to the key of a
  `boundaryconditions` sub-block (below) that drives it.
- Required parameters depend on the chosen function(s) (e.g., `T`, `S`, `C`, `L`).
- `bcfac` (optional) scales a named boundary-condition time series before it drives this aquifer
  — useful for splitting recharge between multiple conceptual sub-models (e.g. a drainage
  component and a groundwater component), as in the [case studies](examples.md). Default is
  `1.0`. It is generally not relevant for a head boundary condition.

Since Flowsim does not simulate exchange of water between aquifers, aquifers with independent
flow can just as well be run in separate Flowsim runs. It is beneficial to include more than one
aquifer in the same run only when they discharge to the *same* head boundary condition (e.g. a
stream gaining flow from both a shallow and a deep aquifer) — in that case Flowsim's output also
includes the *total* flux summed across the aquifers.

### Boundary condition information (mandatory)

```yaml
boundaryconditions:
  recharge:
    file: "./Data/DMI-10065_JB1_MZ_Evacrop_mp_wb_S10.out"
    header: 0
    date: date
    val: Dsum
    dtformat: "%Y-%m-%d"
    convfact: 0.001
  stream_stage:
    file: "./Data/DMI-10065_JB1_MZ_Evacrop_mp_wb_S10.out"
    header: 0
    date: date
    val: S10_WatLvl
    dtformat: "%Y-%m-%d"
    type: head
    convfact: 1.0
```

There must be one sub-block per boundary-condition key referenced under `func` in the `aquifers`
block. Each sub-block is read as a CSV file (via `pandas.read_csv(..., engine='python')`) and
requires:

- `file`: path and name of the CSV file
- `header` (row number for column names) **or** `colnames` (explicit list of column names)
- `date`: name of the datetime column
- `val`: name of the boundary-condition value column
- `dtformat`: format string for parsing the date column (e.g. `"%Y-%m-%d"`)
- `convfact`: multiplier applied to the series to convert it to the units used by the model
  parameters

#### Optional parameters

- `sep`: column delimiter (default `r'\s+|,\s*|;\s*'`)
- `decimal`: decimal character (default `"."`)
- `skiprows`: number of rows to skip (default `0`)

---

## Boundary condition time series requirements

Before simulating, Flowsim checks that:

- every simulation date exists in each boundary-condition time series,
- the time step is constant within each series, and
- all boundary-condition series share the same temporal resolution.

If any check fails, Flowsim aborts with an error message (also written to `flowsim.log`).

To coarsen or refine a Flowsim simulation temporally, coarsen or refine the boundary-condition
input files accordingly.

---

## Simulation of heat conduction

Flowsim can also simulate **heat conduction** in place of groundwater flow, using an analogous
input structure. See the package documentation for details.
