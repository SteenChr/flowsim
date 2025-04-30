# Flowsim Input Instructions and Output

## Overview

To use Flowsim, all input must be in consistent units. Correspondingly, all output from Flowsim will be presented in these consistent units. The simulation is carried out at discrete time steps with constant duration, defined by the time step of the boundary condition input.

All input must be provided in a required text file named `flowsim.yaml`. This file follows the [YAML](https://yaml.org/) format, a human-readable data serialization standard. It must be located in the working directory when Flowsim is executed.

- Section 5.1: Describes the structure and content of `flowsim.yaml`
- Section 5.2: Describes the output from Flowsim

---

## 5.1 Structure of Input File `flowsim.yaml`

The `flowsim.yaml` file may include:

- Blank lines
- Comments (start with `#`)
- Keys and key-value pairs
- Value continuations across multiple lines

### Value Types

A value in the file can be:

- Text string (e.g., `S10_WatLvL`)
- Integer (e.g., `23`)
- Float (e.g., `7.831` or `3.52e-7`)
- Date (`%Y-%m-%d`, e.g., `1998-08-31`)
- Boolean (`True` / `False`)
- List (e.g., `[3.2, 7.32, 7.5e-3]`)
- Dictionary (e.g., `{key_1: 1.3, key_2: 2.7}`)

### Syntax Rules

- Lists begin with `[` and end with `]`, entries separated by commas
- Dictionaries begin with `{` and end with `}`, key-value pairs separated by commas
- Use **spaces** for indentation (not tabs!)

---

## 5.1.1 Input to Make a Simulation

Only one key-value pair is mandatory:

```yaml
simulation_periods:
  - {begin: "2010-01-01", end: "2015-12-31"}
  - {begin: "2015-09-03", end: "2016-11-20", early_begin: "2005-09-03"}
```

### Optional Setup Keys

```yaml
dtformat: "%Y-%m-%d"
response: "head"  # Default is "flux"
x: [0.0, 58.03, 122.89, 247.76, 375.01]
```

### Aquifer Information (Mandatory)

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

- Key names (e.g., `Upp`) must be unique
- Function names must match those in Table 1
- Required parameters depend on chosen functions (e.g., `T`, `S`, `C`, `L`)

### Boundary Condition Information (Mandatory)

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

#### Optional Parameters

- `sep`: column delimiter (default `r'\s+|,\s*|;\s*'`)
- `decimal`: decimal character (default `"."`)
- `skiprows`: number of rows to skip (default `0`)

---

## 5.1.2 Boundary Condition Time Series Requirements

Flowsim checks:

- Every simulation date exists in the time series
- Time steps are constant across the series
- All series have the same temporal resolution

If any check fails, Flowsim aborts with an error.

---

## 5.2 Output

### 5.2.1 Standard Output File and Log File

Outputs:

- `flowsim-res.csv`: Simulation results
- `flowsim.log`: Errors and warnings

#### `flowsim-res.csv` Columns

- `date`: Date/time
- Boundary conditions: Named after `boundaryconditions` keys
- Head: `"h_" + aquifer_key + "_x=value"`
- Flux: `"q_" + aquifer_key + "_x=value"`
- Total flux: `"q_tot_x=value"` (sum across aquifers)

#### Flux Units

- Infinite aquifers: `[L²/T]`
- Finite aquifers: `[L/T]` (per unit aquifer area)

> ⚠️ Output files are **overwritten** each run.

---

### 5.2.2 Plotting Simulation Results

Add a `plot` block to `flowsim.yaml`:

#### Example: Plot with Observation

```yaml
plot:
  N4:
    plotseries: ["h_Upp_x=58.03"]
    ylim: [30.0, 33.0]
    obs:
      file: "./Data/N4-corrected_WatLvl.csv"
      header: 0
      date: Time
      dtformat: "%Y-%m-%d"
      val: ['WatLvl']
      convfact: 1.0
      dividewith:
        WatLvl: 1.0
```

#### Example: Plot Flux with Observation

```yaml
plot:
  S2:
    plotseries: ["q_tot_x=0.0", "q_Low_x=0.0"]
    obs:
      file: "../../Q_lokal/HymerOut_S1-S5.txt"
      colnames: ["date", "S1", "S2", "S3", "S4", "S5"]
      date: date
      dtformat: "%Y%m%d%H%M%S"
      val: ['S2']
      convfact: 86.4
      dividewith:
        S2: 6957381.0
        S3: 4586094.0
        S4: 2468278.0
      sep: ","
      decimal: "."
      skiprows: 7
```

To disable plotting without deleting the block, rename `plot` to e.g. `noplot`.

---

## 5.3 Simulation of Heat Conduction

Flowsim can also simulate **heat conduction** in place of groundwater flow. (Details omitted.)

