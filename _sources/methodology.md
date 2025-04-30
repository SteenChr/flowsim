# Solution methodology

Figure 1 illustrates the flow conditions that Flowsim can simulate by using the suite of Carslaw and Jaeger (1959) solutions implemented in Flowsim. Thus, it is assumed that the aquifer is homogenous, with both transmissivity (T) and storativity (S) being constant in space and time. Also, the flow in the aquifer is one-dimensional, and the aquifer is either semi-infinite or finite. The aquifer is connected to a hydraulic head boundary at location x=0, and if it is spatially finite, it meets a no-flux boundary at distance x=L. Finally, the aquifer may receive spatially uniform groundwater recharge (R).

In Flowsim, a Carslaw and Jaeger solution simulates how hydraulic head or groundwater flux in the aquifer responds temporally and spatially to an instantaneous unit change in either groundwater recharge, or in the boundary hydraulic head (bh) at 0x=. Because it simulates the head or flux response to a unit change in a boundary condition, a Carslaw and Jaeger solution is also called a unit response function.

Figure 1 Semi-infinite, or finite, homogeneous aquifer, which is connected to a head dependent boundary at location x=0. For the finite aquifer, it meets a no-flux boundary at distance x=L. The aquifer receives spatially uniform recharge.

Because all the response functions in Flowsim are linear, their solutions can be superimposed to simulate flow response to time series of recharge, time series of boundary head, or time series of both recharge and boundary head.

The following subchapters 3.1 and 3.2 document how superposition of response function solutions is used by Flowsim to simulate time series of hydraulic head and groundwater flux, respectively. Subchapter 3.3 summarizes the Carslaw and Jaeger solutions implemented in Flowsim, while Appendix A gives a more detailed description of each solution.

## Simulation of hydraulic head

Assume that the initial condition of the aquifer is hydraulic head 
```
h(x, 0) = 0
```
and that hydraulic head at the boundary is 

# Documentation of Flowsim

## 3.1 Head Response to Boundary and Recharge Stress

Assume that the initial condition of the aquifer is hydraulic head:

```
h(x, 0) = 0
```

and that hydraulic head at the boundary is:

```
h_b(0) = 0
```

Furthermore, assume that the aquifer at time `t = 0` is stressed by an **instantaneous unit change** of hydraulic head at the boundary, so:

```
h_b(t ≥ 0) = 1
```

In this case, the hydraulic head response in the aquifer equals the **head unit response function**, symbolized by `h_unit(x, t)`. Appendix A gives a number of response functions, each valid for a particular combination of boundary conditions. The available boundary condition combinations are summarized in **Table 1** below.

If a time interval is subdivided into a series of `n` smaller intervals of same duration, Δt, then the series of discrete times can be symbolized by the vector:

```
t = [0, Δt, 2Δt, ..., nΔt]
```

The corresponding **head unit response vector** is defined as:

```
h_unit = [h_unit(x, 0), h_unit(x, Δt), ..., h_unit(x, nΔt)]      (1)
```

Now suppose that the hydraulic head at the `x = 0` boundary changes instantaneously at each time in `t`, and that the corresponding series of boundary head values is:

```
h_b = [h_b,0, h_b,1, h_b,2, ..., h_b,n],   where h_b,0 = 0
```

In other words, the boundary head changes stepwise:
- from `h_b,0` to `h_b,1` at `t = 0`
- from `h_b,1` to `h_b,2` at `t = Δt`
- from `h_b,i` to `h_b,i+1` at `t = i·Δt`, etc.

In this case, at the end of time step `i`, the hydraulic head at distance `x` can be computed by:

```
h(x, iΔt) = ∑_{j=1}^{i} (h_b,j - h_b,j-1) · h_unit(x, (i-j)Δt)     (2)
```

Flowsim uses (2) to simulate the head response of an aquifer stressed by a variation of hydraulic head at the `x = 0` boundary.

If the aquifer is instead stressed by **time-step variation in recharge**, then the head response of the aquifer is computed by:

```
h(x, iΔt) = ∑_{j=1}^{i} (R_j - R_{j-1}) · h_unit(x, (i-j)Δt)      (3)
```

where:

```
R = [R_0, R_1, ..., R_n]
```

Note: The unit response function `h_unit` used in (3) is **different** from the one used in (2).

Flowsim can simulate the response of an aquifer that is stressed both by variation in `h_b` and variation in recharge `R`. This is done by aggregating (2) and (3), giving:

```
h(x, iΔt) = h_b(x, iΔt) + h_R(x, iΔt)                            (4)
```

## 3.2 Simulation of Groundwater Flux

For the type of aquifer illustrated in Figure 1, the **flux of groundwater** at a given location and time is computed as:

```
q(x, t) = -T · ∂h(x, t)/∂x                                      (5)
```

Besides computing the head unit response vector (1), Flowsim can compute the corresponding **hydraulic gradient vector**:

```
∂h_unit/∂x = [∂h_unit(x, 0)/∂x, ∂h_unit(x, Δt)/∂x, ..., ∂h_unit(x, nΔt)/∂x]   (6)
```

For flow driven by time-step changing **boundary head**, Flowsim computes a time series of groundwater flux as:

```
q(x, iΔt) = ∑_{j=1}^{i} (h_b,j - h_b,j-1) · (-T) · ∂h_unit(x, (i-j)Δt)/∂x   (7)
```

For flow driven by time-step changing **recharge**, the flux time series is calculated as:

```
q(x, iΔt) = ∑_{j=1}^{i} (R_j - R_{j-1}) · (-T) · ∂h_unit(x, (i-j)Δt)/∂x     (8)
```

The combined, aggregated flux (driven by both change in boundary head and change in recharge) is computed as:

```
q(x, iΔt) = q_b(x, iΔt) + q_R(x, iΔt)                              (9)
```

## 3.3 Summary of Unit Response Functions Available in Flowsim

Flowsim contains eight solutions (unit response functions) from **Carslaw and Jaeger (1959)**. The eight solutions are presented in **Appendix A** and summarized in Table 1.

The first four solutions are valid for a **semi-infinite aquifer**, while the last four are valid for an aquifer of **finite extension**. For each solution:

1. Either `h_b` or `R` is time-variable (the other is constant).
2. The hydraulic connection between the head boundary and the aquifer at `x = 0` is either **perfect** or **leaky**.

For a **perfect head boundary**:

```
h(0, t) = h_b(t)                                                 (10)
```

For a **leaky connection**:

```
∂h(x, t)/∂x |_(x=0) = (h_b(t) - h(0, t)) · C / T                 (11)
```

### Table 1: Unit Response Functions from Carslaw and Jaeger (1959)

| Function name       | Aquifer type  | Head at boundary | Recharge        | Connection       | Parameters         | Reference |
|---------------------|---------------|------------------|------------------|------------------|---------------------|-----------|
| `sinf_head_perf`    | Semi-infinite | Variable         | Constant = 0     | Perfect          | T, S               | (A-1)     |
| `sinf_head_leak`    | Semi-infinite | Variable         | Constant = 0     | Leaky            | T, S, C            | (A-6)     |
| `sinf_rech_perf`    | Semi-infinite | Constant = 0     | Variable         | Perfect          | T, S               | (A-7)     |
| `sinf_rech_leak`    | Semi-infinite | Constant = 0     | Variable         | Leaky            | T, S, C            | (A-9)     |
| `fin_head_perf`     | Finite        | Variable         | Constant = 0     | Perfect          | T, S, L            | (A-10)    |
| `fin_head_leak`     | Finite        | Variable         | Constant = 0     | Leaky            | T, S, C, L         | (A-11)    |
| `fin_rech_perf`     | Finite        | Constant = 0     | Variable         | Perfect          | T, S, L            | (A-12)    |
| `fin_rech_leak`     | Finite        | Constant = 0     | Variable         | Leaky            | T, S, C, L         | (A-13)    |

**Units:**
- **Transmissivity (T):** \[L²/T\]
- **Storativity (S):** \[-\]
- **Conductance (C):** \[L/T\]
- **Length (L):** \[L\]
