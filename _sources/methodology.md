# Solution methodology

Figure 1 in the report illustrates the flow conditions simulated by Flowsim.
The aquifer is assumed homogeneous, one-dimensional, and either semi-infinite
or finite. Transmissivity (T) and storativity (S) are constant in space and time.
The aquifer is connected to a head boundary at x = 0, and for finite aquifers a
no-flux boundary is present at x = L. Spatially uniform recharge R may be applied.

Flowsim uses Carslaw and Jaeger (1959) unit response functions. Because these
responses are linear, they can be superimposed to simulate time-variable forcing.

## 3.1 Simulation of hydraulic head

Assume initial head h(x, 0) = 0 and initial boundary head hb(0) = 0.
For a unit instantaneous boundary change hb(t > 0) = 1, the aquifer response is the
head unit response function hunit(x, t).

For discrete times t = [0, Dt, 2Dt, ..., nDt], the head unit response vector is:

$$
h_{unit} = [h_{unit}(x,0), h_{unit}(x,\Delta t), ..., h_{unit}(x,n\Delta t)]
$$

If boundary head changes stepwise, the simulated head at time step i is:

$$
h(x,i\Delta t)=\sum_{j=1}^{i}\left(h_{b,j}-h_{b,j-1}\right) h_{unit}(x,(i-j)\Delta t)
$$

If recharge changes stepwise, head is simulated as:

$$
h(x,i\Delta t)=\sum_{j=1}^{i}\left(R_j-R_{j-1}\right) h_{unit}(x,(i-j)\Delta t)
$$

The combined response is the sum of boundary-head-driven and recharge-driven
responses.

## 3.2 Simulation of groundwater flux

Groundwater flux is computed from the hydraulic gradient:

$$
q(x,t) = -T\,\frac{\partial h(x,t)}{\partial x}
$$

Flowsim computes gradient unit response vectors and superimposes them exactly as
for head, using either boundary-head increments or recharge increments. Combined
flux is the sum of the two contributions.

## 3.3 Unit response functions available in Flowsim

Flowsim includes eight Carslaw and Jaeger (1959) solutions, summarized below.

| Function name | Aquifer type | Head at boundary | Recharge | Connection at head boundary | Parameters |
|---|---|---|---|---|---|
| sinf_head_perf | Semi-infinite | Variable | Constant = 0 | Perfect | T, S |
| sinf_head_leak | Semi-infinite | Variable | Constant = 0 | Leaky | T, S, C |
| sinf_rech_perf | Semi-infinite | Constant = 0 | Variable | Perfect | T, S |
| sinf_rech_leak | Semi-infinite | Constant = 0 | Variable | Leaky | T, S, C |
| fin_head_perf | Finite | Variable | Constant = 0 | Perfect | T, S, L |
| fin_head_leak | Finite | Variable | Constant = 0 | Leaky | T, S, C, L |
| fin_rech_perf | Finite | Constant = 0 | Variable | Perfect | T, S, L |
| fin_rech_leak | Finite | Constant = 0 | Variable | Leaky | T, S, C, L |

Boundary relation for a perfect connection:

$$
h(0,t) = h_b(t)
$$

Boundary relation for a leaky connection:

$$
T\frac{\partial h}{\partial x}\bigg|_{x=0}=C\left(h_b(t)-h(0,t)\right)
$$

Parameter units:

- Transmissivity, T: [L^2/T]
- Storativity, S: [-]
- Conductance, C: [L/T]
- Aquifer length, L: [L]
