# Flowsim



[Flowsim Online Documentation](https://steenchr.github.io/flowsim/home.html)

Flowsim is a Python package, which uses simple Carslaw and Jaeger (1959) solutions for 
one-dimensional (1D-) transient groundwater flow to simulate time-dependent variation 
in either hydraulic head or groundwater flux at specified locations in the aquifer. 
All the solutions in Flowsim were derived for linear flow problems, which means that 
they can be superimposed (aggregated). The aggregation of solutions can for example 
make a model where flow is driven through a layered system by a mixture of time-dependent
boundary conditions. Thus it is possible to use these simple solutions to design
conceptual, physically meaningful models, which simulate time series of hydraulic 
head or groundwater flux comparable to what can be observed in the field. Compared to 
conducting numerical modeling, Flowsim has the advantage that it takes no time to design, 
set up a model, and run a simulation. It is therefore beneficial to use Flowsim at an 
early stage of a hydrological modeling project, before developing and using more a complex
numerical model that may eventually be needed to finalize the study and make the decisions.
Flowsim can for example be used to support early decisions about numerical model design,
parameter estimation strategy, and how to model generation of groundwater recharge from 
climatic data.

Flowsim can be imported and used in user's own script, or it can be executed from the command 
line as a script.
