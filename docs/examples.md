# Example input files

The Flowsim GitHub repository has a folder named Examples, which has a number of subfolders. Each subfolder contains the input files required by Flowsim to run a specific flow simulation example. A very short description of each example follows below.

---

## Sillerup

This is an example where Flowsim is used to simulate discharge to a stream from a catchment with glacial clayey deposits covering sandy deposits forming a deep aquifer. The land use is mainly agricultural, and due to the shallow clayey deposits the fields contain a dense network of drain systems. The streamflow, which is driven by time-variable recharge, is therefore conceptualized to conisist of both fast-responding drain discharge and slow-reponding deep groundwater discharge (the latter forming the stream's base flow). Both types of discharge are modelled by use of the Flowsim response function named “fin_rech_leak”.

The Sillerup example is used in the report's Section 4.1 to demonstrate that it can be very useful to use a simple tool like Flowsim prior to complex numerical modeling to support early decisions about numerical model design and parameter estimation strategy, including decision about objective function. It can also be very useful to use Flowsim to make early decision about how to model generation of groundwater recharge from climatic data. The advantage of using Flowsim is that it takes no time to set up alternative simple models, and it takes no time (just seconds) to run such a model. Therefore a study like this for Sillerup can be made by an experienced modeler in a few hours, providing a lot of useful information and support for decisions that can be used in a following advanced study based on numerical modeling. 

The Sillerup example files can be found [here](https://github.com/SteenChr/flowsim/tree/main/Examples/Sillerup).

---

## Knivsbaek

The Knivsbaek example is quite similar to the Sillerup example, but for a different kind of hydrogeologic setting. It is used in the report's Section 4.2 to demonstrate that Flowsim can be used to do a quick analysis of a streamflow time series, which separates streamflow into a slow baseflow component and a flashy component. Together with geologic information from boreholes, the analysis supports a conceptual understanding (a hypothesis) that baseflow originates from groundwater recharge in areas without shallow clay, while the flashy part of stream flow mainly originates from groundwater recharge in areas with shallow clay.

The Knivsbaek example files can be found [here](https://github.com/SteenChr/flowsim/tree/main/Examples/Knivsbaek).

---

## Abild aa

This is an example where Flowsim is used to simulate hydraulic head variation in piezometers along a transect of an unconfined aquifer. The transect is perpendicular to a stream. The head variation is driven by both time-variable groundwater recharge and time-variable stream stage. Using Flowsim, in this example the head response due to variation in recharge is simulated by the “fin_rech_leak” response function, while the head response due to stream stage variation is simulated by the “sinf_head_leak” response function.

The example is used in the report's Section 4.3 to demonstrate that when the boundary conditions can be specified, various hydrogeologic parameter values can be estimated from head time series observed at various locations. Transmissivity could mainly be estimated from mean hydraulic head observations, and particularly at distance from the stream; storativity could mainly be estimated from the temporal fluctuations of hydraulic head, and particularly at distance from the stream; and streambed conductance could also mainly be estimated from mean hydraulic head observations, and in particular from observations close to the stream.

The Abild aa example files can be found [here](https://github.com/SteenChr/flowsim/tree/main/Examples/Abild_aa).
