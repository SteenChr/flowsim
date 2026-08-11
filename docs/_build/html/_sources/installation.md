# Installing and running Flowsim

## Installation

Flowsim (v. 1.0.0) requires Python 3.8.8 or higher. It also requires the
packages numpy, pandas, matplotlib, and pyyaml.

Flowsim is available from the Python Package Index (PyPI.org):

	pip install flowsim

To install without dependencies:

	pip install --no-deps flowsim

For more information on pip installation, see:
https://pip.pypa.io/en/stable/cli/pip_install/

Flowsim source code and example releases are available at:
https://github.com/SteenChr/flowsim

A brief documentation website is available at:
https://steenchr.github.io/flowsim

## Running Flowsim

There are two ways to run Flowsim.

The first way is from your own Python script, by including:

	from flowsim import flowsim

and calling:

	flowsim.run_model()

The function call supports optional arguments:

	yaml=<name of yaml input file>
	log=<name of log output file>

with defaults:

	yaml='flowsim.yaml'
	log='flowsim.log'

The second way is from the command line:

	python -m flowsim

with optional arguments:

	--yaml <name of yaml input file>
	--log <name of log output file>
