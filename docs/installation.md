# Installation and Quickstart
## Installation

Flowsim (v. 1.0.0) requires Python 3.8.8 or higher. Furthermore it requires the numpy (v. 1.20.1), pandas (v. 1.2.4), matplotlib (v. 3.3.4), and yaml (pyyaml v. 5.4.1) packages.

Flowsim is available from the Python Package Index (PyPI.org) repository. It is installed by executing from the command line:

&ensp;&ensp;&ensp;&ensp; **pip install flowsim**

This installation uses the following dependencies of other Python packages:

```
dependencies = [
		"numpy >=1.20.1",
		"matplotlib >= 3.3.4",
		"pandas >=1.2.4",
		"pyyaml >=5.4.1"
	]
```

To install Edcrop without caring for dependencies, use instead:

&ensp;&ensp;&ensp;&ensp; **pip install --no-deps flowsim**

The above commands install Flowsim into the Python site-packages directory.

For more information on pip installation, consult the pip documentation at https://pip.pypa.io/en/stable/cli/pip_install/.

Flowsim and example releases can also be found on https://github.com/SteenChr/flowsim.

## Quickstart

Running Flowsim can be run from the user’s own script by including in the script both the statement:

&ensp;&ensp;&ensp;&ensp; `from flowsim import flowsim`

and the function call:

&ensp;&ensp;&ensp;&ensp; `flowsim.run_model()`

The function call may include two optional arguments:

&ensp;&ensp;&ensp;&ensp; `yaml=<name of yaml input file>`

&ensp;&ensp;&ensp;&ensp; `log=<name of log output file>`

with the defaults:
&ensp;&ensp;&ensp;&ensp; `yaml='flowsim.yaml'`

&ensp;&ensp;&ensp;&ensp; `log='flowsim.log'`
