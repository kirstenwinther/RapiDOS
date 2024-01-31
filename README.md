Code for plotting DOS and PDOS for VASP.

## Using the rapidos CLI

Go to directory with VASP output files and run `rapidos` from the	command line:

		$ rapidos plot --help

		$ rapidos bandcenter --help


## Using in a Python script:

```py
from rapidos import RapiDOS
R = RapiDOS()
# Get dos arrays:
column_names, dos_data = R.get_total_dos()

# Alternatively, plot dos directly in pylab:
R.plot_dos(xlim, title)

# Get entire projected (pdos) data:
column_names, pdos_data = R.get_pdos()

# Or plot pdos directly for selected atoms and/or orbitals:
R.plot_pdos(elements={'Fe': ['d'],
	          					'O': ['p']},
		     		xlim,
		     		title)

# To calculate the band center of a specific projection:

R.get_pdos_center(elements={'Fe': ['d'],
	          								'O': ['p']},
		     					xlim)
```
