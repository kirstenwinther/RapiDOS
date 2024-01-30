Code for plotting DOS and PDOS for VASP.

## Getting dos from VASP output

Run in directory with VASP output files:


$ rapidos plot --help
$ rapidos bandcenter --help

Or load in Python

```py
from rapidos import RapiDOS
R = RapiDOS()
# Get dos arrays:
columns, dos = R.get_total_dos()

# Plot dos directly in pylab
R.plot_dos(xlim, title)


# Get entire pdos data:
R.get_pdos()

# Plot projected dos for selected atoms and/or orbitals:
R.plot_pdos(elements={'Fe': ['d'],
	             'O': ['p']},
		     xlim,
		     title)

```

