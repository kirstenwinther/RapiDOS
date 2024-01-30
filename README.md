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

```

