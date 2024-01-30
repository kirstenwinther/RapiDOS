import subprocess
import os
import yaml
from yaml import Dumper
import click
import six

from rapidos import RapiDOS


@click.group()
@click.version_option()
def cli():
    pass


@cli.command()
@click.option('--orbitals', '-o',
              nargs=2,
              type=str,
              multiple=True,
              show_default=True,
              help="""Selected orbitals for plotting the projected DOS.
Input should be two values: "element" and "orbital type"
and additional entries can be given to include several
orbitals in the same plot. For example:

"rapidos plot -o Ru s -o Ru d -o O p",

will plot Ru(s), Ru(d), and O(p)
projections as separate lines.

Accepted input: \n
  element: [Atom symbol, ASE atoms integer] \n
  orbital: [s, p, d, f, t2g, eg, all] \n
""")
@click.option('--limit', '-l',
              nargs=2,
              type=float,
              default=(-20, 10),
              show_default=True,
              help='Limits for energy grid: lower and upper bound should be given')

@click.option('--structure-file', '-s',
              type=str,
              default='CONTCAR',
              show_default=True)
def plot(orbitals, structure_file, limit):
    """Plot total or projected dos. To plot total dos:

    "rapidos plot".

    To plot projected dos use "-o/--orbitals" option

    """

    RP = RapiDOS(structure_file=structure_file)
    if not orbitals:
        RP.plot_dos(xlim=limit)
    else:
        elements_dict = {}
        for entry in orbitals:
            if not entry[0] in elements_dict:
                elements_dict[entry[0]] = []
            elements_dict[entry[0]] += [entry[1]]
        RP.plot_pdos(elements_dict, xlim=limit)


@cli.command()
@click.option('--orbitals', '-o',
              nargs=2,
              type=str,
              show_default=True,
              help="""Selected orbitals for getting pdos center.
Input should be two values: "element" and "orbital type"
For example:

"rapidos plot -o Ru d",

will return Ru d band center

Accepted input: \n
  element: [Atom symbol, ASE atoms integer] \n
  orbital: [s, p, d, f, t2g, eg, all] \n
"""
              )
@click.option('--limit', '-l',
              nargs=2,
              type=float,
              default=(-10, 3),
              help='Limits for energy grid integration: lower and upper bound should be given',
              show_default=True)
@click.option('--structure-file', '-s',
              type=str,
              default='CONTCAR',
              show_default=True)
def bandcenter(orbitals, limit, structure_file):
    """
    Calculate pdos band center:

    "rapidos bandcenter -o Ru d -l -10 5"

    """

    elements_dict = {}

    elements_dict[orbitals[0]] = [orbitals[1]]

    RP = RapiDOS(structure_file=structure_file)
    center = RP.get_pdos_center(elements_dict, xlim=limit)

    print('{}({}) Band center: '.format(orbitals[0], orbitals[1]), center)
