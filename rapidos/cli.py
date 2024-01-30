import subprocess
import os
import yaml
from yaml import Dumper
import click
import six

from .rapidos import RapiDOS

@click.group()
@click.version_option()
def cli():
    pass


@cli.command()
@click.option('--orbitals','-o',
              nargs=2,
              type=str,
              multiple=True,
              help=
"""Selected orbitals for plotting the projected DOS.
Input should be two values: "element" and "orbital type"
and additional entries can be given to include several
orbitals in the same plot. For example:

"rapidos plot -o Ru s -o Ru d -o O p",

will plot Ru(s), Ru(d), and O(p)
projections as separate lines.

Accepted input: \n
  element: [Atom symbol, ASE atoms integer] \n
  orbital: [s, p, d, f, t2g, eg, all] \n
""",
              show_default=True)

def plot(orbitals):
    """Plot total or projected dos. To plot total dos:

    "rapidos plot".

    To plot projected dos use "-o/--orbitals" option

    """

    RP = RapiDOS()
    if not orbitals:
        RP.plot_dos()
    else:
        elements_dict = {}
        for entry in orbitals:
            if not entry[0] in elements_dict:
                elements_dict[entry[0]] = []
            elements_dict[entry[0]] += [entry[1]]
        RP.plot_pdos(elements_dict)


@cli.command()
@click.option('--orbitals','-o',
              nargs=2,
              type=str,
              #multiple=True,
              help=
"""Selected orbitals for getting pdos center.
Input should be two values: "element" and "orbital type"
For example:

"rapidos plot -o Ru d",

will return Ru d band center

Accepted input: \n
  element: [Atom symbol, ASE atoms integer] \n
  orbital: [s, p, d, f, t2g, eg, all] \n
""",
              show_default=True)

@click.option('--limit','-l',
              nargs=2,
              type=float,
              #multiple=True,
              default=(-10,5),
              help=
"""
Limits for energy grid: lower and upper bound should be given.
""",
              show_default=True)

def bandcenter(orbitals, limit):
    """
    Calculate pdos band center:

    "rapidos bandcenter -o Ru d -l -10 5"

    """

    elements_dict = {}

    elements_dict[orbitals[0]] = [orbitals[1]]
    print(elements_dict)
    RP = RapiDOS()
    center = RP.get_pdos_center(elements_dict, limit)

    print('{}({}) Band center: '.format(orbitals[0], orbitals[1]), center)
