from rapidos import RapiDOS
from sys import argv

xlim = [-15, 15]

R = RapiDOS()

R.plot_pdos(elements={'Mn': ['d'],
                      'O': ['p'],
                      'Au': ['d']},
            xlim=xlim)
