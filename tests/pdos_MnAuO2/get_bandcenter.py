from rapidos import RapiDOS
from sys import argv

xlim = [-10, 3]

R = RapiDOS()

center_d = R.get_pdos_center(elements={'Mn': ['d']})
center_p = R.get_pdos_center(elements={'O': ['p']})

print(center_d, center_p)
