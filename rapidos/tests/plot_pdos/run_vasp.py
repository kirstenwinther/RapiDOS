import os
exitcode = os.system('srun -n 660 -c 4 --cpu_bind=cores /project/projectdirs/m2997/special_cori_knl_latest_default')
