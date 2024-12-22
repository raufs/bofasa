#!/usr/bin/env python
from setuptools import setup
import os
import traceback

if __name__ == "__main__":
    setup()
    conda_dir = os.path.abspath(os.environ['CONDA_PREFIX']) + '/'
    try:
        bofasa_data_path = conda_dir + 'share/bofasa/db/'
        conda_act_dir = conda_dir + 'etc/conda/activate.d/'
        conda_dea_dir = conda_dir + 'etc/conda/deactivate.d/'
    
        os.system('BOFASA_DB_PATH=' + bofasa_data_path)
            
        os.system('mkdir -p ' + bofasa_data_path)
        os.system('mkdir -p ' + conda_act_dir + ' ' + conda_dea_dir) 
        os.system("echo 'Default conda space for downloading annotation databases.\n' > " + bofasa_data_path + "README.txt")
        
        act_outf = open(conda_act_dir + 'bofasa.sh', 'w')
        act_outf.write('#!/usr/bin/env bash\n\n')
        act_outf.write('export BOFASA_DB_PATH=' + bofasa_data_path + '\n')
        act_outf.close()

        dea_outf = open(conda_dea_dir + 'bofasa.sh', 'w')
        dea_outf.write('#!/usr/bin/env bash\n\n')
        dea_outf.write('unset BOFASA_DB_PATH\n')
        dea_outf.close()
    except:
        outf = open('installation_error.txt', 'w')
        outf.write('conda dir is: ' + conda_dir + '\n')
        outf.write('Had issues setting up space in conda environment for databases!\n')
        outf.write(traceback.format_exc() + '\n')
        outf.close()