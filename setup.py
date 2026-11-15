#!/usr/bin/env python
from setuptools import setup
import os
import traceback

if __name__ == "__main__":
    setup()
    
    # Only setup conda environment if we're in a conda environment
    if 'CONDA_PREFIX' in os.environ:
        conda_dir = os.path.abspath(os.environ['CONDA_PREFIX']) + '/'
        try:
            bofasa_data_path = conda_dir + 'share/bofasa/db/'
            conda_act_dir = conda_dir + 'etc/conda/activate.d/'
            conda_dea_dir = conda_dir + 'etc/conda/deactivate.d/'
        
            os.makedirs(bofasa_data_path, exist_ok=True)
            os.makedirs(conda_act_dir, exist_ok=True)
            os.makedirs(conda_dea_dir, exist_ok=True)
            with open(bofasa_data_path + 'README.txt', 'w') as readme_outf:
                readme_outf.write('Default conda space for downloading annotation databases.\n')
            
            act_outf = open(conda_act_dir + 'bofasa.sh', 'w')
            act_outf.write('#!/usr/bin/env bash\n\n')
            act_outf.write('export BOFASA_DB_PATH=' + bofasa_data_path + '\n')
            act_outf.close()

            dea_outf = open(conda_dea_dir + 'bofasa.sh', 'w')
            dea_outf.write('#!/usr/bin/env bash\n\n')
            dea_outf.write('unset BOFASA_DB_PATH\n')
            dea_outf.close()
        except Exception as e:
            outf = open('installation_error.txt', 'w')
            outf.write('conda dir is: ' + conda_dir + '\n')
            outf.write('Had issues setting up space in conda environment for databases!\n')
            outf.write(traceback.format_exc() + '\n')
            outf.close()