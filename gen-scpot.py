#!/usr/bin/env python

import os 
import sys
import logging
import argparse
import os.path
import numpy as np
import shutil
import pickle
import re
import subprocess
class SurfaceChargingRun():
    @staticmethod
    def isint(str):
        try:
            return int(str)
        except:
            return False


    @staticmethod
    def isfloat(str):
        try:
            return float(str)
        except:
            return False

    def fetch_ref_nb_ele(self, outcarfilename='OUTCAR'):
        with open(outcarfilename) as outcar:
            outcar = outcar.readlines()
        outcar = ''.join(outcar)
        nelect_pattern = re.compile(r'NELECT\s+=\s+(\d+\.\d+)\s+total')
        nelect = re.findall(nelect_pattern, outcar)[-1]
        nelect = float(nelect)
        return nelect
        
        

    def __init__(self, rootdir, logginglevel, ref_nb_ele, bugged):
        # logging setting
        intlogginglevel = SurfaceChargingRun.isint(logginglevel)
        if intlogginglevel:
            logginglevel = intlogginglevel

        logger = logging.getLogger('SC Calculator')
        logger.setLevel(logginglevel)
        handler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter('%(levelname)s: %(message)s')
        handler.formatter = formatter
        logger.addHandler(handler)

        rootdir = os.path.abspath(rootdir)
        self.rootdir = rootdir
        self.logger = logger
        if ref_nb_ele:
            self.ref_nb_ele = ref_nb_ele
        else:
            self.ref_nb_ele = self.fetch_ref_nb_ele()
            self.logger.debug('nelect fetched from OUTCAR:{}'.format(self.ref_nb_ele))
        self.bugged = bugged
        self.logger.debug('initializing, rootdir={}'.format(self.rootdir))
    
    def setup(self):
        os.chdir(self.rootdir)
        nelectlist = np.linspace(-1.5, 1.5, 7) + self.ref_nb_ele
        for nelect in nelectlist:
            os.chdir(self.rootdir)
            dirname = '{}'.format(nelect)
            nelect_to_update = {'nelect':nelect}
            print(nelect_to_update)
            #dirname_new = '3_'+dirname
            os.mkdir(dirname)
            with open('nelect.pickle', 'wb') as picklefile:
                pickle.dump(nelect_to_update, picklefile)
            files_to_copy_list = ['job.sh','ase-bfgs.py','OUTCAR', 'nelect.pickle', 'POSCAR']
            for file_to_copy in files_to_copy_list:
                shutil.copy(file_to_copy, dirname)
            os.chdir(dirname)
            subprocess.call(['sbatch', 'job.sh'])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='set up surface charging run')
    parser.add_argument('-l', '--logging', dest='logginglevel', action='store',
                        help='logging level', default=20)
    parser.add_argument('-r', '--rootdir', dest='rootdir', action='store',
                        help='rootdir', default='.')
    parser.add_argument('-n', '--number', dest='n', action='store',
                        help='reference number of electron', default=None, type=float)
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    parser.add_argument('-b', '--bugged', dest='bugged', action='store', type=str2bool,
                        help='is vaspsol bugged? currently yes, and an extra correction is needed', default=True)
    args = parser.parse_args()
    scrun = SurfaceChargingRun(args.rootdir, args.logginglevel, args.n, args.bugged) 
    scrun.setup()
