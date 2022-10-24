#!/usr/bin/env python

import os 
import sys
import logging
import argparse
import re
import os.path
from numpy import arange
import numpy as np
from scipy.optimize import curve_fit
from glob import glob
from adjustText import adjust_text
import matplotlib
from scipy import interpolate
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
from scipy.optimize import curve_fit

class PotentialCalculator():
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

    def __init__(self, rootdir, logginglevel, plot, ref_nb_ele, bugged, savepickle):
        # logging setting
        intlogginglevel = PotentialCalculator.isint(logginglevel)
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
        self.ref_nb_ele = ref_nb_ele
        self.plot = plot
        self.bugged = bugged
        self.savepickle = savepickle
        self.logger.debug('initializing, rootdir={}'.format(self.rootdir))
    
    def fetch_data(self, n_electronstr, outcarfilename='OUTCAR', stdoutfilename='vasp.out'):
        n_electron = self.isfloat(n_electronstr)        

        cwd = os.getcwd()
        with open(outcarfilename) as outcar:
            outcar = outcar.readlines()
        outcar = ''.join(outcar)

        try: 
            e_pattern = re.compile(r'entropy=\s+-?\d+.\d+\s+energy\(sigma->0\)\s+=\s+(-?\d+\.\d+)')
            e = re.findall(e_pattern, outcar)[-1]
            e = float(e)
            e_to_plot = e
        except: 
            self.logger.warning('no energy found in OUTCAR under current dir\ncwd:{}'.format(cwd))
            e = None

        try:
            e_fermi_pattern = re.compile(r'E-fermi\s+:\s+(\-\d+\.\d+)')
            e_fermi = re.findall(e_fermi_pattern, outcar)[-1]
            e_fermi = float(e_fermi)
        except:
            self.logger.warning('no E-fermi found in OUTCAR under current dir\ncwd:{}'.format(cwd))
            e_fermi = None
        
        if len(glob(stdoutfilename)) > 1: 
            self.logger.warning('more than one stdoutput file under current dir, keeping the latest one only is recommended\ncwd:{}'.format(cwd))

        with open(glob(stdoutfilename)[-1]) as stdoutput:
            stdoutput = stdoutput.readlines()
        stdoutput = ''.join(stdoutput)
        
        try:
            fermi_shift_pattern = re.compile(r'FERMI_SHIFT\s+=\s+(\d+\.\d+)')
            fermi_shift = re.findall(fermi_shift_pattern, stdoutput)[-1]
            fermi_shift = float(fermi_shift)
        except:
            self.logger.warning('no FERMI_SHIFT found in stdoutput under current dir\ncwd:{}'.format(cwd))
            fermi_shift = None

        self.logger.debug('N:{} e:{} E-fermi:{} FERMI_SHIFT:{}'.format(n_electron, e, e_fermi, fermi_shift))
        return [n_electron, e, e_fermi, fermi_shift]
        

    def fit_g_pot(self):
        os.chdir(self.rootdir)
        subdirs = glob('*/')
        subdirs = map(lambda x:re.sub('/', '', x), subdirs)
        subdirs = [subdir for subdir in subdirs if self.isfloat(subdir)]
        self.logger.debug('subdirlist:{}'.format(subdirs))

        datalist = []
        for subdir in subdirs:
            os.chdir(self.rootdir)
            os.chdir(subdir)
            print(subdir)
            datalist.append(self.fetch_data(subdir))
            os.chdir(self.rootdir)
        datalist = map(np.array, zip(*datalist))
        n_electron, e, e_fermi, fermi_shift = datalist

        vacpot = 4.44
        e_to_plot = e
        print(e_to_plot)
        work_function = -e_fermi - fermi_shift
        SHEpot = work_function - vacpot
        estimated_ref_nb_ele = n_electron.mean()
        if np.where(abs(n_electron-estimated_ref_nb_ele) < 0.01)[0].size == 0:
            self.logger.critical('make sure input ref N_ELECT is calculated')
        if self.ref_nb_ele:
            if not abs(self.ref_nb_ele - estimated_ref_nb_ele) < 0.01:
                self.logger.critical('make sure input ref N_ELECT is correct')
        else:
            self.ref_nb_ele = estimated_ref_nb_ele        
        if self.bugged:
            e = e + fermi_shift * (n_electron - self.ref_nb_ele)
        g = e + work_function*(n_electron - self.ref_nb_ele)
        self.logger.debug('pot:{}\nG:{}'.format(SHEpot, g))

        # fitting
        quadratic = lambda x, a, b, c: a * x ** 2 + b * x + c
        parameter, _ = curve_fit(quadratic, SHEpot, g)
        at_SHE_0 = parameter[0] * 0  ** 2 + parameter[1] * 0  + parameter[2] 
        top = -parameter[1]/(2*parameter[0])
        top_y = parameter[0] * ( -parameter[1]/(2*parameter[0])) ** 2 + parameter[1] * (-parameter[1]/(2*parameter[0])) +  parameter[2]
        curvature = 2*parameter[0]/((1+parameter[1]**2 - 4* parameter[2]*parameter[0] + 4*parameter[0]*top_y)**1.5)
        #print(curvature)
        fitted = lambda x: parameter[0] * x ** 2 + parameter[1] * x + parameter[2]
        if self.savepickle:
           pickle.dump(parameter, open('sc.pickle', 'wb'))
        print('a,b,c:', parameter[0],' * x ** 2 +', parameter[1],'* x +' , parameter[2])
        # plotting
        os.chdir(self.rootdir)
        self.logger.debug('figure filename:{}'.format(self.plot))
        x = np.linspace(SHEpot.min()-0.5, SHEpot.max()+1.6, 100)
        y = list(map(fitted, x))
        #print('Energy at 1.23V',(list(map(fitted, [1.6]))))
        for i in x:
            if i == 0:
                print(list(map(fitted, [0])))
            else: continue
        n_electron = np.array(n_electron)
        n_electron = - n_electron + self.ref_nb_ele
        plt.title('a={}\nb={}\nc={}'.format(*parameter) + '\ncurvature={}'.format( curvature))
        fig, ax = plt.subplots()
        plt.scatter(SHEpot, g, label='original data')
        x_ann =  np.linspace(min(x), max(x) , len(n_electron))
        y_ann = np.linspace(min(y), max(g) , len(n_electron)) 
        for i, txt in enumerate(n_electron):
            print(g[i]*(g[i]-g[i-1]))
            ax.annotate(round(txt,2), xy = (SHEpot[i], g[i]), xytext= (SHEpot[i]*1.2 , y_ann[i]) , arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3")) 
        plt.plot([top,top],[min(y),top_y] , ':', color = 'thistle')
        plt.plot([min(x),0],[at_SHE_0,at_SHE_0], '-.', color = 'teal')
        plt.plot([min(x),0],[at_SHE_0,at_SHE_0], '-.', color = 'teal')
        plt.plot([0,0],[min(y),at_SHE_0], '-.', color = 'teal')
        plt.ylabel('G (eV)')
        plt.xlabel('U vs SHE (V)')
        plt.plot(x, y, label='fitted', color = 'lightcoral')
        plt.legend()
        plt.savefig(self.plot, bbox_inches='tight')
        plt.clf()
        ####THIS IS FOR PLOTTING THE CHARGE vs POTENTIAL###
        popt, _ = curve_fit(quadratic, SHEpot, n_electron)
        a, b, c = popt
        x_line = arange(min(SHEpot), 0.2, 0.1)
        y_line = objective(x_line, a, b, c)
        at_0 = objective(1.23, a, b, c) ##CHARGE AT 0 POTENTIAL: HER CONDITION
        #plotting the point (0,at_0)
        print('CHARGE AT 1.23V POTENTIAL', at_0)
        plt.scatter(SHEpot, n_electron)
        plt.plot([min(SHEpot),1.5],[at_0,at_0], '--', color = 'blue')
        plt.plot([1.5,1.5],[min(n_electron),at_0], '--', color = 'blue')
        plt.plot(x_line, y_line, '--', color = 'red', label='fitted')
        plt.legend()
        plt.savefig('charge.png', bbox_inches='tight')
        plt.clf()
        ###THIS IS FOR PLOTTING CHARGE vs ENERGY ###
        print(n_electron)
        plt.scatter(e_to_plot, n_electron)
        plt.plot(e_to_plot, n_electron)
        plt.savefig('energy.png',  bbox_inches='tight')


            
        
def objective(x,a,b,c):
    return a * x ** 2 + b * x + c

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='fit the G-pot curve')
    parser.add_argument('-l', '--logging', dest='logginglevel', action='store',
                        help='logging level', default=20)
    parser.add_argument('-r', '--rootdir', dest='rootdir', action='store',
                        help='rootdir', default='.')
    parser.add_argument('-n', '--number', dest='n', action='store',
                        help='reference number of electron', default=None, type=float)
    parser.add_argument('-p', '--plot', dest='plot', action='store',
                        help='figure name', default='g-pot.png')
    def str2bool(v):
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')
    parser.add_argument('-b', '--bugged', dest='bugged', action='store', type=str2bool,
                        help='is vaspsol bugged? currently yes, and an extra correction is needed', default=True)
    parser.add_argument('-s', '--savepickle', dest='savepickle', action='store', type=str2bool,
                        help='whether save a pickle file which contains fitted parabola parameters', default=True)
    args = parser.parse_args()
    pot = PotentialCalculator(args.rootdir, args.logginglevel, args.plot, args.n, args.bugged, args.savepickle) 
    pot.fit_g_pot()
