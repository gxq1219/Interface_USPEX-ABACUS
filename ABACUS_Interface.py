"""
USPEX.Stages.ABACUS_Interface
===========================

"""
import logging
import numpy as np
import shutil
import os
from os.path import join as pj
from typing import List
from ase.io import write,read
from .KPoints import KPoints, BadKPoints
from ase.atoms import Atoms
import glob
from ase.io.abacus import AbacusOutCalcChunk, AbacusOutHeaderChunk

logger = logging.getLogger(__name__)


class ABACUS_Interface:
    '''
    Calculator for ABACUS.
    Local running
    '''


    inputFile, outputFile, errorFile = 'INPUT', 'output', 'error'

    # working input files
    input_file = 'INPUT'
    kpoints_file = 'KPT'
    stru_file = 'STRU'
    potentials_file = 'Specific/ATOMIC_SPECIES'
    orbital_file = 'Specific/NUMERICAL_ORBITAL'
    sub_file='Specific/run.sh'

    with open("Specific/INPUT_1", "r") as file:
        lines = file.readlines()
        suffix_line = None
        for line in lines:
            if 'suffix' in line:
                suffix_line = line
                break
        if suffix_line:
            word = suffix_line.rstrip('\n').split('suffix')[0]
        else:
            word = ['USPEX']

    #### Reading potential
    pseudopotentials = dict()

    with open('Specific/ATOMIC_SPECIES', 'r') as pse:
        Lines =  pse.readlines()
        for line in Lines:
            element = line.rstrip('\n').split()[0]
            mass = line.rstrip('\n').split()[1]
            pot = line.rstrip('\n').split()[2]
            pseudopotentials[element] = pot
    
    ### Reading orbital
    basis = dict()
    with open('Specific/NUMERICAL_ORBITAL','r') as orb:
        Lines_orb = orb.readlines()
        for line_orb in Lines_orb:
            element = line_orb.rstrip('\n').split()[0]
            Orb = line_orb.rstrip('\n').split()[1]
            basis[element] = Orb


    # working output files
    log1_file = 'OUT.USPEX/running_cell-relax.log'
    log2_file = 'OUT.USPEX/running_relax.log'
    log3_file = 'OUT.USPEX/running_scf.log'
    cif_file = 'OUT.USPEX/STRU_NOW.cif'

    DEFAULT_SLEEP_TIME = 30

    aseAdapterType = None

    #print(pseudopotentials)
    @classmethod
    def registerTypes(cls, aseAdapterType):
        cls.aseAdapterType = aseAdapterType

    def __init__(self, tag: str, kresol: float, input: str = None, pseudopotentials = pseudopotentials, basis =  basis ,targetProperties: list = None, **kwargs):
        '''
        :param params: dictionary with parameters:
                * commandExecutable: str of executable command
                * kresol: float of K-points resolution
                * remote: dict of remote server params
                * taskManager: dict of task managers params
        :param step: int of current step
        '''
        self.tag = tag
        if input is not None:
            self.input = input
        else:
            self.input = pj(os.getcwd(), f'Specific/INPUT_{tag}')
        assert os.path.exists(self.input)
        self.output_file='OUT.USPEX'
        self.sub_file='Specific/run.sh'
        self.log_file1='OUT.USPEX/running_cell-relax.log'
        self.log_file2='OUT.USPEX/running_relax.log'
        self.log_file3='OUT.USPEX/running_scf.log'
        #print(pseudopotentials)
        self.adapter = self.aseAdapterType()
        self.kPoints = KPoints(kresol)
        #pseudopotentials_path = 'Specific/'
        self.pseudopotentials = {s : p for s, p in pseudopotentials.items()}
        #print(self.input)

        #for x, p in self.pseudopotentials.items():
        #    assert p.exists()
        with open(self.input, 'r') as file:
            file_content = file.read()
            if 'lcao' not in file_content:
                self.basis = None
            else:
                self.basis = {s: orbital for s, orbital in basis.items()}
        self.failedSystems = []
        self.targetProperties = targetProperties if targetProperties is not None else ['structure', 'enthalpy']


    def prepareLocalCalculation(self, system, calcFolder: str):
        '''
        :param system: our system
        :return:
        '''
        with open(pj(calcFolder, self.input_file), 'wt') as f:
            pass
       
        structure = system['structure']

        with open('Specific/ATOMIC_SPECIES', 'r') as pse:
             Lines =  pse.readlines()
             for line in Lines:
                 upf = line.rstrip('\n').split()[2]
                 source_file = 'Specific' + "/" + upf
                 shutil.copy(source_file, calcFolder)

        if self.basis is not None:
             with open('Specific/NUMERICAL_ORBITAL', 'r') as orb:
                 Lines =  orb.readlines()
                 for line in Lines:
                     ls = line.rstrip('\n').split()[1]
                     source_file = 'Specific' + "/" + ls
                     shutil.copy(source_file, calcFolder)


        ############################## INPUT ################################
        source = self.input 
        dest =  pj(calcFolder, self.input_file)
        shutil.copy2(source, dest)

        if system['externalPressure']:
            with open(pj(calcFolder, self.input_file), 'a') as myfile:
                myfile.write(f"press1 \t {10*system['externalPressure']:10f}\n")
                myfile.write(f"press2 \t {10*system['externalPressure']:10f}\n")
                myfile.write(f"press3 \t {10*system['externalPressure']:10f}\n")
        if calcFolder in self.failedSystems:
            with open(pj(calcFolder, self.input_file), 'a') as myfile: myfile.write('symmetry -1\n')
        ############################# KPT #################################
        try:
            kPoints = self.kPoints.build(structure.getCell())
        except BadKPoints:
            # This LATTICE is extremely wrong, let's skip it from now
            logger.info('K-points cannot be built, so it\'s set as   [1, 1, 1]')
            kPoints = [1, 1, 1]

        with open(pj(calcFolder, self.kpoints_file), 'w') as fp:
            fp.write('K_POINTS\n0\nGamma\n')
            fp.write('%4d %4d %4d   0   0   0\n' % tuple(kPoints))

    ############################# STRU ##################################
        #system['ase'] = self.adapter.write( structure, self.pseudopotentials, self.basis, calcFolder)
        cell = structure.getCell()
        symbols = np.asarray([el.short_name for el in structure.getAtomTypes()])
        order = np.argsort(symbols)
        atoms = Atoms(symbols[order], structure.getCartesianCoordinates()[order], cell=cell.getCellVectors())

        if self.pseudopotentials is not None:
            pse = {element : poten for element, poten in self.pseudopotentials.items() if element in set(symbols)}
        else:
            pse = None

        if self.basis is not None:
            bas = {element : basi for element, basi in self.basis.items() if element in set(symbols)}
        else:
            bas = None
        system['ase']=self.adapter.write(pj(calcFolder,self.stru_file), structure, pp = pse, basis = bas)

        return ''

    # Checking whether the SCF has converged
    def isConverged(self, calcFolder: str):
        filename  = pj(calcFolder, self.output_file)
        for i in os.listdir(filename):
            if i.startswith('running') and i.endswith('log'):
                log = i
                with open(pj(filename,log), 'r') as file:
                    content = file.read()
                    header = AbacusOutHeaderChunk(content)
                    result= AbacusOutCalcChunk(content,header=header)
                break
            else:
                logger.info('No log file.')
 
        return result

############   Read
    def readOutput(self, system: dict, calcFolder: str):
        filename=pj(calcFolder, self.output_file)
        aseResults = self.adapter.read(filename, **system.pop('ase'))
        results = {}
        logger.info(f"aseResults['results'].results.keys()")
        if 'structure' in self.targetProperties:
            results['structure'] = aseResults['structure']
        if 'enthalpy' in self.targetProperties:
            volume = aseResults['volume']
            results['enthalpy'] = aseResults['results'].results['energy'] + system['externalPressure'] * volume
        if 'energy' in self.targetProperties:
            results['energy'] = aseResults['results'].results['energy']
        if 'forces' in self.targetProperties:
            results['forces'] = aseResults['results'].results['forces']
        
        return results 
