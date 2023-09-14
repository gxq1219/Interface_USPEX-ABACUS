import logging
import numpy as np
from os.path import join as pj, exists as ex
from ase.io.vasp import iread_vasp_out, read_vasp_xml, write_vasp
from ase.io.espresso import read_fortran_namelist, read_espresso_out, write_espresso_in
from ase.io.abacus import write_abacus, read_abacus_results, read_abacus_out
from ase.io import ParseError, read, write
from ase.atoms import Atoms
from ase.constraints import FixAtoms
import os

logger = logging.getLogger(__name__)
class ASEInterfaceAdapter:
    structureType = None
    atomType = None
    cellType = None

    @classmethod
    def registerTypes(cls, structureType, atomType, cellType):
        cls.structureType = structureType
        cls.atomType = atomType
        cls.cellType = cellType

    class Results:

        EV_PER_CUBIC_ANGSTREM_PER_GPA = 1 / 160.21766208

        def __init__(self, atoms):
            self.results = atoms.get_calculator().results
            self.volume = atoms.get_volume()

        def getEnthalpy(self, pressure):
            return self.results['enthalpy'] if 'enthalpy' in self.results else \
                self.results['energy'] + self.volume * pressure * self.EV_PER_CUBIC_ANGSTREM_PER_GPA

    class VASP:

        # VASP output files
        outcar_file = 'OUTCAR'
        oszicar_file = 'OSZICAR'
        contcar_file = 'CONTCAR'
        xml_file = 'vasprun.xml'

        # VASP input files
        incar_file = 'INCAR'
        kpoints_file = 'KPOINTS'
        poscar_file = 'POSCAR'
        potcar_file = 'POTCAR'

        def write(self, structure, fixedIndices, label, calcFolder):
            cell = structure.getCell()
            symbols = np.asarray([el.short_name for el in structure.getAtomTypes()])
            order = np.argsort(symbols)
            atoms = Atoms(symbols[order], structure.getCartesianCoordinates()[order], cell=cell.getCellVectors())
            if len(fixedIndices) > 0:
                atoms.set_constraint(FixAtoms(indices=fixedIndices))
            write_vasp(pj(calcFolder, self.poscar_file), atoms, label=label, direct=True, vasp5=True, long_format=False)
            return {'pbc': cell.getPBC(), 'symbolsOrder': order}

        def read(self, calcFolder, pbc, symbolsOrder):
            try:
                with open(pj(calcFolder, self.outcar_file)) as f:
                    trajectoryAtoms = list(iread_vasp_out(f, None))
            except (KeyError, ParseError):
                with open(pj(calcFolder, self.xml_file)) as f:
                    trajectoryAtoms = list(read_vasp_xml(f))
            trajectory = []
            for atoms in trajectoryAtoms:
                size = len(atoms)
                positions = np.empty((size, 3), dtype=float)
                atomTypes = np.empty(size, dtype=ASEInterfaceAdapter.atomType)
                for i, symbol, position in zip(symbolsOrder, atoms.get_chemical_symbols(), atoms.get_positions()):
                    positions[i] = position
                    atomTypes[i] = ASEInterfaceAdapter.atomType(symbol)
                cell = ASEInterfaceAdapter.cellType(atoms.get_cell().array, pbc)
                structure = ASEInterfaceAdapter.structureType(atomTypes, positions, cell=cell)
                trajectory.append(dict(
                    structure=structure,
                    results=ASEInterfaceAdapter.Results(atoms)
                ))
            return trajectory

    class LAMMPS:

        # LAMMPS files
        data_file = 'STRUC'
        dump_file = 'lammps.dump'

        def write(self, structure, fixedIndices, label, specorder, calcFolder):
            cell = structure.getCell()
            atoms = Atoms([el.short_name for el in structure.getAtomTypes()], structure.getCartesianCoordinates(),
                          cell=cell.getCellVectors())
            filename = pj(calcFolder, self.data_file)
            atoms.write(filename, format='lammps-data', specorder=specorder)
            with open(filename, 'rt') as f:
                content = f.readlines()
            content[0] = f'{label}\n'
            with open(filename, 'wt') as f:
                f.writelines(content)
            return {'pbc': cell.getPBC()}

        def read(self, calcFolder, specorder, pbc):
            if ex(pj(calcFolder, self.dump_file)):
                atoms = read(pj(calcFolder, self.dump_file), format='lammps-dump-text')
                atomTypes = np.array(
                    [ASEInterfaceAdapter.atomType(specorder[i - 1]) for i in atoms.get_atomic_numbers()])
                structure = ASEInterfaceAdapter.structureType(atomTypes, atoms.get_positions(),
                                                              cell=ASEInterfaceAdapter.cellType(atoms.get_cell().array,
                                                                                                pbc))
                return dict(
                    structure=structure,
                    results=ASEInterfaceAdapter.Results(atoms)
                )
            else:
                return None

    class QE:

        inputFile, outputFile = 'input', 'output'

        def __init__(self, options):
            with open(options) as fp:
                data, card_lines = read_fortran_namelist(fp)
            if 'system' not in data:
                raise KeyError('Required section &SYSTEM not found.')
            self.data = data

        def write(self, structure, fixedIndices, kPoints, pseudopotentials, calcFolder):
            cell = structure.getCell()
            atoms = Atoms(symbols=[el.short_name for el in structure.getAtomTypes()],
                          positions=structure.getCartesianCoordinates(),
                          cell=cell.getCellVectors())
            if len(fixedIndices) > 0:
                atoms.set_constraint(FixAtoms(indices=fixedIndices))
            with open(calcFolder / self.inputFile, 'wt') as f:
                write_espresso_in(f,
                                  atoms=atoms, input_data=self.data,
                                  pseudopotentials={s: p.name for s, p in pseudopotentials.items()},
                                  kpts=kPoints,
                                  crystal_coordinates=True)
            return {'pbc': cell.getPBC()}

        def read(self, calcFolder, pbc):
            with open(pj(calcFolder, self.outputFile)) as f:
                atoms = next(read_espresso_out(f, index=slice(None, -2, -1)))
            atomTypes = np.array([ASEInterfaceAdapter.atomType(s) for s in atoms.get_chemical_symbols()])
            structure = ASEInterfaceAdapter.structureType(atomTypes, atoms.get_positions(),
                                                          cell=ASEInterfaceAdapter.cellType(atoms.get_cell().array,
                                                                                            pbc))
            return dict(
                structure=structure,
                results=ASEInterfaceAdapter.Results(atoms)
            )

    class ABACUS:
        # ABACUS input files
        input_file = 'INPUT'
        kpoints_file = 'KPT'
        stru_file = 'STRU'

        # ABACUS output files
        cif_file = 'OUT.USPEX/STRU_NOW.cif'
        log1_file = 'OUT.USPEX/running_cell-relax.log'
        log2_file = 'OUT.USPEX/running_relax.log'
        log3_file = 'OUT.USPEX/running_scf.log'
        output_file = 'OUT.USPEX'        

        def write(self, filename, structure, pp, basis):
            cell = structure.getCell()
            symbols = np.asarray([el.short_name for el in structure.getAtomTypes()])
            order = np.argsort(symbols)
            atoms = Atoms(symbols[order], structure.getCartesianCoordinates()[order], cell=cell.getCellVectors())
            #if len(fixedIndices) > 0:
            #    atoms.set_constrain(FixAtoms(indices=fixedIndices))
            with open(filename, 'wt') as f:
                write_abacus(f, atoms, pp, basis)
            
            return {'pbc': cell.getPBC()}

        def read(self, foldername, pbc):
            for i in os.listdir(foldername):
               if i.startswith('running') and i.endswith('log'):
                   log = i
                   break
               else:
                   logger.info('No log file.')
            atoms = read(pj(foldername,log), format='abacus-out')
            atomTypes = np.array([ASEInterfaceAdapter.atomType(s) for s in atoms.get_chemical_symbols()])
            structure = ASEInterfaceAdapter.structureType(atomTypes, atoms.get_positions(), \
                        cell=ASEInterfaceAdapter.cellType(atoms.get_cell().array, pbc))
            volumeline = os.popen('grep Volume {} |tail -n 1'.format(pj(foldername,log))).read()
            volume = volumeline.strip().split('=')[1]
            return dict(structure=structure, results=ASEInterfaceAdapter.Results(atoms), volume = float(volume))
            
