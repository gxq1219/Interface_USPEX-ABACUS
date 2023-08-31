from .Atomistic.Element import Element
from .Atomistic.CellUtility import Cell
from .Atomistic.AtomicPrimitives import AtomicStructure, AtomicDisassembler
from .IO.AtomisticRepresentation import AtomisticRepresentation
AtomisticRepresentation.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Optimizers.GlobalOptimizer import GlobalOptimizer
from .Fitness.Fitness import Fitness
GlobalOptimizer.setFitnessType(Fitness)
from .Selection.USPEXClassic import USPEXClassic
GlobalOptimizer.registerSelection(USPEXClassic)
from .Atomistic.CompositionSpace import CompositionSpace
from .Atomistic.EnvironmentUtility import EnvironmentUtility
EnvironmentUtility.setRepresentation(AtomisticRepresentation)
EnvironmentUtility.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Atomistic.RadialDistributionUtility import RadialDistributionUtility
RadialDistributionUtility.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Atomistic.CellUtility import CellUtility
CellUtility.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Atomistic.SimpleMoleculeUtility import SimpleMoleculeUtility
SimpleMoleculeUtility.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Atomistic.Conditions import Conditions
from .Atomistic.BondUtility import BondUtility
BondUtility.registerTypes(Element, AtomicDisassembler)
from .Atomistic.ElasticML import ElasticML
ElasticML.registerTypes(AtomicDisassembler)
from .Atomistic.Constraints import Constraints
from .XRay.PowderSpectrumAnalyzer import PowderSpectrumAnalyzer
PowderSpectrumAnalyzer.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .XRay.SingleCrystalSpectrumAnalyzer import SingleCrystalSpectrumAnalyzer
SingleCrystalSpectrumAnalyzer.registerTypes(AtomicStructure, Element, Cell, AtomicDisassembler)
from .Atomistic.Operators.Heredity import Heredity
from .Atomistic.Operators.RandTop import RandTop
from .Atomistic.Operators.RandSym import RandSym
from .Atomistic.Operators.RandSymPyXtal import RandSymPyXtal
from .Atomistic.Operators.Softmodemutation import Softmodemutation
from .Atomistic.Operators.Permutation import Permutation
from .Atomistic.Operators.Transmutation import Transmutation
from .Atomistic.Operators.AddAtom import AddAtom
from .Atomistic.Operators.RemoveAtom import RemoveAtom
from .Atomistic.Operators.TeleportAtom import TeleportAtom
from .Atomistic.Operators.Seeds import Seeds
Seeds.registerTypes(AtomisticRepresentation)
GlobalOptimizer.registerTarget('Atomistic',
                      utilities=[CompositionSpace, RadialDistributionUtility, CellUtility, EnvironmentUtility,
                                 SimpleMoleculeUtility, Conditions, BondUtility, Constraints, ElasticML,
                                 PowderSpectrumAnalyzer, SingleCrystalSpectrumAnalyzer, ],
                      hybridizations=[Heredity],
                      mutations=[Softmodemutation, Permutation, Transmutation, AddAtom, RemoveAtom, TeleportAtom],
                      creations=[RandTop, RandSym, RandSymPyXtal],
                      seeds=Seeds)
from .Stages.Executor import Executor
from .Stages.Interfaces.ASEInterfaceAdapter import ASEInterfaceAdapter
ASEInterfaceAdapter.registerTypes(AtomicStructure, Element, Cell)
from .Stages.Interfaces.ABINIT_Interface import ABINIT_Interface
ABINIT_Interface.registerTypes(AtomicStructure, Element, Cell)
Executor.registerInterface('abinit', ABINIT_Interface)
from .Stages.Interfaces.GULP_Interface import GULP_Interface
GULP_Interface.registerTypes(AtomicStructure, Element, Cell)
Executor.registerInterface('gulp', GULP_Interface)
from .Stages.Interfaces.LAMMPS_Interface import LAMMPS_Interface
LAMMPS_Interface.registerTypes(AtomisticRepresentation, ASEInterfaceAdapter.LAMMPS)
Executor.registerInterface('lammps', LAMMPS_Interface)
from .Stages.Interfaces.MLIP_Interface import MLIP_Interface
MLIP_Interface.registerTypes(AtomisticRepresentation, AtomicDisassembler)
Executor.registerInterface('mlip', MLIP_Interface)
from .Stages.Interfaces.PWmat_Interface import PWmat_Interface
PWmat_Interface.registerTypes(AtomicStructure, Element, Cell)
from .Stages.Interfaces.QE_Interface import QE_Interface
QE_Interface.registerTypes(ASEInterfaceAdapter.QE)
Executor.registerInterface('qe', QE_Interface)

from .Stages.Interfaces.ABACUS_Interface import ABACUS_Interface
ABACUS_Interface.registerTypes(ASEInterfaceAdapter.ABACUS)
Executor.registerInterface('abacus', ABACUS_Interface)

from .Stages.Interfaces.VASP_Interface import VASP_Interface
VASP_Interface.registerTypes(ASEInterfaceAdapter.VASP)
Executor.registerInterface('vasp', VASP_Interface)

from .Stages.Interfaces.MOPAC_Interface import MOPAC_Interface
MOPAC_Interface.registerTypes(AtomicStructure, Element, Cell)
Executor.registerInterface('mopac', MOPAC_Interface)
from .Stages.Interfaces.FHIaims_Interface import FHIaims_Interface
FHIaims_Interface.registerTypes(AtomicStructure, Element, Cell)
Executor.registerInterface('aims', FHIaims_Interface)
from .Stages.TaskManagers.BSUB import BSUB
Executor.registerTaskManager('BSUB', BSUB)
from .Stages.TaskManagers.QSUB import QSUB
Executor.registerTaskManager('QSUB', QSUB)
from .Stages.TaskManagers.SBATCH import SBATCH
Executor.registerTaskManager('SBATCH', SBATCH)
from .Stages.TaskManagers.SHELL import SHELL
Executor.registerTaskManager('SHELL', SHELL)
from .Optimizers.ModelOptimizer import ModelOptimizer, External
External.setExecutorType(Executor)
ModelOptimizer.registerModel(External)
ModelOptimizer.registerTarget('Atomistic',
                      utilities=[CompositionSpace, RadialDistributionUtility, CellUtility, EnvironmentUtility,
                                 SimpleMoleculeUtility, Conditions, BondUtility, Constraints, ElasticML,
                                 PowderSpectrumAnalyzer, SingleCrystalSpectrumAnalyzer, ],
                      hybridizations=[Heredity],
                      mutations=[Softmodemutation, Permutation, Transmutation, AddAtom, RemoveAtom, TeleportAtom],
                      creations=[RandTop, RandSym, RandSymPyXtal],
                      seeds=Seeds)
from .Stages.AtomisticStage import AtomisticStage
AtomisticStage.registerTypes(Executor, AtomicDisassembler)
from .Stages.PopulationProcessor import PopulationProcessor
from .Stages.GenerationController import GenerationController
GenerationController.registerOptimizer(GlobalOptimizer)
GenerationController.registerOptimizer(ModelOptimizer)
from .Stages import Stages
Stages.registerStage('execute', Executor)
Stages.registerStage('atomistic', AtomisticStage)
Stages.registerStage('populationProcessor', PopulationProcessor)
GenerationController.setPopulationProcessor(PopulationProcessor)
PopulationProcessor.setStages(Stages)
