import click

import openmm
from forcebalance.openmmio import LocalEnergyMinimizer
import MDAnalysis as mda
import numpy as np
import pathlib
from forcebalance.molecule import Molecule

from forcebalance.smirnoffio import OptGeoTarget_SMIRNOFF
from forcebalance.parser import parse_inputs
from forcebalance.forcefield import FF
from forcebalance.objective import Objective
from forcebalance.optimizer import Optimizer

def get_positions(target_name, other_xml, options, tgt_opts, forcefield):
    from openmm import unit
    from copy import deepcopy
    import numpy as np
    from forcebalance.nifty import printcool_dictionary
    from forcebalance.openmmio import energy_components
    
    target = OptGeoTarget_SMIRNOFF(options, tgt_opts[0], forcefield)
    
    smirnoff_target = target.engines[target_name]
    smirnoff_target.update_simulation()
    self = smirnoff_target
    
    with open(other_xml) as input:
        self.system = openmm.XmlSerializer.deserialize(input.read())
    
    self.create_simulation(**self.simkwargs)
    smirnoff_target.set_positions(0)
    
    self = smirnoff_target
    if self.restraint_frc_index is not None:
        self.set_restraint_positions(0)

    crit = 1e-4
    # if 'force' in LocalEnergyMinimizer.minimize.__doc__:
    #     crit = max(crit, 1e-2)
    # else:
    #     crit = max(crit, 1e-8)
    shot = 0
    steps = int(max(1, -1*np.log10(crit)))
    
    # Get the previous geometry.
    X0 = self.simulation.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.angstrom)[self.realAtomIdxs]
    
    printcool_dictionary(energy_components(self.simulation), title='Energy component analysis before minimization, shot %i' % shot)
    # Minimize the energy.  Optimizer works best in "steps".
    for logc in np.linspace(0, np.log10(crit), steps):
        e_minimized = self.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        print(e_minimized)
        self.simulation.minimizeEnergy(tolerance=10**logc * unit.kilojoule_per_mole, maxIterations=100000)
        
    # check if energy minimization is successful
    # try 1000 times with 10 steps each as openmm minimizer is not very stable at the tolerance
    for _ in range(1000):
        e_minimized = self.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        print(e_minimized)
        self.simulation.minimizeEnergy(tolerance=crit * unit.kilojoule_per_mole, maxIterations=10)
        e_new = self.simulation.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        if abs(e_new - e_minimized) < crit * 10:
            break
    else:
        logger.error("Energy minimization did not converge")
        raise RuntimeError("Energy minimization did not converge")
    # Remove the restraint energy from the total energy if desired.
    groups = set(range(32))
    if self.restraint_frc_index is not None and not include_restraint_energy:
        frc = self.simulation.system.getForce(self.restraint_frc_index)
        groups.remove(frc.getForceGroup())
    printcool_dictionary(energy_components(self.simulation), title='Energy component analysis after minimization, shot %i' % 0)
    S = self.simulation.context.getState(getPositions=True, getEnergy=True, groups=groups)
    # Get the optimized geometry.
    X1 = S.getPositions(asNumpy=True).value_in_unit(unit.angstrom)[self.realAtomIdxs]
    
    M = deepcopy(self.mol[0])
    M += deepcopy(M)
    M.xyzs = [X0, X1]
    if not self.pbc:
        M.align(center=False)
    X1 = M.xyzs[1]
    
    self._update_positions(X1, False)
    pos = self.getContextPosition()
    
    return pos


@click.command()
@click.option(
    "--batch",
    type=str,
    help="The name of the batch to optimize.",
)
@click.option(
    "--target-name",
    type=str,
    help="The name of the target to optimize.",
)
@click.option(
    "--other-environment",
    type=str,
    help="The name of the environment to compare to.",
)
def compare(
    batch: str,
    target_name: str,
    other_environment: str,
):

    other_xml = f"001_fit-all/{other_environment}/rep1/optimize.tmp/{batch}/iter_0000/{target_name}_mminit.xml"
    other_xyz = f"001_fit-all/{other_environment}/rep1/optimize.tmp/{batch}/iter_0000/{target_name}_mmopt.xyz"

    data_directory = pathlib.Path("targets") / batch

    options, all_tgt_opts = parse_inputs("optimize.in")

    tgt_opts = []
    for opts in all_tgt_opts:
        if opts["name"] == batch:
            tgt_opts.append(opts)
    
    forcefield = FF(options)

    opt_positions = get_positions(target_name, other_xml, options, tgt_opts, forcefield)
    print(opt_positions)

    u = mda.Universe(other_xyz)
    print(batch, target_name, np.allclose(u.atoms.positions, opt_positions))


if __name__ == "__main__":
    compare()