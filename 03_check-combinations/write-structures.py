import click
import pathlib

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
    "--environment",
    "environments",
    type=str,
    help="The name of the environments to compare to.",
    multiple=True,
)
@click.option(
    "--output-directory",
    type=str,
    help="The directory to write the structures to.",
)
def get_structures(
    batch: str,
    target_name: str,
    environments: list[str],
    output_directory: str
):

    from openff.toolkit.topology.molecule import Molecule, unit
    import MDAnalysis as mda

    batch_directory = pathlib.Path("targets") / batch
    qm = batch_directory / f"{target_name}.xyz"

    output_directory = pathlib.Path(output_directory) / f"{batch}_{target_name}"
    output_directory.mkdir(parents=True, exist_ok=True)

    u_qm = mda.Universe(qm)
    u_qm.atoms.write(str(output_directory / "qm.xyz"))

    for environment in environments:
        xyz = f"001_fit-all/{environment}/rep1/optimize.tmp/{batch}/iter_0000/{target_name}_mmopt.xyz"

        u = mda.Universe(xyz)
        u.atoms.write(str(output_directory / f"{environment}.xyz"))


if __name__ == "__main__":
    get_structures()