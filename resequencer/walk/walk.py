import copy
import logging
from dataclasses import dataclass
from pathlib import Path

from resequencer.addition.add import Addition, _run_addition
from resequencer.pdb import PDB, Chain
from resequencer.pdb.residue import Residue


@dataclass
class ProteinWalk:
    """
    Represents a single protein walk in a biomolecular sequence.

    Attributes
    ----------
    chains: tuple[str, ...]
        The chain type of each nucleic acid strand, expects 2.
    target_chain: str
        Name of target strain that interacts with protein
    steps: int
        Number of walk steps, or number of structures to generate
    start_pos: int
        Position of base interacting with protein
    direction: str
        Direction of walk, either '-' or '+',
        where '+' is towards 3' end and '-' is towards 5' end
    original_geometry: str
        Type of nucleic acid geometry to add.
    """

    chains: tuple[str, ...]
    target_chain: str
    steps: int
    start_pos: int
    direction: str
    original_geometry: str

    @classmethod
    def load_walk_file(cls, file: Path, form: str) -> list:
        """Load and parse a protein walk file to create ProteinWalk objects.
        This class method reads a file containing protein walk specifications,
        where each line defines a single ProteinWalk object.
        Each line mus contain exactly 7 whitespace-separated values.
        Empty lines and comments (lines starting with '#') are ignored.

        Parameters
        ----------
        file : Path
            Path object pointing to the protein walk file to be read
        form: str
            The original geometry form to assign to all additions during the walk.

        Returns
        -------
        list[ProteinWalk]
            A list of ProteinWalk objects created from each line in the file.
            Each ProteinWalk object contains:
            - chains: tuple of two chain identifiers
            - target_chain: single chain
            - steps: number of steps in protein walk
            - start_pos: start position of protein walk
            - direction: direction of protein walk, either towards 3' or towards 5'
            - original_geometry: the geometry to add during the protein walk

        Raises
        ------
        ValueError
            If any line does not contain exactly 6 whitespace-separated values
        ValueError
            If the step column contains a non-integer value
        ValueError
            If the start position column contains a non-integer value
        """
        with open(file, "r") as f:
            lines: list[str] = f.readlines()

        chains: list[str] = []
        walks: list[ProteinWalk] = []

        for idx, line in enumerate(lines):
            line: str = line.strip()
            if not line or line.startswith("#"):
                continue

            if len(line.split()) != 6:
                raise ValueError(
                    f"Error reading Line #{idx + 1}: Expected 6 values separated by whitespace, got {len(line.split())}"
                )

            split_elements = line.split(None, 5)
            chains.append(split_elements[0].lower())
            chains.append(split_elements[1].lower())
            target_chain: str = split_elements[2].lower()
            try:
                steps: int = int(split_elements[3])
            except ValueError:
                raise ValueError(
                    f"Expected numeric value for steps in row {idx + 1}, column 4, got {split_elements[3]}!"
                )
            try:
                start_pos: int = int(split_elements[4])
            except ValueError:
                raise ValueError(
                    f"Expected numeric value for start position in row {idx + 1}, column 5, got {split_elements[3]}!"
                )
            direction: str = split_elements[5]

            walks.append(
                cls(
                    tuple(chains),
                    target_chain,
                    steps,
                    start_pos,
                    direction,
                    form,
                )
            )

        return walks


def run_protein_walk(
    pdb: PDB, input_path: Path, walk_path: Path, output_path: Path, form: str
):
    if not walk_path.is_file():
        raise FileNotFoundError(f"Addition file '{walk_path}' does not exist!")

    walks: list[ProteinWalk] = ProteinWalk.load_walk_file(walk_path, form)

    # ---------------------------- Create path objects --------------------------- #

    output_path = Path(output_path)

    for idx, walk in enumerate(walks):
        logging.info(f"Performing Protein Walk ({idx + 1} of {len(walks)})")

        # Identify which chain is target and which is complementary
        if walk.chains[0] == walk.target_chain:
            target_chain_id: str = walk.chains[0]
            other_chain_id: str = walk.chains[1]
        elif walk.chains[1] == walk.target_chain:
            target_chain_id: str = walk.chains[1]
            other_chain_id: str = walk.chains[0]
        else:
            raise ValueError(
                f"Target chain not found within list of provided chains in walk file: {walk_path}"
            )

        target_chain: Chain = pdb[target_chain_id]  # type: ignore
        other_chain: Chain = pdb[other_chain_id]  # type: ignore

        # Extract residues from target and other chain
        target_residues: list[Residue] = copy.deepcopy(target_chain._residues)
        other_residues: list[Residue] = copy.deepcopy(other_chain._residues)

        # Ensure the start position is a valid residue
        if walk.start_pos - 1 >= len(target_residues):
            raise ValueError(
                f"Start position {walk.start_pos} outside of chain length {len(target_residues)}!"
            )

        # Create a right shift or left shift based on the walk direction
        if walk.direction == "+":
            walked_target_residues: list[Residue] = (
                target_residues[1:] + target_residues[:1]
            )
            walked_other_residues: list[Residue] = (
                other_residues[-1:] + other_residues[:-1]
            )
        elif walk.direction == "-":
            walked_target_residues: list[Residue] = (
                target_residues[-1:] + target_residues[:-1]
            )
            walked_other_residues: list[Residue] = (
                other_residues[1:] + other_residues[:1]
            )
        else:
            # Raise an error if the walk direction is neither "+" or "-"
            raise ValueError(
                f"Expected '+' or '-' for walk direction; Got {walk.direction}"
            )

        if len(walked_target_residues) != len(walked_other_residues):
            raise ValueError(
                "Target Chain '{}' of length {} does not match Other Chain '{}' of length {} ".format(
                    target_chain_id.upper(),
                    len(walked_target_residues),
                    other_chain_id.upper(),
                    len(walked_other_residues),
                )
            )
        logging.info(
            "Original Target Chain: "
            + "-".join(res.residue_name for res in pdb[target_chain_id])  # type: ignore
        )
        logging.info(
            "     New Target Chain: "
            + "-".join(res.residue_name for res in walked_target_residues)
        )
        logging.info(
            "Original Other Chain: "
            + "-".join(res.residue_name for res in pdb[other_chain_id])  # type: ignore
        )
        logging.info(
            "     New Other Chain: "
            + "-".join(res.residue_name for res in walked_other_residues)
        )

        forward_range = list(range(len(walked_target_residues)))
        reverse_range = list(reversed(range(len(walked_other_residues))))
        for idx, reverse_idx in zip(forward_range, reverse_range):
            # ----------------------------- Substitute Target ---------------------------- #
            logging.info(
                "Substituting Target Residue #{}: {} for {}".format(
                    idx + 1,
                    target_chain[idx].residue_name,  # type: ignore
                    walked_target_residues[idx].residue_name,
                )
            )
            target_chain[idx].substitute(  # type: ignore
                target_chain[idx].residue_name,  # type: ignore
                walked_target_residues[idx].residue_name,
            )
            # ----------------------------- Substitute Other ----------------------------- #
            logging.info(
                "Substituting Other Residue #{}: {} for {}".format(
                    reverse_idx + 1,
                    other_chain[reverse_idx].residue_name,  # type: ignore
                    walked_other_residues[reverse_idx].residue_name,
                )
            )
            other_chain[reverse_idx].substitute(  # type: ignore
                other_chain[reverse_idx].residue_name,  # type: ignore
                walked_other_residues[reverse_idx].residue_name,
            )
            # if idx == 0:
            #     logging.warning(
            #         "Running first iteration of loop replaces end of new chains"
            #     )

        if walk.direction == "+":
            idx_a = 0
            idx_b = -1
        else:
            idx_a = -1
            idx_b = 0

        target_residue: str = target_chain[idx_a].residue_name  # type: ignore
        other_residue: str = other_chain[idx_b].residue_name  # type: ignore
        target_chain.remove_at(idx_b)
        other_chain.remove_at(idx_a)
        target_seq: list[str] = [res.residue_name for res in target_chain]
        other_seq: list[str] = [res.residue_name for res in target_chain]

        if walk.direction == "+":
            new_target_seq: list[str] = [target_residue] + [
                res.residue_name for res in target_chain
            ]
            new_other_seq: list[str] = [res.residue_name for res in other_chain] + [
                other_residue
            ]
        else:
            new_target_seq: list[str] = [res.residue_name for res in target_chain] + [
                target_residue
            ]
            new_other_seq: list[str] = [other_residue] + [
                res.residue_name for res in other_chain
            ]

        add = Addition(
            (target_chain.chain_id, other_chain.chain_id),
            (target_seq, other_seq),
            (new_target_seq, new_other_seq),
            walk.original_geometry,
            target_chain_id,
        )

        print(
            "{} {} {}".format(
                add.chains[0],
                "".join([a.strip("D") for a in add.old_seq[0]]),
                "".join([a.strip("D") for a in add.new_seq[0]]),
            )
        )
        print(
            "{} {} {}".format(
                add.chains[1],
                "".join([a.strip("D") for a in add.old_seq[1]]),
                "".join([a.strip("D") for a in add.new_seq[1]]),
            )
        )

        logging.info("Protein Walk Substitution Complete!")
        pdb[target_chain_id] = target_chain
        pdb[other_chain_id] = other_chain

        logging.info("Performing Protein Walk Addition...")
        _run_addition(pdb, [add], input_path, output_path)
