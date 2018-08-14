# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import numpy as np
from pymatgen.core.structure import Molecule
from pymatgen.analysis.bond_dissociation import BondDissociationEnergies
from pymatgen.analysis.local_env import OpenBabelNN
from atomate.qchem.database import QChemCalcDb
from atomate.qchem.firetasks.fragmenter import FragmentMolecule
from pymatgen.analysis.graphs import build_MoleculeGraph


db_file = "/global/homes/s/sblau/config/db.json"
# mol = Molecule.from_file("../test_files/top_11/BF4-.xyz")
# mol = Molecule.from_file("../test_files/top_11/PF6-.xyz")
# mol = Molecule.from_file("../test_files/top_11/FSI-.xyz")
mol = Molecule.from_file("../test_files/top_11/TFSI-.xyz")
mol.set_charge_and_spin(charge=-1)

# build the MoleculeGraph
mol_graph = build_MoleculeGraph(mol,
                                strategy=OpenBabelNN,
                                reorder=False,
                                extend_structure=False)

FM = FragmentMolecule()
FM.mol = mol
FM.unique_fragments = mol_graph.build_unique_fragments()
FM._build_unique_relevant_molecules()

unique_formulae = []
for molecule in FM.unique_molecules:
    if molecule.composition.reduced_formula not in unique_formulae:
        unique_formulae.append(molecule.composition.reduced_formula)


mmdb = QChemCalcDb.from_db_file(db_file, admin=True)

target_entries = list(
    mmdb.collection.find({
        "formula_pretty": mol.composition.reduced_formula
    }, {
        "formula_pretty": 1,
        "input": 1,
        "output": 1,
        "calcs_reversed.input.rem": 1
    }))

num_good_entries = 0
for entry in target_entries:
    initial_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["input"]["initial_molecule"]),
                                            strategy=OpenBabelNN,
                                            reorder=False,
                                            extend_structure=False)
    final_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["output"]["optimized_molecule"]),
                                          strategy=OpenBabelNN,
                                          reorder=False,
                                          extend_structure=False)
    if mol_graph.isomorphic_to(initial_mol_graph) and mol_graph.isomorphic_to(final_mol_graph) and mol_graph.molecule.charge == final_mol_graph.molecule.charge and mol_graph.molecule.spin_multiplicity == final_mol_graph.molecule.spin_multiplicity and entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"] == "gdm":
        # print(entry)
        num_good_entries += 1
        target_entry = entry

if num_good_entries > 1:
    print("WARNING: There are " + str(num_good_entries) + " valid entries to choose from! Currently using the last one...")

fragment_entries = list(
    mmdb.collection.find({
        "formula_pretty": {
            "$in": unique_formulae
        },
        "calcs_reversed.input.rem.method": target_entry["calcs_reversed"][-1]["input"]["rem"]["method"],
        "calcs_reversed.input.rem.basis": target_entry["calcs_reversed"][-1]["input"]["rem"]["basis"],
        "calcs_reversed.input.rem.scf_algorithm": target_entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"]
    }, {
        "formula_pretty": 1,
        "input": 1,
        "output": 1,
        "calcs_reversed.input.rem": 1,
        "task_id": 1,
        "dir_name": 1
    }))

print(len(fragment_entries))

# missing_tasks = [2376, 2617, 2427, 2613]

# for entry in fragment_entries:
#     # print(entry["task_id"])
#     if entry["task_id"] in missing_tasks:
#         print("Found missing task " + str(entry["task_id"]) + "!")

unique_fragment_entries = []
for entry in fragment_entries:
    found_equivalent = False
    for unique_entry in unique_fragment_entries:
        if entry["output"] == unique_entry["output"]:
            found_equivalent = True
    if not found_equivalent:
        unique_fragment_entries += [entry]

print(len(unique_fragment_entries))

# for entry in unique_fragment_entries:
#     # print(entry["task_id"])
#     if entry["task_id"] in missing_tasks:
#         print("Found missing task " + str(entry["task_id"]) + "!")

bond_dissociation = BondDissociationEnergies(target_entry, unique_fragment_entries)
print(bond_dissociation.bond_dissociation_energies)

