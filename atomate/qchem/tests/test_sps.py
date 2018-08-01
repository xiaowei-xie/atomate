# coding: utf-8
# Copyright (c) Materials Virtual Lab.
# Distributed under the terms of the BSD License.

from __future__ import division, print_function, unicode_literals, absolute_import

import os
import unittest
import numpy as np
from pymatgen.core.structure import Molecule
from pymatgen.analysis.bond_dissociation import BondDissociationEnergies
from atomate.qchem.database import QChemCalcDb
from atomate.qchem.firetasks.fragmenter import build_MoleculeGraph, build_unique_fragments, build_unique_molecules, is_isomorphic
from atomate.qchem.firetasks.parse_outputs import QChemToDb

# db_file = "/global/homes/s/sblau/config/db.json"
# # mol = Molecule.from_file("../test_files/top_11/BF4-.xyz")
# # mol = Molecule.from_file("../test_files/top_11/PF6-.xyz")
# mol = Molecule.from_file("../test_files/top_11/FSI-.xyz")
# mol.set_charge_and_spin(charge=-1)

# # build the MoleculeGraph
# mol_graph = build_MoleculeGraph(mol)

# # find all unique fragments
# unique_fragments = build_unique_fragments(mol_graph)

# # build three molecule objects for each unique fragment:
# # original charge, original charge +1, original charge -1
# unique_molecules = build_unique_molecules(unique_fragments, mol.charge)

# unique_formulae = []
# for molecule in unique_molecules:
#     if molecule.composition.reduced_formula not in unique_formulae:
#         unique_formulae.append(molecule.composition.reduced_formula)

# # attempt to connect to the database to later check if a fragment has already been calculated
# mmdb = QChemCalcDb.from_db_file(db_file, admin=True)

# target_entries = list(
#     mmdb.collection.find({
#         "formula_pretty": mol.composition.reduced_formula
#     }, {
#         "formula_pretty": 1,
#         "input": 1,
#         "output": 1,
#         "calcs_reversed.input.rem": 1
#     }))

# num_good_entries = 0
# for entry in target_entries:
#     initial_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["input"]["initial_molecule"]))
#     final_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["output"]["optimized_molecule"]))
#     if is_isomorphic(mol_graph.graph, initial_mol_graph.graph) and is_isomorphic(mol_graph.graph, final_mol_graph.graph) and mol_graph.molecule.charge == final_mol_graph.molecule.charge and mol_graph.molecule.spin_multiplicity == final_mol_graph.molecule.spin_multiplicity:
#         num_good_entries += 1
#         target_entry = entry

# if num_good_entries > 1:
#     print("WARNING: There are " + str(num_good_entries) + " entries to choose from! Currently using the last one...")

# fragment_entries = list(
#     mmdb.collection.find({
#         "formula_pretty": {
#             "$in": unique_formulae
#         },
#         "calcs_reversed.input.rem.method": target_entry["calcs_reversed"][-1]["input"]["rem"]["method"],
#         "calcs_reversed.input.rem.basis": target_entry["calcs_reversed"][-1]["input"]["rem"]["basis"],
#         "calcs_reversed.input.rem.scf_algorithm": target_entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"],
#         "input.job_type": "sp"
#     }, {
#         "formula_pretty": 1,
#         "input": 1,
#         "output": 1,
#         "calcs_reversed": 1,
#         "dir_name": 1
#     }))

# print(len(fragment_entries))
# print()

# for entry in fragment_entries:
#     print(entry)
#     print()

# # unique_fragment_entries = []
# # for entry in fragment_entries:
# #     found_equivalent = False
# #     for unique_entry in unique_fragment_entries:
# #         if entry["output"] == unique_entry["output"]:
# #             found_equivalent = True
# #     if not found_equivalent:
# #         unique_fragment_entries += [entry]

# # print(len(unique_fragment_entries))

# # bond_dissociation = BondDissociationEnergies(target_entry, unique_fragment_entries)
# # print(bond_dissociation.bond_dissociation_energies)


db_file = "/global/homes/s/sblau/config/db.json"
mmdb = QChemCalcDb.from_db_file(db_file, admin=True)
target_entries = list(
    mmdb.collection.find({
        "input.job_type": "sp"
    }, {
        "calc_dir": 1
    }))

print(len(target_entries))
for entry in target_entries:
    print()
    print(entry)

