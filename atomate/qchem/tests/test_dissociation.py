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


# db_file = "/global/homes/s/sblau/config/db.json"
db_file = "/Users/samuelblau/Desktop/db.json"
# xyz_file = "../test_files/top_11/PC.xyz"
# charge = 0
xyz_file = "../test_files/top_11/TFSI-.xyz"
charge = -1
allow_additional_charge_separation = True
multibreak = False

if not allow_additional_charge_separation:
    if charge > 0:
        valid_charges = [charge, charge-1]
    elif charge < 0:
        valid_charges = [charge, charge+1]
    else:
        valid_charges = [-1, 0, 1]
else:
    valid_charges = [charge-2, charge-1, charge, charge+1, charge+2]

mol = Molecule.from_file(xyz_file)
mol.set_charge_and_spin(charge=charge)

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

print(len(target_entries))

num_good_entries = 0
for entry in target_entries:
    if "optimized_molecule" in entry["output"]:
        initial_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["input"]["initial_molecule"]),
                                                strategy=OpenBabelNN,
                                                reorder=False,
                                                extend_structure=False)
        final_mol_graph = build_MoleculeGraph(Molecule.from_dict(entry["output"]["optimized_molecule"]),
                                              strategy=OpenBabelNN,
                                              reorder=False,
                                              extend_structure=False)
        if mol_graph.isomorphic_to(initial_mol_graph) and mol_graph.isomorphic_to(final_mol_graph) and mol_graph.molecule.charge == final_mol_graph.molecule.charge and mol_graph.molecule.spin_multiplicity == final_mol_graph.molecule.spin_multiplicity and entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"] == "gdm":
            num_good_entries += 1
            target_entry = entry

if num_good_entries > 1:
    print("WARNING: There are " + str(num_good_entries) + " valid entries to choose from! Currently using the last one...")

fragment_entries = list(
    mmdb.collection.find({
        "formula_pretty": {
            "$in": unique_formulae
        },
        "input.initial_molecule.charge": {
            "$in": valid_charges
        },
        "calcs_reversed.input.rem.method": target_entry["calcs_reversed"][-1]["input"]["rem"]["method"],
        "calcs_reversed.input.rem.basis": target_entry["calcs_reversed"][-1]["input"]["rem"]["basis"],
        "calcs_reversed.input.rem.scf_algorithm": target_entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"],
        "state": "successful"
    }, {
        "formula_pretty": 1,
        "input": 1,
        "output": 1,
        "calcs_reversed.input.rem": 1,
        "task_id": 1,
        "smiles": 1
    }))

print(len(fragment_entries))

def agnostize(entry):
    to_return = {}
    to_return["formula_pretty"] = entry["formula_pretty"]
    if "smiles" in entry:
        to_return["smiles"] = entry["smiles"]
    to_return["final_energy"] = entry["output"]["final_energy"]
    if "orig" in entry:
        to_return["initial_molecule"] = entry["orig"]["molecule"]
    else:
        to_return["initial_molecule"] = entry["input"]["initial_molecule"]
    if "optimized_molecule" not in entry["output"]:
        if entry["calcs_reversed"][-1]["input"]["rem"]["job_type"] != "sp":
            raise AssertionError("Should only fail to find an optimized_molecule entry from a single point calculation!")
        to_return["final_molecule"] = entry["output"]["initial_molecule"]
    else:
        to_return["final_molecule"] = entry["output"]["optimized_molecule"]
    return to_return

unique_fragment_entries = []
agnostic_entries = []
for entry in fragment_entries:
    found_equivalent = False
    for unique_entry in unique_fragment_entries:
        if entry["output"] == unique_entry["output"]:
            found_equivalent = True
    if not found_equivalent:
        unique_fragment_entries += [entry]
        agnostic_entries += [agnostize(entry)]

print(len(unique_fragment_entries))

bond_dissociation = BondDissociationEnergies(agnostize(target_entry), agnostic_entries, allow_additional_charge_separation, multibreak)
print(bond_dissociation.bond_dissociation_energies)

