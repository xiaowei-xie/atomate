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


db_file = "/global/homes/s/sblau/config/db.json"
mol = Molecule.from_file("../test_files/top_11/BF4-.xyz")
mol.set_charge_and_spin(charge=-1)

# build the MoleculeGraph
mol_graph = build_MoleculeGraph(mol)

# find all unique fragments
unique_fragments = build_unique_fragments(mol_graph)

# build three molecule objects for each unique fragment:
# original charge, original charge +1, original charge -1
unique_molecules = build_unique_molecules(unique_fragments, mol.charge)

unique_formulae = []
for molecule in unique_molecules:
    if molecule.composition.reduced_formula not in unique_formulae:
        unique_formulae.append(molecule.composition.reduced_formula)

# attempt to connect to the database to later check if a fragment has already been calculated
mmdb = QChemCalcDb.from_db_file(db_file, admin=True)

target_entries = list(
    mmdb.collection.find({
        "formula_pretty": mol.composition.reduced_formula
    }, {
        "formula_pretty": 1,
        "input.initial_molecule": 1,
        "output.optimized_molecule": 1,
        "output.final_energy": 1
    }))

for entry in target_entries:
    initial_mol_graph = build_MoleculeGraph(entry["input"]["initial_molecule"])
    final_mol_graph = build_MoleculeGraph(entry["output"]["optimized_molecule"])
    if is_isomorphic(mol_graph.graph, initial_mol_graph.graph) and is_isomorphic(mol_graph.graph, final_mol_graph.graph):
        target_entry = entry

fragment_entries = list(
    mmdb.collection.find({
        "formula_pretty": {
            "$in": unique_formulae
        }
    }, {
        "formula_pretty": 1,
        "input.initial_molecule": 1,
        "output.optimized_molecule": 1,
        "output.final_energy": 1
    }))

dissociation_energies = BondDissociationEnergies(target_entry, fragment_entries)
print(dissociation_energies)

