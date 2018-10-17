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
from pymatgen.analysis.fragmenter import Fragmenter
from pymatgen.analysis.graphs import MoleculeGraph

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# First, define the db_file for previously calculated fragment entries,
# the XYZ file and the charge of the principle molecule of which we aim
# to compute the bond dissociation energies

db_file = "/global/homes/s/sblau/config/db.json"
# db_file = "/Users/samuelblau/Desktop/db.json"
# xyz_file = "../test_files/top_11/PC.xyz"
# charge = 0
# xyz_file = os.path.join(module_dir, "..", "test_files", "top_11", "TFSI-.xyz")
xyz_file = os.path.join(module_dir, "..", "test_files", "top_11", "EC.xyz")
charge = -1
# pcm_dielectric = 65.0
pcm_dielectric = 40.0

# By default, we will consider one level of charge separation. If additional
# charge separation should be considered, set the following parameter to True:
allow_additional_charge_separation = False

# By default, we will only examine bond dissociation energies when breaking
# one bond at a time. If two ring bonds breaking simultaneously should be
# considered, set the following parameter to True:
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

# Build the principle Molecule and MoleculeGraph
mol = Molecule.from_file(xyz_file)
mol.set_charge_and_spin(charge=charge)
mol_graph = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN(), reorder=False, extend_structure=False)

# Connect to the database
mmdb = QChemCalcDb.from_db_file(db_file, admin=True)

# Find all entries in the database for our principle
target_entries = list(
    mmdb.collection.find({
        "formula_pretty": mol.composition.reduced_formula
    }, {
        "formula_pretty": 1,
        "input": 1,
        "output": 1,
        "calcs_reversed.input": 1
    }))

# Then narrow those to an optimized entry which has our desired charge, multiplicity, SCF strategy, and PCM dielectric.
num_good_entries = 0
for entry in target_entries:
    if "optimized_molecule" in entry["output"]:
        initial_mol_graph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["input"]["initial_molecule"]),
                                                                  OpenBabelNN(), reorder=False, extend_structure=False)
        final_mol_graph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["optimized_molecule"]),
                                                                OpenBabelNN(), reorder=False, extend_structure=False)
        if mol_graph.isomorphic_to(initial_mol_graph) and mol_graph.isomorphic_to(final_mol_graph) and mol_graph.molecule.charge == final_mol_graph.molecule.charge and mol_graph.molecule.spin_multiplicity == final_mol_graph.molecule.spin_multiplicity and entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"] == "gdm":
            if pcm_dielectric != 0:
                if "solvent_method" in entry["calcs_reversed"][-1]["input"]["rem"]:
                    if entry["calcs_reversed"][-1]["input"]["rem"]["solvent_method"] == "pcm" and entry["calcs_reversed"][-1]["input"]["solvent"]["dielectric"] == str(pcm_dielectric):
                        num_good_entries += 1
                        target_entry = entry
            else:
                if "solvent_method" not in entry["calcs_reversed"][-1]["input"]["rem"]:
                    num_good_entries += 1
                    target_entry = entry

if num_good_entries == 0:
    print("No good principle entries found! Exiting...")
    raise RuntimeError

# There should only be one principle entry!
if num_good_entries > 1:
    print("WARNING: There are " + str(num_good_entries) + " valid entries to choose from! Currently using the last one...")

# Use the fragmenter in pymatgen to get a list of unique fragments relevant for BDE calculations:
fragments = Fragmenter(molecule=mol, depth=1).unique_fragments

# Convert the list of fragments into a list of formulae which can then
# be used to search our database:
unique_formulae = []
for mol_graph in fragments:
    if mol_graph.molecule.composition.reduced_formula not in unique_formulae:
        unique_formulae.append(mol_graph.molecule.composition.reduced_formula)

# Find all fragment entries in our database using our unique formulae

find_dict = {"formula_pretty": {"$in": unique_formulae},
             "input.initial_molecule.charge": {"$in": valid_charges},
             "calcs_reversed.input.rem.method": target_entry["calcs_reversed"][-1]["input"]["rem"]["method"],
             "calcs_reversed.input.rem.basis": target_entry["calcs_reversed"][-1]["input"]["rem"]["basis"],
             "calcs_reversed.input.rem.scf_algorithm": target_entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"],
             "state": "successful"
             }

if pcm_dielectric != 0:
    find_dict["calcs_reversed.input.solvent.dielectric"] = str(pcm_dielectric)

fragment_entries = list(
    mmdb.collection.find(find_dict, {
        "formula_pretty": 1,
        "input": 1,
        "output": 1,
        "calcs_reversed.input": 1,
        "task_id": 1,
        "smiles": 1
    }))
# where we've only projected out the relevant information for BDE calculations
# and eventual BDE database entries

# The actual BDE analysis is done in Pymatgen. However, the function we are 
# going to call must have principle and fragment entries in a general format
# that is agnostic of how they were obtained. Thus we define the following
# function to simplify and "agnostize" the data we have taken from our database
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
    if pcm_dielectric != 0:
        to_return["pcm_dielectric"] = int(float(entry["calcs_reversed"][-1]["input"]["solvent"]["dielectric"]))
    return to_return

# We now go through and apply this to our fragment entries while also removing
# any duplicates we find:
# NOTE: Is there a better and more "pythonic" way to do this?
agnostic_entries = []
for entry in fragment_entries:
    agnostic_entry = agnostize(entry)
    found_equivalent = False
    for ag_entry in agnostic_entries:
        if agnostic_entry == ag_entry:
            found_equivalent = True
    if not found_equivalent:
        agnostic_entries += [agnostize(entry)]

# Finally, we call the pymatgen analysis BDE function and print the result:
bond_dissociation = BondDissociationEnergies(agnostize(target_entry), agnostic_entries, allow_additional_charge_separation, multibreak)
print(bond_dissociation.bond_dissociation_energies)

