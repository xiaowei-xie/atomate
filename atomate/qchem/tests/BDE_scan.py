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


# remove unnecessary pcm_dielectric param when less lazy
def agnostize(entry, pcm_dielectric):
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


def call_BDE_analysis(molecule, db_file, pcm_dielectric, allow_additional_charge_separation=False, multibreak=False):

    if not allow_additional_charge_separation:
        if charge > 0:
            valid_charges = [charge, charge-1]
        elif charge < 0:
            valid_charges = [charge, charge+1]
        else:
            valid_charges = [-1, 0, 1]
    else:
        valid_charges = [charge-2, charge-1, charge, charge+1, charge+2]

    # Build the principle MoleculeGraph
    mol_graph = MoleculeGraph.with_local_env_strategy(molecule, OpenBabelNN(), reorder=False, extend_structure=False)

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
        print("No good principle entries found! Will look for entries that match everything but final structure in case the initial FF yielded a modified connectivity...")
        for entry in target_entries:
            if "optimized_molecule" in entry["output"]:
                initial_mol_graph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["input"]["initial_molecule"]),
                                                                          OpenBabelNN(), reorder=False, extend_structure=False)
                final_mol_graph = MoleculeGraph.with_local_env_strategy(Molecule.from_dict(entry["output"]["optimized_molecule"]),
                                                                        OpenBabelNN(), reorder=False, extend_structure=False)
                if mol_graph.isomorphic_to(initial_mol_graph) and mol_graph.molecule.charge == final_mol_graph.molecule.charge and mol_graph.molecule.spin_multiplicity == final_mol_graph.molecule.spin_multiplicity and entry["calcs_reversed"][-1]["input"]["rem"]["scf_algorithm"] == "gdm":
                    if pcm_dielectric != 0:
                        if "solvent_method" in entry["calcs_reversed"][-1]["input"]["rem"]:
                            if entry["calcs_reversed"][-1]["input"]["rem"]["solvent_method"] == "pcm" and entry["calcs_reversed"][-1]["input"]["solvent"]["dielectric"] == str(pcm_dielectric):
                                num_good_entries += 1
                                target_entry = entry
                    else:
                        if "solvent_method" not in entry["calcs_reversed"][-1]["input"]["rem"]:
                            num_good_entries += 1
                            target_entry = entry
        if num_good_entries != 0:
            print("Found a potential entry, but beware that connectivity has changed!")

    if num_good_entries == 0:
        print("No good principle entries found! Exiting...")
        raise RuntimeError

    # There should only be one principle entry! But in practice, principles often get duplciated when rerunning a WF that previously had fizzled FWs.
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

    # Remove duplicates and agnostize entries
    agnostic_entries = []
    for entry in fragment_entries:
        agnostic_entry = agnostize(entry, pcm_dielectric)
        found_equivalent = False
        for ag_entry in agnostic_entries:
            if agnostic_entry == ag_entry:
                found_equivalent = True
        if not found_equivalent:
            agnostic_entries += [agnostize(entry, pcm_dielectric)]

    # Finally, we call the pymatgen analysis BDE function and print the result:
    bond_dissociation = BondDissociationEnergies(agnostize(target_entry, pcm_dielectric), agnostic_entries, allow_additional_charge_separation, multibreak)
    return bond_dissociation.bond_dissociation_energies


module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))
db_file = "/global/homes/s/sblau/config/db.json"
xyz_file = os.path.join(module_dir, "..", "test_files", "top_11", "EC.xyz")
charge = -2
mol = Molecule.from_file(xyz_file)
mol.set_charge_and_spin(charge=charge)
eps_vals = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0]

for eps in eps_vals:
    print(eps)
    try:
        vals = call_BDE_analysis(mol, db_file, eps)
        print([entry[0] for entry in vals])
        print()
    except RuntimeError:
        print('ERROR')
        print()

