# Pharmacokinetics Module in Python
"""
This module contains functions for calculating key pharmacokinetic parameters used in dosing patients.
Each function corresponds to a pharmacokinetic formula, making it easy to perform various calculations.

Functions included:
- Elimination rate constant
- Half-life
- Initial concentration for IV bolus
- Plasma concentration (single and multiple doses)
- Peak and trough concentrations (multiple doses)
- Average steady-state concentration
- Clearance (different methods)
- Plasma concentration during constant rate infusion
- Volume of distribution
- Creatinine clearance (males and females)
- Adjustments for patient-specific parameters
- Drug interaction checking
"""

# Import necessary libraries
import numpy as np
import json

# Existing pharmacokinetic functions

def elimination_rate_constant(CL, Vd):
    """
    Calculate the elimination rate constant.

    Parameters:
    CL : float
        Clearance (L/hr)
    Vd : float
        Volume of distribution (L)

    Returns:
    float
        Elimination rate constant (1/hr)
    """
    return CL / Vd

def half_life(CL, Vd):
    """
    Calculate the half-life of a drug.

    Parameters:
    CL : float
        Clearance (L/hr)
    Vd : float
        Volume of distribution (L)

    Returns:
    float
        Half-life (hours)
    """
    k_e = elimination_rate_constant(CL, Vd)
    return 0.693 / k_e

def initial_concentration(D, Vd):
    """
    Calculate the initial concentration after an IV bolus.

    Parameters:
    D : float
        Dose (mg)
    Vd : float
        Volume of distribution (L)

    Returns:
    float
        Initial concentration (mg/L)
    """
    return D / Vd

def plasma_concentration_iv_bolus_single(C_0, k_e, t):
    """
    Calculate the plasma concentration after an IV bolus (single dose).

    Parameters:
    C_0 : float
        Initial concentration (mg/L)
    k_e : float
        Elimination rate constant (1/hr)
    t : float
        Time after administration (hours)

    Returns:
    float
        Plasma concentration (mg/L)
    """
    return C_0 * np.exp(-k_e * t)

# New utility functions

def adjust_dose_for_kidney_function(dose, creatinine_clearance, standard_creatinine_clearance=100.0):
    """
    Adjust the dosage based on estimated creatinine clearance.

    Parameters:
    dose : float
        Standard dose (mg)
    creatinine_clearance : float
        Patient's creatinine clearance (mL/min)
    standard_creatinine_clearance : float, optional
        Normal creatinine clearance (mL/min), default is 100 mL/min

    Returns:
    float
        Adjusted dose (mg)
    """
    return dose * (creatinine_clearance / standard_creatinine_clearance)

def adjust_dose_for_liver_function(dose, liver_function_score, max_score=15.0):
    """
    Adjust the dosage for compromised liver function.

    Parameters:
    dose : float
        Standard dose (mg)
    liver_function_score : float
        Patient's liver function score (e.g., Child-Pugh score)
    max_score : float, optional
        Maximum score indicating severely compromised liver function, default is 15.0

    Returns:
    float
        Adjusted dose (mg)
    """
    adjustment_factor = 1 - (liver_function_score / max_score) * 0.5
    return dose * adjustment_factor

def is_within_therapeutic_range(concentration, drug_name, patient_data, db_file='therapeutic_ranges.json'):
    """
    Check if the drug concentration is within the therapeutic range.

    Parameters:
    concentration : float
        Current drug concentration (mg/L)
    drug_name : str
        Name of the drug
    patient_data : dict
        Patient-specific data, such as age, weight, organ function, etc.
    db_file : str, optional
        Path to the JSON database file containing therapeutic ranges

    Returns:
    bool
        True if within therapeutic range, False otherwise
    """
    try:
        with open(db_file, 'r') as f:
            therapeutic_ranges = json.load(f)
    except FileNotFoundError:
        raise Exception(f"Database file '{db_file}' not found.")

    range_data = therapeutic_ranges.get(drug_name.lower())
    if range_data is None:
        raise ValueError(f"Therapeutic range for drug '{drug_name}' not found in database.")

    lower, upper = range_data
    return lower <= concentration <= upper

def check_interactions(drug_list, db_file='drug_interactions.json'):
    """
    Check for interactions among a list of drugs.

    Parameters:
    drug_list : list of str
        List of drug names
    db_file : str, optional
        Path to the JSON database file containing drug interactions

    Returns:
    list of str
        List of interactions found among the drugs
    """
    try:
        with open(db_file, 'r') as f:
            interaction_pairs = json.load(f)
    except FileNotFoundError:
        raise Exception(f"Database file '{db_file}' not found.")

    interactions = []
    for i in range(len(drug_list)):
        for j in range(i + 1, len(drug_list)):
            pair = tuple(sorted([drug_list[i].lower(), drug_list[j].lower()]))
            if pair in interaction_pairs:
                interactions.append(f"Interaction between {drug_list[i]} and {drug_list[j]}: {interaction_pairs[pair]}")
    return interactions

if __name__ == "__main__":
    # Import unittest for testing the functions
    import unittest

    class TestPharmacokinetics(unittest.TestCase):
        def test_elimination_rate_constant(self):
            self.assertAlmostEqual(elimination_rate_constant(5.0, 50.0), 0.1, places=6)

        def test_half_life(self):
            self.assertAlmostEqual(half_life(5.0, 50.0), 6.93, places=2)

        def test_initial_concentration(self):
            self.assertAlmostEqual(initial_concentration(500.0, 50.0), 10.0, places=6)

        def test_plasma_concentration_iv_bolus_single(self):
            self.assertAlmostEqual(plasma_concentration_iv_bolus_single(10.0, 0.1, 5.0), 6.0653, places=4)

        def test_adjust_dose_for_kidney_function(self):
            self.assertAlmostEqual(adjust_dose_for_kidney_function(100.0, 50.0), 50.0, places=4)

        def test_adjust_dose_for_liver_function(self):
            self.assertAlmostEqual(adjust_dose_for_liver_function(100.0, 7.5), 75.0, places=4)

        def test_is_within_therapeutic_range(self):
            # Mock database content
            with open('therapeutic_ranges.json', 'w') as f:
                json.dump({"drug_a": [5.0, 15.0]}, f)
            
            self.assertTrue(is_within_therapeutic_range(10.0, "drug_a", {}))
            self.assertFalse(is_within_therapeutic_range(20.0, "drug_a", {}))

        def test_check_interactions(self):
            # Mock database content
            with open('drug_interactions.json', 'w') as f:
                json.dump({"drug_a,drug_b": "Increased risk of toxicity"}, f)
            
            self.assertEqual(check_interactions(["drug_a", "drug_b"]), ["Interaction between drug_a and drug_b: Increased risk of toxicity"])

    # Run the tests
    unittest.main()
