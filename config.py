"""
Global Configuration & Reference Rules for Population Pharmacogenomics Analyzer.
Contains reference mapping dictionaries for Indian geographic regions,
SNP valid allele definitions, CYP2C19 star allele rules, diplotype calls, and CPIC phenotypes.
"""

from typing import Dict, Set, List

# -----------------------------------------------------------------------------
# 1. GEOGRAPHIC & DEMOGRAPHIC MAPPINGS (INDIA)
# -----------------------------------------------------------------------------
STATE_TO_REGION: Dict[str, str] = {
    # South India
    'Andhra Pradesh': 'South India',
    'Telangana': 'South India',
    'Karnataka': 'South India',
    'Tamil Nadu': 'South India',
    'Kerala': 'South India',
    'Puducherry': 'South India',
    'Goa': 'South India',
    'Lakshadweep': 'South India',
    'Andaman and Nicobar Islands': 'South India',
    
    # North India
    'Punjab': 'North India',
    'Haryana': 'North India',
    'Uttar Pradesh': 'North India',
    'Uttarakhand': 'North India',
    'Himachal Pradesh': 'North India',
    'Delhi': 'North India',
    'Chandigarh': 'North India',
    'Jammu & Kashmir': 'North India',
    'Jammu and Kashmir': 'North India',
    'Ladakh': 'North India',
    'Rajasthan': 'North India',
    
    # East India
    'Odisha': 'East India',
    'West Bengal': 'East India',
    'Bihar': 'East India',
    'Jharkhand': 'East India',
    'Assam': 'East India',
    'Meghalaya': 'East India',
    'Tripura': 'East India',
    'Arunachal Pradesh': 'East India',
    'Manipur': 'East India',
    'Mizoram': 'East India',
    'Nagaland': 'East India',
    'Sikkim': 'East India',
    
    # West India
    'Maharashtra': 'West India',
    'Gujarat': 'West India',
    'Dadra and Nagar Haveli': 'West India',
    'Daman and Diu': 'West India',
    
    # Central India
    'Madhya Pradesh': 'Central India',
    'Chhattisgarh': 'Central India'
}

# City/District to State gazetteer for native place prediction when State is missing/ambiguous
CITY_TO_STATE: Dict[str, str] = {
    'west godavari': 'Andhra Pradesh',
    'east godavari': 'Andhra Pradesh',
    'visakhapatnam': 'Andhra Pradesh',
    'guntur': 'Andhra Pradesh',
    'vijayawada': 'Andhra Pradesh',
    'hyderabad': 'Telangana',
    'warangal': 'Telangana',
    'bangalore': 'Karnataka',
    'bengaluru': 'Karnataka',
    'mysore': 'Karnataka',
    'chennai': 'Tamil Nadu',
    'coimbatore': 'Tamil Nadu',
    'madurai': 'Tamil Nadu',
    'trivandrum': 'Kerala',
    'thiruvananthapuram': 'Kerala',
    'kochi': 'Kerala',
    'lucknow': 'Uttar Pradesh',
    'kanpur': 'Uttar Pradesh',
    'varanasi': 'Uttar Pradesh',
    'agra': 'Uttar Pradesh',
    'noida': 'Uttar Pradesh',
    'ludhiana': 'Punjab',
    'amritsar': 'Punjab',
    'jalandhar': 'Punjab',
    'jaipur': 'Rajasthan',
    'jodhpur': 'Rajasthan',
    'pune': 'Maharashtra',
    'mumbai': 'Maharashtra',
    'nagpur': 'Maharashtra',
    'nashik': 'Maharashtra',
    'ahmedabad': 'Gujarat',
    'surat': 'Gujarat',
    'vadodara': 'Gujarat',
    'cuttack': 'Odisha',
    'bhubaneswar': 'Odisha',
    'puri': 'Odisha',
    'kolkata': 'West Bengal',
    'howrah': 'West Bengal',
    'patna': 'Bihar',
    'ranchi': 'Jharkhand',
    'guwahati': 'Assam',
    'shillong': 'Meghalaya',
    'bhopal': 'Madhya Pradesh',
    'indore': 'Madhya Pradesh'
}

# -----------------------------------------------------------------------------
# 2. SNP VALIDATION RULES & ALLELE DEFINITIONS
# -----------------------------------------------------------------------------
SNP_CONFIG: Dict[str, Dict] = {
    'CYP2C19*2': {
        'column_patterns': ['CYP2C19*2', 'rs4244285', 'cyp2c19_2'],
        'rsid': 'rs4244285',
        'ref_allele': 'G',
        'var_allele': 'A',
        'valid_genotypes': {'GG', 'GA', 'AA'},
        'wildtype_genotype': 'GG',
        'het_genotype': 'GA',
        'hom_var_genotype': 'AA'
    },
    'CYP2C19*3': {
        'column_patterns': ['CYP2C19*3', 'rs4986893', 'cyp2c19_3'],
        'rsid': 'rs4986893',
        'ref_allele': 'G',
        'var_allele': 'A',
        'valid_genotypes': {'GG'},
        'wildtype_genotype': 'GG',
        'het_genotype': 'GA',
        'hom_var_genotype': 'AA'
    },
    'CYP2C19*17': {
        'column_patterns': ['CYP2C19*17', 'rs12248560', 'cyp2c19_17'],
        'rsid': 'rs12248560',
        'ref_allele': 'C',
        'var_allele': 'T',
        'valid_genotypes': {'CC', 'CT', 'TT'},
        'wildtype_genotype': 'CC',
        'het_genotype': 'CT',
        'hom_var_genotype': 'TT'
    }
}

# -----------------------------------------------------------------------------
# 3. CPIC CYP2C19 DIPLOTYPE TO PHENOTYPE MAPPING MATRIX
# -----------------------------------------------------------------------------
DIPLOTYPE_PHENOTYPE_MAP: Dict[str, str] = {
    '*1/*1': 'Normal Metabolizer',
    '*1/*17': 'Rapid Metabolizer',
    '*17/*17': 'Ultrarapid Metabolizer',
    '*1/*2': 'Intermediate Metabolizer',
    '*1/*3': 'Intermediate Metabolizer',
    '*17/*2': 'Indeterminate',
    '*2/*17': 'Indeterminate',
    '*3/*17': 'Indeterminate',
    '*17/*3': 'Indeterminate',
    '*2/*2': 'Poor Metabolizer',
    '*2/*3': 'Poor Metabolizer',
    '*3/*3': 'Poor Metabolizer',
    'Indeterminate': 'Indeterminate'
}

PHENOTYPE_ORDER: List[str] = [
    'Intermediate Metabolizer',
    'Normal Metabolizer',
    'Indeterminate',
    'Poor Metabolizer',
    'Rapid Metabolizer',
    'Ultrarapid Metabolizer'
]

ALPHA_SIGNIFICANCE: float = 0.05
