# Bioinformatics Data Mining with Python APIs

This repository contains a complete Python workflow for mining biological and chemical data using public APIs. The project involves retrieving structural protein data, accessing metadata via UniProt and PubChem, and representing molecules with RDKit and Biopython.

## Project Structure

- `scripts/protein_data_mining.py`: Python script implementing the full pipeline
- `heteromolecules.sdf`: Output file with molecular data in SDF format
- `README.md`: Project documentation

## Objective

To extract and analyze biological and chemical data through programmatic access to public APIs, processing molecular data with Biopython and RDKit to build a usable chemical dataset.

---

## API Integration and Data Collection

- Downloaded `.cif` files for 10 PDB IDs using the RCSB PDB API, with error handling
- Queried UniProt’s API to retrieve metadata for the PDB ID `1tup`, including:
  - Entry publication and update dates
  - Review status (Swiss-Prot or Trembl)
  - Gene name, synonyms, organism
  - Protein name and sequence
  - Associated PDB IDs
- Extracted cofactor information from the UniProt entry
- Queried PubChem API to retrieve compound details:
  - CID, molecular weight, InChI, InChIKey, IUPAC name

---

## Molecular Data Processing

### Biopython

- Parsed the `4ogq.cif` file with `MMCIFParser` to extract heteromolecules (excluding water)
- Parsed the same file with `MMCIF2Dict` to extract and tabulate heteroatom names and 3-letter codes

### RDKit

- Queried PubChem for SMILES strings of each heteromolecule
- Converted SMILES to molecules, computed 2D coordinates
- Added molecular weight as a property
- Saved valid molecules to an `.sdf` file (`heteromolecules.sdf`)

---

## Output Example

The final `.sdf` file contains:
- Molecular structures from PubChem
- Standardized 2D coordinates
- Molecular weights stored as properties

All molecules were verified to be readable and complete using RDKit.

---

## Key Takeaways

- Demonstrated integration of **RESTful APIs** (RCSB, UniProt, PubChem)
- Applied **pandas** for structured data aggregation
- Used **Biopython** for structural data parsing
- Leveraged **RDKit** to convert and process molecular formats
- Automated end-to-end bioinformatics data retrieval and transformation

---

## Author

**Oscar Saiz Gutierrez**  
MSc in Bioinformatics  

---

**Note:** This project was developed as part of the *Python Programming* course in the MSc in Bioinformatics.
