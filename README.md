![Melodia](Melodia_logo.png)
# Melodia_py
## Protein Structure Analysis

**Melodia_py** is a Python library for computing Differential Geometry
and Knot Theory descriptors of protein structures. 

## Installation [Anaconda Python](https://www.anaconda.com/products/individual)

The first step is to clone Melodia's repository.
```shell
git clone https://github.com/rwmontalvao/Melodia.git
cd ./Melodia_py
```
We recommend using [Miniforge](https://github.com/conda-forge/miniforge) and Mamba for installation (optional). 
Miniforge is an Anaconda Python-compatible distribution with a faster and more reliable package manager (Mamba).
It is as simple to install as the Anaconda distribution. 

We start with creation of a new environment for **Melodia**.

```shell
conda env create -f environment.yml
```
or (optionally, but highly recommended)

```shell
mamba env create -f environment.yml
```
Next step is to activate the Melodia environment

```shell
conda activate melodia_py
```
To build and install **Melodia**:

```shell
python setup.py install
```

## Documentation
The *examples* folder contains Jupyter Notebooks, with short tutorials explaining **Melodia's** functionalities. 
1. Getting Started
2. Alignment Basics
3. Basic Similarity Analysis
4. Advanced Similarity Analysis
5. Machine Leaning Ensemble Analysis
6. Alignment Clustering and PDB Superimposition

### Authors
- Rinaldo W. Montalvão, PhD
- Antonio Marinho da Silva Neto, PhD
- William R. Pitt, PhD

### References
- Montalvão R, Smith R, Lovell S, Blundell T: CHORAL: a differential geometry approach to the prediction of the cores of protein structures. Bioinformatics. 2005, 21: 3719-3725.
- Chang PL, Rinne AW, Dewey TG: Structure alignment based on coding of local geometric measures. BMC Bioinformatics. 2006, 7:346.
- Leung H, Montaño B, Blundell T, Vendruscolo M, Montalvão R: ARABESQUE: A tool for protein structural comparison using differential geometry and knot theory. World Res J Peptide Protein. 2012, 1: 33-40.
- Pitt WR, Montalvão R, Blundell T: Polyphony: superposition independent methods for ensemble-based drug discovery. BMC Bioinformatics. 2014, 15:324 
- Marinho da Silva Neto A, Reghim Silva S, Vendruscolo M, Camilloni C, Montalvão R: A Superposition Free Method for Protein Conformational Ensemble Analyses and Local Clustering Based on a Differential Geometry Representation of Backbone. Proteins: Structure, Function, and Bioinformatics. 2018, 87(4):302-312
- Marinho da Silva Neto A, Montalvão R, Gondim Martins DB, Lima Filho JL, Madeiros Castelletti CH: A model of key residues interactions for HPVs E1 DNA binding domain-DNA interface based on HPVs residues conservation profiles and molecular dynamics simulations, Journal of Biomolecular Structure and Dynamics. 2019, 38(12):3720-3729.
