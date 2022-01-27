# Melodia
## Differential Geometry of Protein Backbones

**Melodia** is a python based library for computing Differential Geometry
of protein structures (backbones). 

## Installation
It is recommended the creation of a new environment for **Melodia**.

```shell
conda create --name melodia python=3.8
conda activate melodia
```

**Melodia** requires [Jupyterlab](https://jupyter.org/) and [NGL Viewer](https://github.com/nglviewer/nglview).

```shell
conda install jupyterlab
conda install nglview -c conda-forge
jupyter-nbextension enable nglview --py --sys-prefix
```

The examples need [seaborn](https://seaborn.pydata.org/) and [scikit-learn](https://scikit-learn.org/stable).

```shell
conda install -c conda-forge seaborn
conda install -c conda-forge scikit-learn
```

To build and install **Melodia**:

```shell
python setup.py install
```
### Authors
- Rinaldo W. Montalvão, PhD
- William R. Pitt, PhD

### References
- Montalvão R, Smith R, Lovell S, Blundell T: CHORAL: a differential geometry approach to the prediction of the cores of protein structures. Bioinformatics. 2005, 21: 3719-3725.
- Leung H, Montaño B, Blundell T, Vendruscolo M, Montalvão R: ARABESQUE: A tool for protein structural comparison using differential geometry and knot theory. World Res J Peptide Protein. 2012, 1: 33-40.
- Pitt WR, Montalvão R, Blundell T: Polyphony: superposition independent methods for ensemble-based drug discovery. BMC Bioinformatics. 2014, 15:324 
- Marinho da Silva Neto A, Reghim Silva S, Vendruscolo M, Camilloni C, Montalvão R: A Superposition Free Method for Protein Conformational Ensemble Analyses and Local Clustering Based on a Differential Geometry Representation of Backbone. Proteins: Structure, Function, and Bioinformatics. 2018, 87(4):302-312
- Marinho da Silva Neto A, Montalvão R, Gondim Martins DB, Lima Filho JL, Madeiros Castelletti CH: A model of key residues interactions for HPVs E1 DNA binding domain-DNA interface based on HPVs residues conservation profiles and molecular dynamics simulations, Journal of Biomolecular Structure and Dynamics. 2019, 38(12):3720-3729.
