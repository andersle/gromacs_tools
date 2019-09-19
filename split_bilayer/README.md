# Splitting a bilayer

This folder contains some scripts which can be used to split a bilayer.

## Usage

### Extract the upper/lower parts:

```bash
python extract_layers.py POPC_303K.gro POPC
```

### Join and translate two layers:

```bash
python translate.py upper_DPPC_293K.gro lower_POPC_303K.gro 1.5
```

### Remove some water molecules after solvation:

```bash
python remove_water.py with-water.gro
```

## Requirements

The following packages are required:

* matplotlib
* numpy
