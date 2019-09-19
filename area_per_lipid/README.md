# Area per lipid.

Some scripts to help obtain the area per lipid for bilayers.

## Usage ``get_area.sh``

```
./get_area.sh
```

Will extract the box lengths from all .edr files found.


## Usage ``get_area.py`

```
get_area.py run*/box-*.xvg
```

will read the box lengths from the given .xvg files and calculate
an area per lipid assuming:

* The area is given by the box length in ``x`` and ``y`` directions.
* A number of lipid molecules as hard-coded in ``get_area.py``.
