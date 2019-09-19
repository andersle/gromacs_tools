# Area per lipid.

Some scripts to obtain the area per lipid for bilayers.
Note: The box we consider here is assumed to be **rectangular**.

## Usage ``get_area.sh``

```
./get_area.sh
```

Will extract the box lengths from all ``.edr`` files found.


## Usage ``get_area.py``

```
get_area.py run*/box-*.xvg
```

will read the box lengths from the given ``.xvg`` files and calculate
an area per lipid assuming:

* The area is given by the box lengths in ``x`` and ``y`` directions.
* A certain number of lipid molecules as hard-coded in ``get_area.py``.
