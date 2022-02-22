# krovak05

Geodetic package for transformation ETRS89 (ETRF2000) coordinates to S-JTSK (Czech national coordinate system)
and heights to Bpv system (Baltic vertical datum After Adjustment).

## Installation

Run the following to install:

```python
pip install krovak05 
```

## Usage

```python
krovak = Krovak05()

output = krovak.etrs_jtsk(50.000,15.000,100)

--> (703011.898, 1058147.296, 55.562)

```

[Repository](https://github.com/SteveeH/krovak05)
