# krovak05

Geodetic package for transformation ETRS89 (ETRF2000) coordinates to S-JTSK (Czech national coordinate system)
and heights to Bpv system (Baltic vertical datum After Adjustment).

## Installation

<code> 
pip -m install krovak05 
</code>

## Example of use:

<code> 
krovak = Krovak05()

output = krovak.etrs_jtsk(50.000,15.000,100)

--> (703011.898, 1058147.296, 55.562)

</code>

[Repository](https://github.com/SteveeH/krovak05)
