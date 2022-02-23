# krovak05

Geodetic package for transformation ETRS89 (ETRF2000) coordinates to S-JTSK (Czech national coordinate system)
and heights to Bpv system (Baltic vertical datum After Adjustment).

## Installation

Run the following to install:

```python
pip install krovak05
```

## Methods

- _get_available_diff_tables()_ -> table_names - List[str]
- _interpolate_dydx(Y,X)_ -> dy, dx - float,float
- _interpolate_undulation(B,L)_ -> undulation - float
- _bicubic_dotr(Y,X)_ -> dy,dx - float, float
- _etrs_jtsk05(B,L,H)_ -> Y_jtsk05,X_jtsk05,H_bpv - float, float, float
- _etrs_jtsk(B,L,H)_ -> Y_jtsk,X_jtsk,H_bpv - float, float, float
- _jtsk05_jtsk(Y_jtsk05,X_jtsk05)_ -> Y_jtsk,X_jtsk - float, float
- _jtsk_jtsk05(Y_jtsk,X_jtsk)_ -> Y_jtsk05,X_jtsk05 - float, float
- _jtsk_etrs(self,Y,X,H)_ -> B,L,H - float, float

## Usage

```python
import krovak05

krovak = krovak05.Transformation()

## Undulation of kvasigeoid
undulation = krovak.interpolate_undulation(50, 15)
print(undulation)
# --> 44.438

## Differences between S-JTSK and S-JTSK/05
dy, dx = krovak.interpolate_dydx(750000, 1050000)
print(dy, dx)
# --> 0.072 -0.037

## Get list of possible dydx grid data
grids = krovak.get_available_diff_tables()
print(grids)
# --> ['table_yx_3_v1710', 'table_yx_3_v1202', 'table_yx_3_v1005']

## Transform ETRS89 (ETRF2000) coordinates to S-JTSK/05
B_etrs_in = 50
L_etrs_in = 15
H_etrs_in = 100

Y_sjtsk05, X_sjtsk05, H_bpv = krovak.etrs_jtsk05(
    B_etrs_in, L_etrs_in, H_etrs_in)
print(Y_sjtsk05, X_sjtsk05, H_bpv)
# --> 5703011.866856858 6058147.235673166 55.562

## Transform ETRS89 (ETRF2000) coordinates to S-JTSK
Y_sjtsk, X_sjtsk, H_bpv = krovak.etrs_jtsk(B_etrs_in, L_etrs_in, H_etrs_in)
print(Y_sjtsk, X_sjtsk, H_bpv)
# --> 703011.8997768582 1058147.294883166 55.562

## Reverse transformation S-JTSK coordinate to ETRS89
B_etrs_out, L_etrs_out, H_etrs_out = krovak.jtsk_etrs(
    Y_sjtsk, X_sjtsk, H_bpv)
print(B_etrs_out, L_etrs_out, H_etrs_out)
# --> 50.00000000579285 15.000000005855975 100.0

print("Differences:")
print(f"dB = {(B_etrs_in-B_etrs_out)*(math.pi/180)*6378000*1000} mm")
print(f"dL = {(B_etrs_in-B_etrs_out)*(math.pi/180)*(6378000*math.cos(B_etrs_in*(math.pi/180)))*1000} mm")
print(f"dH = {(H_etrs_in-H_etrs_out)*1000} mm")
# --> dB = -0.64 mm
# --> dL = -0.41 mm
# --> dH = 0.0 mm
```

### Set different grid table:

```python
krovak = krovak05.Transformation("table_yx_3_v1005")
```

# TODO

- Better documentation
- Rewrite to more time efficient code - constants, less calculation of goniometric functions

# Data validation

validation of data accuracy was performed using of the [CUZK transformation service](<https://geoportal.cuzk.cz/(S(idlg1tno0nodmoby14poaa1d))/Default.aspx?mode=TextMeta&text=wcts&menu=19>)

---

[Repository](https://github.com/SteveeH/krovak05)
