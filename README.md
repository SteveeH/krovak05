# krovak05

Geodetic package for transformation ETRS89 (ETRF2000) coordinates to S-JTSK (Czech national coordinate system)
and heights to Bpv system (Baltic vertical datum After Adjustment).

## Installation

Run the following to install:

```python
pip install krovak05 
```

## Methods

- *get_available_diff_tables()* -> table_names  - List[str]
- *interpolate_dydx(Y,X)* -> dy, dx - float,float
- *interpolate_undulation(B,L)* -> undulation - float
- *bicubic_dotr(Y,X)* -> dy,dx - float, float
- *etrs_jtsk05(B,L,H)* -> Y_jtsk05,X_jtsk05,H_bpv - float, float, float
- *etrs_jtsk(B,L,H)* -> Y_jtsk,X_jtsk,H_bpv - float, float, float
- *jtsk05_jtsk(Y_jtsk05,X_jtsk05)* -> Y_jtsk,X_jtsk - float, float
- *jtsk_jtsk05(Y_jtsk,X_jtsk)* -> Y_jtsk05,X_jtsk05 - float, float

---
<span style="color:red">
jtsk_etrs(self,Y,X,H) -> raise Exception (!!Not Implemented yet!!)
</span>


## Usage

```python

krovak = Transformation()

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
Y_sjtsk05, X_sjtsk05, H_bpv = krovak.etrs_jtsk05(50, 15, 100)
print(Y_sjtsk05, X_sjtsk05, H_bpv)
# --> 5703011.866856858 6058147.235673166 55.562

## Transform ETRS89 (ETRF2000) coordinates to S-JTSK
Y_sjtsk, X_sjtsk, H_bpv = krovak.etrs_jtsk(50, 15, 100)
print(Y_sjtsk, X_sjtsk, H_bpv)
# --> 703011.8997768582 1058147.294883166 55.562
```

### Set different grid table:

```python
krovak = Transformation("table_yx_3_v1005")
```

# TODO
- Better documentation
- Code reverse transformation from S-JTSK to ETRS89

# Data validation
validation of data accuracy was performed using of the [ČUZK transformation service](https://geoportal.cuzk.cz/(S(idlg1tno0nodmoby14poaa1d))/Default.aspx?mode=TextMeta&text=wcts&menu=19)


---
[Repository](https://github.com/SteveeH/krovak05)
