from sfincsOutputLib import sfincsRadialAndErScan


ds = sfincsRadialAndErScan('/u/lebra/src/stelloptPlusSfincs/outsideTest/sixthObj', verbose=0)

ErQ=ds.Ersearch(ErQuantity='Er', verbose=1, launch='no')
