# Data

This directory contains the relevant data used for modelling in the serialized native R data forma `rds`:

```console
Alpine_Space_Boundary_4326.rds » Alpine Space outline (sf)
basins_compressed.rds          » Basins with static data (sf) -- polygons with attributes
df_falls.rds                   » Training data for fall-type (RF) model
df_flows.rds                   » Training data for flow-type (DF) model
df_slides.rds                  » Training data for slide-type (SL) model
myDFcv_5f_10r.rds              » Cross validation results for flow type (DF) model
myRFcv_5f_10r.rds              » Cross validation results for fall type (RF) model
mySLcv_5f_10r.rds              » Cross validation results for slide type (SL) model
varimpo_raw.rds                » Results (raw) from variable importance assessment (raw permute outputs)
varimpoDF.rds                  » Plot-ready variable importance data for flow-type (DF)
varimpoRF.rds                  » Plot-ready variable importance data for fall-type (RF)
varimpoSL.rds                  » Plot-ready variable importance data for slide-type (SL)
```
