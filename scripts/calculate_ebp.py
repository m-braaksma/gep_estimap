import numpy as np
import pygeoprocessing
from osgeo import gdal

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Unpack scores & params from config
lulc_scores = params.scores["lulc"]
pa_score = params.scores["pa_score"]
kba_score = params.scores["kba_score"]
ebp_classify = params.scores["ebp_classify"]
nodata_value = params.nodata

# Raster calculator function
def rec_potential_op(crop, forest, grass, othernat, urban, water, pa, kba):
    # Weighted LULC composite
    lulc_potential = (
        lulc_scores['cropland'] * crop +
        lulc_scores['forest']   * forest +
        lulc_scores['grassland']* grass +
        lulc_scores['othernat'] * othernat +
        lulc_scores['urban']    * urban +
        lulc_scores['water']    * water
    )
    # KBA mask & combined score
    kba_mask = np.where(kba != -9999, 1, 0)
    combined = lulc_potential + pa_score * pa + kba_score * kba_mask * kba

    # Classify into EBP bins (0–3)
    ebp_class = np.select(
        condlist=[
            combined > ebp_classify['high'],
            combined > ebp_classify['med_high'],
            combined > ebp_classify['med_low'],
            combined >= ebp_classify['low']
        ],
        choicelist=[0, 1, 2, 3],
        default=0
    )
    return ebp_class

# Prepare band list for calculator
base_list = [
    (inputs.cropland, 1),
    (inputs.forest, 1),
    (inputs.grass, 1),
    (inputs.othernat,1),
    (inputs.urban, 1),
    (inputs.water, 1),
    (inputs.pa, 1),
    (inputs.kba, 1),
]

# Run raster calculator
pygeoprocessing.raster_calculator(
    base_list,
    rec_potential_op,
    outputs.ebp,
    gdal.GDT_Int32,
    nodata_value,
    calc_raster_stats=True,
)