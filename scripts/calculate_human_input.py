import numpy as np
import pygeoprocessing
from osgeo import gdal

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Unpack scores & params from config
urban_classify = params.scores["urban_classify"]
road_classify = params.scores["road_classify"]
urban_road_matrix = np.array(params.scores['urban_road_matrix'])
ebp_classify = params.scores["ebp_classify"]
nodata_value = params.nodata

# Raster calculator function
def human_input_op(urban, roads):
    urban_class = np.select(
        condlist=[
            urban > urban_classify['high'],
            urban > urban_classify['med_high'],
            urban > urban_classify['med'],
            urban > urban_classify['med_low'],
            urban >= urban_classify['low']
        ],
        choicelist=[0, 1, 2, 3, 4],  # 0 = most urban, 4 = least
        default=4
    )

    road_class = np.select(
        condlist=[
            roads > road_classify['high'],
            roads > road_classify['med_high'],
            roads > road_classify['med'],
            roads > road_classify['med_low'],
            roads >= road_classify['low']
        ],
        choicelist=[0, 1, 2, 3, 4],  # 0 = closest to roads, 4 = farthest
        default=4
    )

    return urban_road_matrix[urban_class, road_class].astype(np.int32)

# Prepare band list for calculator
base_list = [
    (inputs.urban, 1),
    (inputs.road, 1),
]

print("Urban input:", inputs.urban)
print("Road input:", inputs.road)
# Run raster calculator
pygeoprocessing.raster_calculator(
    base_list,
    human_input_op,
    outputs.human_input,
    gdal.GDT_Int32,
    nodata_value,
)