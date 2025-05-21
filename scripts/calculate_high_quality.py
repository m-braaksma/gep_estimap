import numpy as np
import pygeoprocessing
from osgeo import gdal

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Unpack scores & params from config
ros_matrix = np.array(params.scores['ros_matrix'])
nodata_value = params.nodata

# Raster calculator function
def hq_op(ros_class):
    return np.where(ros_class == 9, 1, 0).astype(np.int32)

base_raster_path_band_const_list = [
    (inputs.ros, 1),
]

pygeoprocessing.raster_calculator(
    base_raster_path_band_const_list,
    hq_op,
    outputs.hq,
    gdal.GDT_Int32,
    nodata_value,
)
