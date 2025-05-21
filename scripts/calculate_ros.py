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
def ros_op(human_input_class, ebp_class):
    human_input_class_index = human_input_class-1
    return ros_matrix[human_input_class_index, ebp_class].astype(np.int32)

base_raster_path_band_const_list = [
    (inputs.human_input, 1),
    (inputs.ebp, 1),
]

pygeoprocessing.raster_calculator(
    base_raster_path_band_const_list,
    ros_op,
    outputs.ros,
    gdal.GDT_Int32,
    nodata_value,
)
