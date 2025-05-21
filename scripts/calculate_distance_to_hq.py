import numpy as np
import pygeoprocessing
from osgeo import gdal
from scipy.ndimage import distance_transform_edt

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Unpack scores & params from config
nodata_value = params.nodata

# Load high quality area mask
hq_array = pygeoprocessing.raster_to_numpy_array(inputs.hq)

# Create a distance map (distance to nearest 1)
distance = distance_transform_edt(hq_array == 0)

# Bin the distance map into 4 classes
bins = [1, 2, 3, 4]
buffered_array = np.digitize(distance, bins) + 1
buffered_array = buffered_array.astype(np.int16) 

# Save as raster
hq_info = pygeoprocessing.get_raster_info(inputs.hq)
pixel_size = hq_info['pixel_size']
origin = [hq_info['geotransform'][0], hq_info['geotransform'][3]]
projection = hq_info['projection_wkt']
pygeoprocessing.numpy_array_to_raster(
    buffered_array,
    nodata_value,
    pixel_size,
    origin,
    projection,
    target_path=outputs.distance_to_hq
)
