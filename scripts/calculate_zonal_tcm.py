import numpy as np
import pygeoprocessing

inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Unpack scores & params from config
nodata_value = params.nodata

# Load the distance raster and population raster
# Load the distance raster and population raster
distance_array = pygeoprocessing.raster_to_numpy_array(inputs.distance_to_hq)
pop_array = pygeoprocessing.raster_to_numpy_array(inputs.population)

# Mask out nodata values
pop_array[pop_array == -9999.0] = np.nan

# Convert to km
distance_km = distance_array / 1000.0

# Distance bins
# "parameters have been estimated on the MENE outcomes for England 
# and then applied to the other EU LAUs through a transfer function approach. 
# The number of visits is then multiplied by the 52 weeks in a year."
dist_buffers = [1, 2, 3, 4]
k_vals = [0.0132, 0.0267, 0.0518, 0.1067]
alpha_vals = [0.00155, 0.00115, 0.00098, 0.00067]
travel_costs = [1.50, 3.00, 4.50, 6.00]  # € per trip for buffers 1 through 4


# Initialize output arrays
visits_array = np.zeros_like(pop_array, dtype=float)
value_array = np.zeros_like(pop_array, dtype=float)

for i in dist_buffers:
    k = k_vals[i - 1]
    alpha = alpha_vals[i - 1]
    travel_cost = travel_costs[i - 1]

    dist_buffer_mask = (distance_array == i)
    combined_mask = dist_buffer_mask & ~np.isnan(pop_array)

    if np.any(combined_mask):
        pop_vals = pop_array[combined_mask]
        visits = (1 + k) / (k + np.exp(-alpha * pop_vals))
        annual_visits = visits * 52
        visits_array[combined_mask] = annual_visits

        # Multiply visits by cost to get value
        value_array[combined_mask] = annual_visits * travel_cost

        print(f"Buffer {i}: Avg cost = €{travel_cost}, Avg visits = {visits.mean():.2f}, Avg value = {value_array[combined_mask].mean():.2f}")
    else:
        print(f"Buffer {i}: No valid population.")

# Ensure alignment with missing pop data
visits_array[np.isnan(pop_array)] = np.nan
value_array[np.isnan(pop_array)] = np.nan


# Save as Geotiff
ref_raster_info = pygeoprocessing.get_raster_info(inputs.population)
target_nodata = ref_raster_info['nodata'][0]
pixel_size = ref_raster_info['pixel_size']
origin = [ref_raster_info['geotransform'][0], ref_raster_info['geotransform'][3]]
projection = ref_raster_info['projection_wkt']

pygeoprocessing.numpy_array_to_raster(
    visits_array, 
    target_nodata,
    pixel_size, 
    origin,
    projection, 
    target_path=outputs.visits
)
pygeoprocessing.numpy_array_to_raster(
    value_array, 
    target_nodata,
    pixel_size, 
    origin,
    projection, 
    target_path=outputs.value
)

visits_array = np.zeros_like(pop_array, dtype=float)
value_array = np.zeros_like(pop_array, dtype=float)


# Next Steps: add country level variation in travel costs
# for i in dist_buffers:
#     k = k_vals[i - 1]
#     alpha = alpha_vals[i - 1]

#     for country_id in np.unique(country_array[~np.isnan(pop_array)]):
#         if np.isnan(country_id):
#             continue

#         # Get cost for this buffer + country
#         cost = travel_costs.get((int(country_id), i), None)
#         if cost is None:
#             continue  # Skip if no cost defined

#         # Create mask: distance == buffer & country == ID & pop is valid
#         mask = (
#             (distance_array == i) &
#             (country_array == country_id) &
#             ~np.isnan(pop_array)
#         )

#         if np.any(mask):
#             pop_vals = pop_array[mask]
#             visits = (1 + k) / (k + np.exp(-alpha * pop_vals))
#             annual_visits = visits * 52
#             visits_array[mask] = annual_visits
#             value_array[mask] = annual_visits * cost

#             print(f"Buffer {i}, Country {country_id}: €{cost}, mean value = {value_array[mask].mean():.2f}")

