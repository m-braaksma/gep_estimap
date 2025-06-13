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

# Loop over each distance buffer (e.g., 1 km, 2 km, etc.)
for i in dist_buffers:
    # Get the corresponding model parameters for this buffer
    k = k_vals[i - 1]
    alpha = alpha_vals[i - 1]
    travel_cost = travel_costs[i - 1]

    # Create a mask for pixels that are exactly 'i' km from the HQ
    dist_buffer_mask = (distance_array == i)

    combined_mask = dist_buffer_mask & ~np.isnan(pop_array)

    if np.any(combined_mask):
        # Extract population values for these valid pixels
        pop_vals = pop_array[combined_mask]

        # Estimate weekly visits using a nonlinear model of population
        visits = (1 + k) / (k + np.exp(-alpha * pop_vals))

        # Scale up to annual visits (52 weeks per year)
        annual_visits = visits * 52

        # Store the calculated annual visits in the corresponding pixels
        # Each pixel's value represents the total number of recreation trips per year made by people living in that pixel
        # (to a generic recreation site, under the assumed buffer travel cost and access model)
        visits_array[combined_mask] = annual_visits

        # Estimate the economic value of those visits using travel cost
        # Each pixel represents the annual economic value of the visits made by the people in that pixel
        # (based on the assumed travel cost for that zone)
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


# Value per capita
# Replace only NaNs with np.nan; keep 0s for explicit handling
safe_pop_array = np.where(np.isnan(pop_array), np.nan, pop_array)
# Use np.where to handle the three cases:
# - If population is NaN → result is NaN
# - If population is 0 → result is 0
# - Else → value / population
value_per_capita_array = np.where(
    np.isnan(safe_pop_array), np.nan,
    np.where(safe_pop_array == 0, 0, value_array / safe_pop_array)
)
pygeoprocessing.numpy_array_to_raster(
    value_per_capita_array,
    target_nodata,
    pixel_size,
    origin,
    projection,
    target_path=outputs.value_per_capita
)


# Visits per km2
# Extract pixel dimensions in meters
# pixel_width, pixel_height = ref_raster_info['pixel_size']
# pixel_area_km2 = abs(pixel_width * pixel_height) / 1e6  # convert m² to km²
# visits_per_km2_array = visits_array / pixel_area_km2
# pygeoprocessing.numpy_array_to_raster(
#     visits_per_km2_array, 
#     target_nodata,
#     pixel_size, 
#     origin,
#     projection, 
#     target_path=outputs.visits_per_km2
# )