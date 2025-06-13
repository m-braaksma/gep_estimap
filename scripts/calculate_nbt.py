import numpy as np
import pandas as pd
import geopandas as gpd
import pygeoprocessing
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')
from osgeo import gdal
import tempfile
import os
import re


inputs = snakemake.input
outputs = snakemake.output
params = snakemake.params
config = snakemake.config

# Parameters
nodata_value = params.nodata
target_year = params.target_year
tourism_type = params.tourism_type


def load_admin_data(admin_path):
    """Load and filter admin boundaries to ADM2 level"""
    print("Loading admin boundaries...")
    admin_gdf = gpd.read_file(admin_path, layer='ADM_2', engine='pyogrio')
    
    # Ensure we have necessary columns (adjust based on your data structure)
    required_cols = ['geometry']
    if 'ISO_A3' in admin_gdf.columns:
        required_cols.append('ISO_A3')
    if 'NAME_2' in admin_gdf.columns:
        required_cols.append('NAME_2')
    
    return admin_gdf[required_cols + [col for col in admin_gdf.columns if col not in required_cols]]


def extract_clean_overnights(df, tourism_type):
    """Extract and clean overnight stay data from UNWTO format"""
    df = df.copy()
    df.columns = df.columns.map(str)
    
    country_col = "Basic data and indicators"
    indicator_col = "Unnamed: 6"
    type_col = "Unnamed: 5"
    
    # Forward fill the hierarchical columns
    df[country_col] = df[country_col].ffill()
    df[indicator_col] = df[indicator_col].ffill()
    df[type_col] = df[type_col].ffill()
    
    # Filter only overnights rows
    overnight_df = df[df[indicator_col].str.contains("Overnights", na=False)]
    
    # Keep only relevant columns
    year_cols = [col for col in df.columns if col.isdigit()]
    overnight_df = overnight_df[[country_col, type_col] + year_cols]
    
    # Melt to long format
    tidy = overnight_df.melt(
        id_vars=[country_col, type_col],
        var_name="year",
        value_name="overnights"
    )
    
    tidy.rename(columns={country_col: "country_name", type_col: "overnight_type"}, inplace=True)
    tidy["tourism_type"] = tourism_type
    tidy["year"] = tidy["year"].astype(int)
    tidy["overnights"] = pd.to_numeric(tidy["overnights"], errors="coerce")
    
    return tidy


def load_overnight_data(overnight_path_domestic, overnight_path_international):
    """Load and process UNWTO overnight stay data from both domestic and international files"""
    print("Loading overnight stay data...")
    
    # Load both datasets
    # Find starting row (where 'AFGHANISTAN' appears)
    domestic_start = pd.read_excel(overnight_path_domestic, sheet_name='Domestic Tourism-Accommodation')
    inbound_start = pd.read_excel(overnight_path_international, sheet_name='Inbound Tourism-Accommodation')
    domestic_start_idx = domestic_start[domestic_start.iloc[:, 3] == 'AFGHANISTAN'].index[0]
    inbound_start_idx = inbound_start[inbound_start.iloc[:, 3] == 'AFGHANISTAN'].index[0]
    # Load the data from the correct starting row
    domestic_data = pd.read_excel(overnight_path_domestic, sheet_name='Domestic Tourism-Accommodation', skiprows=domestic_start_idx)
    international_data = pd.read_excel(overnight_path_international, sheet_name='Inbound Tourism-Accommodation', skiprows=inbound_start_idx)
    
    # Extract and clean data from both files
    domestic_df = extract_clean_overnights(domestic_data, "domestic")
    international_df = extract_clean_overnights(international_data, "international")
    
    # Combine datasets
    panel_df = pd.concat([domestic_df, international_df], ignore_index=True)
    
    # Normalize types (in case of spacing or case issues)
    panel_df["overnight_type"] = panel_df["overnight_type"].str.lower().str.strip()
    
    # Pivot to wide format
    pivoted = panel_df.pivot_table(
        index=["country_name", "year"],
        columns=["overnight_type", "tourism_type"],
        values="overnights",
        aggfunc="first"
    ).reset_index()
    
    # Flatten MultiIndex columns
    pivoted.columns = ['country_name', 'year'] + [
        f"{otype}_overnights_{ttype}" for otype, ttype in pivoted.columns[2:]
    ]
    
    # Add combined (total + hotel) columns
    for ttype in ['domestic', 'international']:
        total_col = f"total_overnights_{ttype}"
        hotel_col = f"hotels and similar establishments_overnights_{ttype}"
        combined_col = f"combined_overnights_{ttype}"
        
        if total_col in pivoted.columns and hotel_col in pivoted.columns:
            pivoted[combined_col] = pivoted[total_col] + pivoted[hotel_col]
    
    # Rename hotel columns for clarity
    pivoted.rename(columns={
        'hotels and similar establishments_overnights_domestic': 'hotel_overnights_domestic',
        'hotels and similar establishments_overnights_international': 'hotel_overnights_international'
    }, inplace=True)
    
    # Create final column structure
    final_cols = ['country_name', 'year']
    for ttype in ['domestic', 'international']:
        final_cols += [
            f"total_overnights_{ttype}",
            f"hotel_overnights_{ttype}",
            f"combined_overnights_{ttype}"
        ]
    
    # Filter to final cols that exist
    final_cols = [col for col in final_cols if col in pivoted.columns]
    final_df = pivoted[final_cols]
    
    return final_df


def create_country_nights_mapping(overnight_df, target_year, tourism_type):
    """Create mapping from country to overnight stays for target year and tourism type"""
    # Filter to target year
    year_data = overnight_df[overnight_df['year'] == target_year].copy()
    
    if year_data.empty:
        print(f"Warning: No data found for year {target_year}")
        return {}
    
    # Select appropriate column based on tourism type
    if tourism_type == "international":
        # Try different column patterns for international data
        possible_cols = [
            'total_overnights_international',
            'hotel_overnights_international', 
            'combined_overnights_international'
        ]
    elif tourism_type == "domestic":
        # Try different column patterns for domestic data
        possible_cols = [
            'total_overnights_domestic',
            'hotel_overnights_domestic',
            'combined_overnights_domestic'
        ]
    else:  # "total" or combined
        # Try to combine both types or use total
        possible_cols = [
            'combined_overnights_domestic',
            'combined_overnights_international',
            'total_overnights_domestic',
            'total_overnights_international'
        ]
    
    # Find the first available column
    overnight_col = None
    for col in possible_cols:
        if col in year_data.columns:
            overnight_col = col
            break
    
    if overnight_col is None:
        print(f"Warning: No suitable overnight column found for tourism type '{tourism_type}'")
        print(f"Available columns: {list(year_data.columns)}")
        return {}
    
    # Create country mapping (assuming we need to map country names to ISO codes)
    # This is a simplified mapping - you may need a more comprehensive country name to ISO mapping
    country_nights = {}
    for _, row in year_data.iterrows():
        if pd.notna(row[overnight_col]):
            # Use country name as key for now - you may need to map to ISO_A3 codes
            country_nights[row['country_name']] = row[overnight_col]
    
    print(f"Loaded overnight data for {len(country_nights)} countries using column '{overnight_col}'")
    return country_nights


def dms_to_decimal(value):
    """Function to convert DMS and degrees-decimal minutes to decimal degrees."""

    # Clean the value by removing trailing periods
    value = value.rstrip('.')  # Remove any trailing periods

    # Check if the input is in DMS format (with optional spaces)
    dms_pattern = r'([-+]?\d+)°\s*(\d+)?\'?\s*([\d.]+)?\"?'
    match_dms = re.match(dms_pattern, value)
    
    # Check if the input is in degrees and decimal minutes format
    degrees_minutes_pattern = r'([-+]?\d+)\s+([\d.]+)'
    match_dm = re.match(degrees_minutes_pattern, value)
    
    # Check if the input is in "deg" format
    deg_pattern = r'([-+]?\d+)\s*deg\s*(\d+)?\.?(\d+)?'
    match_deg = re.match(deg_pattern, value)
    
    if match_dms:
        degrees = float(match_dms.group(1))
        minutes = float(match_dms.group(2)) if match_dms.group(2) else 0  # Default to 0 if not present
        seconds = float(match_dms.group(3)) if match_dms.group(3) else 0  # Default to 0 if not present
        decimal = degrees + (minutes / 60) + (seconds / 3600)
        return decimal
    elif match_dm:
        degrees = float(match_dm.group(1))
        minutes = float(match_dm.group(2))
        decimal = degrees + (minutes / 60)
        return decimal
    elif match_deg:
        degrees = float(match_deg.group(1))
        minutes = float(match_deg.group(2)) if match_deg.group(2) else 0  # Default to 0 if not present
        seconds = float(match_deg.group(3)) if match_deg.group(3) else 0  # Default to 0 if not present
        decimal = degrees + (minutes / 60) + (seconds / 3600)
        return decimal
    else:
        # If not DMS or degrees-decimal minutes, check for decimal format with comma
        value = value.replace(',', '.')  # Replace comma with period
        try:
            return float(value)  # Try converting to float
        except ValueError:
            return None  # Return None if conversion fails

def allocate_nights_to_admin(admin_gdf, accommodation_gdf, country_nights):
    """Allocate national overnight stays to admin regions based on accommodation locations."""
    print("Allocating nights to admin regions based on accommodation distribution...")

    admin_gdf = admin_gdf.copy()
    admin_gdf['nights'] = 0.0

    # Ensure both GeoDataFrames use the same CRS
    if accommodation_gdf.crs != admin_gdf.crs:
        accommodation_gdf = accommodation_gdf.to_crs(admin_gdf.crs)

    # Spatial join: assign each accommodation to an admin region
    joined = gpd.sjoin(accommodation_gdf, admin_gdf, how='left', predicate='within')

    # Count accommodations per admin region
    accom_counts = joined.groupby(admin_gdf.index.name or 'index_right').size()
    admin_gdf['accom_count'] = admin_gdf.index.map(accom_counts).fillna(0)

    # Loop through each country and allocate based on accommodation counts
    for country_key, nights in country_nights.items():
        # Match country using name
        country_mask = admin_gdf['COUNTRY'].str.upper() == country_key
        # name_mask = admin_gdf['COUNTRY'].str.upper() == country_key
        # name_mask = admin_gdf['NAME_0'].str.contains(str(country_key), case=False, na=False) if 'NAME_0' in admin_gdf.columns else pd.Series(False, index=admin_gdf.index)

        # country_mask = iso_mask | name_mask
        country_admin = admin_gdf[country_mask]

        total_accom = country_admin['accom_count'].sum()
        if total_accom > 0:
            proportions = country_admin['accom_count'] / total_accom
            allocated = proportions * nights
            admin_gdf.loc[country_mask, 'nights'] = allocated.values
        else:
            print(f"No accommodations found for country {country_key}. Skipping allocation.")

    # Drop intermediate accommodation count column
    admin_gdf.drop(columns="accom_count", inplace=True)

    return admin_gdf



def calculate_nature_attribution(distance_raster_path, accommodation_gdf, admin_gdf):
    """Calculate nature-based tourism attribution for accommodation and admin regions."""
    print("Calculating nature attribution...")
    
    # Load distance raster
    distance_array = pygeoprocessing.raster_to_numpy_array(distance_raster_path)
    raster_info = pygeoprocessing.get_raster_info(distance_raster_path)
    
    # Reproject accommodation GeoDataFrame to match raster if needed
    if raster_info['projection_wkt']:
        accommodation_gdf = accommodation_gdf.to_crs(raster_info['projection_wkt'])
    
    # Initialize a list to store distance scores for each accommodation
    distance_scores = []
    
    # Define a mapping from raster values to distance scores
    distance_score_mapping = {
        1: 1.0,     # 0-1 km
        2: 0.75,    # 1-2 km
        3: 0.5,     # 2-3 km
        4: 0.25,    # 3-4 km
        5: 0.0      # >4 km
    }
    
    # Sample distance values at accommodation locations
    for point in accommodation_gdf.geometry:
        try:
            # Convert point to raster coordinates
            x, y = point.x, point.y
            col = int((x - raster_info['geotransform'][0]) / raster_info['geotransform'][1])
            row = int((y - raster_info['geotransform'][3]) / raster_info['geotransform'][5])
            
            # Check bounds
            if 0 <= row < distance_array.shape[0] and 0 <= col < distance_array.shape[1]:
                raster_value = distance_array[row, col]
                # Assign score based on raster value
                score = distance_score_mapping.get(raster_value, 0.0)  # Default to 0 if not found
                distance_scores.append(score)
            else:
                distance_scores.append(np.nan)
        except Exception as e:
            print(f"Error processing point {point}: {e}")
            distance_scores.append(np.nan)
    
    # Add distance scores to accommodation GeoDataFrame
    accommodation_gdf['distance_score'] = distance_scores
    
    # Perform spatial join to associate accommodations with admin regions
    accommodation_with_admin = gpd.sjoin(accommodation_gdf, admin_gdf, how='left', predicate='intersects')
    print(accommodation_with_admin.columns)

    # Calculate average distance score per admin area
    admin_scores = accommodation_with_admin.groupby('GID_2')['distance_score'].mean().to_dict()
    
    # Calculate NBT overnights for each admin region
    admin_gdf['nbt_overnights'] = admin_gdf.index.map(admin_scores).fillna(0.0)  # Default to 0.0 if no hotels
    
    return admin_gdf, accommodation_gdf


# Main processing
def main():
    # Load data
    admin_gdf = load_admin_data(inputs.admin_borders)
    
    # Load overnight data from both domestic and international files
    # Assuming inputs now has both domestic and international paths
    if hasattr(inputs, 'overnight_stays_domestic') and hasattr(inputs, 'overnight_stays_international'):
        overnight_df = load_overnight_data(inputs.overnight_stays_domestic, inputs.overnight_stays_international)
        country_nights = create_country_nights_mapping(overnight_df, target_year, tourism_type)
    else:
        # Fallback to original single file approach
        overnight_df = load_overnight_data(inputs.overnight_stays, inputs.overnight_stays)
        country_nights = create_country_nights_mapping(overnight_df, target_year, tourism_type)
    
    # Load accommodation data
    accommodation_df = pd.read_excel(inputs.accommodation_locations)

    # Convert DMS to decimal for 'lat' and 'lon' columns
    accommodation_df['GEO_LONG'] = accommodation_df['GEO_LONG'].apply(dms_to_decimal)
    accommodation_df['GEO_LAT'] = accommodation_df['GEO_LAT'].apply(dms_to_decimal)
    accommodation_gdf = gpd.GeoDataFrame(
        accommodation_df,
        geometry=gpd.points_from_xy(accommodation_df.GEO_LONG, accommodation_df.GEO_LAT),
        crs='EPSG:4326'
    )

    # Allocate national nights to admin regions
    admin_gdf = allocate_nights_to_admin(admin_gdf, accommodation_gdf, country_nights)
    
    # Calculate nature attribution
    admin_gdf, accommodation_gdf = calculate_nature_attribution(
        inputs.distance_to_hq, accommodation_gdf, admin_gdf
    )
    
    # Save accommodation zones for diagnostics
    admin_gdf.to_file(outputs.admin_gdf, driver='GPKG')
    accommodation_gdf.to_file(outputs.accommodation_zones, driver='GPKG')
    
    # Print summary statistics
    print(f"\nNBT Summary Statistics:")
    print(f"Total admin regions: {len(admin_gdf)}")
    print(f"Total accommodation units: {len(accommodation_gdf)}")
    print(f"Total nights allocated: {admin_gdf['nights'].sum():,.0f}")
    # print(f"Total nature-based nights: {admin_gdf['nature_nights'].sum():,.0f}")
    # print(f"Average nature attribution: {admin_gdf['nature_attribution'].mean():.2%}")
    # print(f"Estimated total value: €{value_array.sum():,.0f}")


if __name__ == "__main__":
    main()