from pathlib import Path

# Load config
configfile: "config.yaml"

# Get data_dir from config
data_dir = config["data_dir"]
USE_S3 = str(data_dir).startswith("s3://") or str(data_dir).startswith("s3:")
if USE_S3:
    S3 = S3RemoteProvider()
    def D(path_str):
        return S3(f"{data_dir.rstrip('/')}/{path_str}")
else:
    DATA_DIR = Path(data_dir).expanduser()
    def D(path_str):
        return str(DATA_DIR / path_str)

rule all:
    input:
        # ESTIMAP outputs
        D(config["estimap"]["ebp"]),
        D(config["estimap"]["human_inputs"]),
        D(config["estimap"]["ros"]),
        D(config["estimap"]["high_quality"]),
        D(config["estimap"]["distance_to_hq"]),
        # Daily Rec ouputs
        D(config["daily_rec"]["visits"]),
        D(config["daily_rec"]["value"]),
        D(config["daily_rec"]["value_per_capita"]),
        # NBT outputs
        D(config["nbt"]["nights"]),
        # D(config["nbt"]["value"]),
        # D(config["nbt"]["nights_per_capita"]),
        D(config["nbt"]["accommodation_zones"])

rule calculate_ebp:
    input:
        cropland = D(config["inputs"]["lulc"]["cropland"]),
        forest = D(config["inputs"]["lulc"]["forest"]),
        grass = D(config["inputs"]["lulc"]["grassland"]),
        othernat = D(config["inputs"]["lulc"]["othernat"]),
        urban = D(config["inputs"]["lulc"]["urban"]),
        water = D(config["inputs"]["lulc"]["water"]),
        pa = D(config["inputs"]["pa"]["wdpa"]),
        kba = D(config["inputs"]["pa"]["kba"])
    output:
        ebp = D(config["estimap"]["ebp"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_ebp.py"

rule calculate_human_input:
    input:
        urban = D(config["inputs"]["lulc"]["urban"]),
        road = D(config["inputs"]["road"]),
    output:
        human_input = D(config["estimap"]["human_inputs"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_human_input.py"

rule calculate_ros:
    input:
        ebp = D(config["estimap"]["ebp"]),
        human_input = D(config["estimap"]["human_inputs"])
    output:
        ros = D(config["estimap"]["ros"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_ros.py"

rule calculate_high_quality:
    input:
        ros = D(config["estimap"]["ros"])
    output:
        hq = D(config["estimap"]["high_quality"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_high_quality.py"

rule calculate_distance_to_hq:
    input:
        hq = D(config["estimap"]["high_quality"])
    output:
        distance_to_hq = D(config["estimap"]["distance_to_hq"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_distance_to_hq.py"

rule calculate_daily_rec:
    input:
        distance_to_hq = D(config["estimap"]["distance_to_hq"]),
        population = D(config["inputs"]["pop"])
    output:
        visits = D(config["daily_rec"]["visits"]),
        value = D(config["daily_rec"]["value"]),
        value_per_capita = D(config["daily_rec"]["value_per_capita"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_zonal_tcm.py"

rule calculate_nbt:
    input:
        distance_to_hq = D(config["estimap"]["distance_to_hq"]),
        admin_borders = D(config["inputs"]["admin_borders"]),
        overnight_stays = D(config["inputs"]["tourism"]["overnight_stays"]),
        accommodation_locations = D(config["inputs"]["tourism"]["accommodation_locations"])
    output:
        admin_gdf = D(config["nbt"]["nights"]),
        # value = D(config["nbt"]["value"]),
        # nights_per_capita = D(config["nbt"]["nights_per_capita"]),
        accommodation_zones = D(config["nbt"]["accommodation_zones"])  # Diagnostic output
    params:
        scores = config["scores"],
        nodata = config["nodata"],
        target_year = config.get("target_year", 2019),  # Default 2019
        tourism_type = config.get("tourism_type", "total")  # "international", "domestic", or "total"
    conda:
        config['conda']
    script:
        "scripts/calculate_nbt.py"