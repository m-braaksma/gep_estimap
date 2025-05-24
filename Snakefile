# Load config
configfile: "config.yaml"

# Define tier 2 storage settings
# note: requires snakemake-storage-plugin-s3
storage:
    provider="s3",
    endpoint_url="https://s3.msi.umn.edu",

# Define flexible paths depending on data dir
data_dir = config["data_dir"]
def D(path_str):
    if str(config["data_dir"]).startswith("s3://"):
        return storage.s3(os.path.join(data_dir, path_str))
    else:
        return str(data_dir.joinpath(path_str))


rule all:
    input:
        D(config["estimap"]["ebp"]),
        D(config["estimap"]["human_inputs"]),
        D(config["estimap"]["ros"]),
        D(config["estimap"]["high_quality"]),
        D(config["estimap"]["distance_to_hq"]),
        D(config["estimap"]["visits"]),
        D(config["estimap"]["value"]),

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

rule calculate_zonal_tcm:
    input:
        distance_to_hq = D(config["estimap"]["distance_to_hq"]),
        population = D(config["inputs"]["pop"])
    output:
        visits = D(config["estimap"]["visits"]),
        value = D(config["estimap"]["value"])
    params:
        scores = config["scores"],
        nodata = config["nodata"]
    conda:
        config['conda']
    script:
        "scripts/calculate_zonal_tcm.py"