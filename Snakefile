from pathlib import Path
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

# Load config
configfile: "config.yaml"


S3 = S3RemoteProvider(
    access_key_id="8CBN25LZUNYWLJZPQSXQ",
    secret_access_key="cbU7On5L3peM5AiPF7Oa7G55aW6FDIIY5Uh0ynyk",
    endpoint_url="https://s3.msi.umn.edu",
    # region="us-east-1",  # or whatever MSI pretends to use
    # config={"signature_version": "s3v4"}
)

# Get data_dir from config
data_dir = config["data_dir"]
USE_S3 = str(data_dir).startswith("s3://") or str(data_dir).startswith("s3:")
# if USE_S3:
#     S3 = S3RemoteProvider()
#     def D(path_str):
#         return S3[f"{data_dir.rstrip('/')}/{path_str}"]
# else:
#     DATA_DIR = Path(data_dir).expanduser()
#     def D(path_str):
#         return str(DATA_DIR / path_str)
if USE_S3:
    S3 = S3RemoteProvider()
    def D(path_str):
        # Remove s3:// prefix if present for remote()
        s3_path = data_dir
        if s3_path.startswith("s3://"):
            s3_path = s3_path[len("s3://"):]
        elif s3_path.startswith("s3:"):
            s3_path = s3_path[len("s3:"):]

        # Join base s3 path and relative path with slash
        full_path = f"{s3_path.rstrip('/')}/{path_str.lstrip('/')}"
        # Return remote path object
        return S3.remote(full_path)
        print(full_path)
else:
    DATA_DIR = Path(data_dir).expanduser()
    def D(path_str):
        return str(DATA_DIR / path_str)

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