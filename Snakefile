import os
from snakemake import shell

# Define the input FASTA file
FASTA = "2024-01-05-reviewed-contam-UP000005640-UP000464024.fas"

# Define the working directory
WORKDIR = "reanalysis"
MZMLDIR = os.path.join(WORKDIR, "mzML")

# Define a wildcard to match all mzML files
MZML = os.listdir(MZMLDIR)

# Create a rule to compute
rule compute:
    input:
        expand(os.path.join(WORKDIR, "bullseye", "{iMZML}"), iMZML=MZML),
        expand(os.path.join(WORKDIR, "tide-index")),
        expand(os.path.join(WORKDIR, "tide-search", "{iMZML}"), iMZML=MZML),
        expand(os.path.join(WORKDIR, "concatenated_search_results", "combined.tsv")),
        expand(os.path.join(WORKDIR, "concatenated_search_results", "combined_fixed.tsv")),
        expand(os.path.join(WORKDIR, "percolator")),

        #expand(os.path.join(WORKDIR, "percolator", "{iMZML}"), iMZML=MZML)


# Create a rule to run Crux bullseye
rule crux_bullseye:
    input:
        mzML=os.path.join(WORKDIR, "mzML", "{iMZML}"),
    output:
        output_dir = os.path.join(WORKDIR, "bullseye", "{iMZML}")
    shell:
        "crux bullseye {input.mzML} {input.mzML} --output-dir {output.output_dir}"

# Create a rule to run Crux tide-index
rule crux_tide_index:
    input:
        fasta="2024-01-05-reviewed-contam-UP000005640-UP000464024.fas"
    params:
        index_name="tide-index"
    output:
        output_dir = os.path.join(WORKDIR, "tide-index")
    shell:
        "crux tide-index {input.fasta} {params.index_name} --output-dir {output.output_dir}"

# Create a rule to run Crux tide-index
rule crux_tide_search:
    input:
        mzML=os.path.join(WORKDIR, "mzML", "{iMZML}"),
        fasta="2024-01-05-reviewed-contam-UP000005640-UP000464024.fas"
    output:
        output_dir = os.path.join(WORKDIR, "tide-search", "{iMZML}")
    shell:
        "crux tide-search {input.mzML} {input.fasta} --output-dir {output.output_dir} --concat T"


rule concatenate_search_result:
    input:
        search_results=os.path.join(WORKDIR, "tide-search"),
    output:
        output = os.path.join(WORKDIR, "concatenated_search_results", "combined.tsv")
    shell:
        "./concatenate_tide_search_io.sh {input.search_results} {output.output}"

rule fix_target_decoy_column:
    input:
        search_results=os.path.join(WORKDIR, "concatenated_search_results", "combined.tsv"),
    output:
        output=os.path.join(WORKDIR, "concatenated_search_results", "combined_fixed.tsv"),
    shell:
        "python fix_target_decoy.py {input.search_results} {output.output}"

# Create a rule to run Crux tide-index
rule crux_percolator:
    input:
        search_results=os.path.join("tmp", "concat.txt"), #manually curated NEG1.txt + NEG2.txt
        #search_results=os.path.join(WORKDIR, "concatenated_search_results", "combined_fixed.tsv"), #maybe txt formatting matters....
        #search_results=os.path.join(WORKDIR, "tide-search", "NEG1.mzML","tide-search.txt"),
    output:
        output_dir = os.path.join(WORKDIR, "percolator")
    shell:
        "crux percolator {input.search_results} --output-dir {output.output_dir}"

# Create a rule to create the output directory if it does not exist
#rule all:
#    params:
#        bullseye_dir = os.path.join(WORKDIR, "bullseye"),
#    shell:
#        "mkdir -p {params.bullseye_dir}"



# Create a workflow
#snakemake.workflow(
#    runners={"shell": "bash"},
#    threads=1,
#    use_conda=True,
#    conda_prefix=".",
#)