# coding=utf-8

from pathlib import Path
import gzip
import shlex
from typing import List


def _open_text_auto(path: str, mode: str):
    """
    Open a text file transparently for .gz or plain files.

    Args:
        path: File path; if it ends with '.gz' a gzip handle is used.
        mode: 'rt' for read-text, 'wt' for write-text.
    """
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def split_vcf(input_vcf: str, output_dir: str, file_suffix: str, chunks: int = 10) -> None:
    """
    Split a VCF file into multiple chunks with an even distribution of variants.

    Parameters:
        input_vcf (str): Path to the input VCF (.vcf or .vcf.gz).
        output_dir (str): Directory to write chunked files.
        file_suffix (str): Prefix for output chunk file names.
        chunks (int): Number of chunks to create (default: 10). Must be >= 1.

    Returns:
        None. Writes files to `output_dir`.

    Output Files:
        {output_dir}/{file_suffix}.chunk_{0..N-1}.vcf  or  .vcf.gz (matches input extension)
    """
    if chunks < 1:
        raise ValueError("`chunks` must be >= 1")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    header: List[str] = []
    variants: List[str] = []

    # Read header and variants
    with _open_text_auto(input_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                header.append(line)
            else:
                variants.append(line)

    total_variants = len(variants)
    if total_variants == 0:
        # Still emit empty chunks with header for reproducibility
        for i in range(chunks):
            ext = ".vcf.gz" if input_vcf.endswith(".gz") else ".vcf"
            output_file = out_dir / f"{file_suffix}.chunk_{i}{ext}"
            with _open_text_auto(str(output_file), "wt") as w:
                w.writelines(header)
        return

    base = total_variants // chunks
    remainder = total_variants % chunks

    start = 0
    for i in range(chunks):
        size = base + (1 if i < remainder else 0)
        end = start + size
        chunk_variants = variants[start:end]

        ext = ".vcf.gz" if input_vcf.endswith(".gz") else ".vcf"
        output_file = out_dir / f"{file_suffix}.chunk_{i}{ext}"

        with _open_text_auto(str(output_file), "wt") as w:
            w.writelines(header)
            w.writelines(chunk_variants)

        start = end


def annotation_vcf(run_sample_id, sample, tool, paths, configure):
    """
    Run VEP annotation inside an Apptainer/Singularity container, post-filter,
    and then split the annotated VCF into chunks.

    Required keys (examples, do not hardcode absolute paths in code):
        paths['path']['neoantigen']['VEP_SIF']                  -> path to .sif image
        paths['database']['neoantigen']['VEP_DATA_DIR']         -> VEP cache directory
        paths['database']['neoantigen']['HG38']['REF_FASTA']    -> reference fasta
        paths['database']['neoantigen']['VEP_PLUGINS_VERSION']  -> plugins dir name (optional)

        configure['path']['output_dir']                         -> base output dir
        configure['args']['thread']                             -> int threads
        configure['step_name']['annotation']                    -> step name string
        configure['step_name']['vqsr']                          -> step name string

    The function expects `tool` to provide:
        - tool.judge_then_exec(run_id, cmd, expect_path)
        - tool.write_log(msg, level)
    """
    # Resolve parameters
    human_vep_sif = paths['path']['neoantigen']['VEP_SIF']
    vep_data_dir = paths['database']['neoantigen']['VEP_DATA_DIR']
    ref_fasta_path = Path(paths["database"]["neoantigen"]["HG38"]["REF_FASTA"]).expanduser()
    hg38_ref_dir = str(ref_fasta_path.parent)
    ref_name = ref_fasta_path.name

    # Optional: plugin version with fallback
    plugins_version = (
        paths.get("database", {})
             .get("neoantigen", {})
             .get("VEP_PLUGINS_VERSION", "VEP_plugins-release-110")
    )

    output_dir = configure['path']['output_dir'].rstrip("/")
    thread = int(configure['args']['thread'])
    hla_binding_threads = int(configure['args']['hla_binding_threads'])

    # Create output directory for this step
    vep_dir = f"{output_dir}/{sample}/{configure['step_name']['annotation']}/"
    cmd_mkdir = f"mkdir -p {shlex.quote(vep_dir)}"
    tool.judge_then_exec(run_sample_id, cmd_mkdir, vep_dir)

    # Input/previous step directory (e.g., VQSR output)
    data_dir = f"{output_dir}/{sample}/{configure['step_name']['vqsr']}/"
    suffix = "mutect.filtered"

    # Container-bind relative paths
    input_file = f"/data_dir/{sample}.{suffix}.vcf.gz"
    output_file = f"/output_dir/{sample}.{suffix}.VEP.vcf"
    abs_output_file = f"{vep_dir}{sample}.{suffix}.VEP.vcf"
    abs_filtered_vcf = f"{vep_dir}{sample}.{suffix}.VEP.filtered.vcf"

    # Notes:
    #   --pick: restrict to top transcript per variant to control downstream runtime.
    #   --transcript_version: append transcript version; useful if expression sources are versioned.
    vep_options = (
        "--offline --cache --format vcf --vcf --symbol --terms SO --tsl --biotype "
        "--hgvs --plugin Frameshift --plugin Wildtype --af --af_1kg --sift b "
        "--transcript_version --pick"
    )

    # Build Apptainer command
    cmd_args = [
        "apptainer", "exec",
        "-B", f"{vep_dir}:/output_dir/",
        "-B", f"{data_dir}:/data_dir/",
        "-B", f"{vep_data_dir}:/vep_data/",
        "-B", f"{hg38_ref_dir}:/fasta_data/",
        human_vep_sif,
        "vep",
        "--fork", str(thread),
        "--input_file", input_file,
        "--output_file", output_file,
        "--dir_cache", "/vep_data/",
        f"--fasta=/fasta_data/{ref_name}",
        f"--dir_plugins=/vep_data/{plugins_version}/",
    ] + vep_options.split()

    cmd1 = " ".join(shlex.quote(x) for x in cmd_args)

    # Run VEP
    tool.judge_then_exec(run_sample_id, cmd1, abs_output_file)

    # Post-process: report ref-transcript mismatches (requires `vatools`)
    cmd2 = f"ref-transcript-mismatch-reporter {shlex.quote(abs_output_file)} -f hard"
    tool.judge_then_exec(run_sample_id, cmd2, abs_filtered_vcf)

    # Keep only FILTER == PASS variants before chunking
    abs_filtered_pass_vcf = f"{vep_dir}{sample}.{suffix}.VEP.filtered.PASS.vcf"
    cmd_keep_pass = (
        "awk 'BEGIN{FS=OFS=\"\\t\"} /^#/ || $7==\"PASS\"' "
        f"{shlex.quote(abs_filtered_vcf)} > {shlex.quote(abs_filtered_pass_vcf)}"
    )
    tool.write_log("Filter VCF to FILTER==PASS for downstream analysis.", "info")
    tool.judge_then_exec(run_sample_id, cmd_keep_pass, abs_filtered_pass_vcf)

    # Split VCF into chunks equal to hla_binding_threads count (or choose another number as needed)
    chunk_dir = f"{vep_dir}chunks"
    cmd_mkdir_chunks = f"mkdir -p {shlex.quote(chunk_dir)}"
    tool.judge_then_exec(run_sample_id, cmd_mkdir_chunks, chunk_dir)

    tool.write_log(f"Split PASS-only VCF into {hla_binding_threads} chunks.", "info")
    split_vcf(
        input_vcf=abs_filtered_pass_vcf,
        output_dir=chunk_dir,
        file_suffix=f"{sample}.{suffix}.VEP.filtered.PASS",
        chunks=hla_binding_threads
    )
    tool.write_log("VCF splitting completed.", "info")
    return vep_dir
