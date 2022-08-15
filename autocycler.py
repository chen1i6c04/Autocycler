#!/usr/bin/env python

import os
import re
import sys
import shutil
import logging
import argparse
import subprocess
from tempfile import TemporaryDirectory


VERSION = {
    'kmc': 'kmc | grep KMC',
    'bwa': 'bwa 2>&1 | grep Version:',
    'rasusa': 'rasusa --version',
    'nanoq': 'nanoq --version',
    'unicycler': 'unicycler --version',
    'masurca': 'masurca --version',
    'polypolish_insert_filter.py': 'polypolish_insert_filter.py --version',
    'polypolish': 'polypolish --version',

}

FMT = "%(asctime)-20s[%(levelname)s] %(message)s"
DATEFMT = "%Y-%m-%d %H:%M:%S"


def run(cmd, catch_output=False):
    if catch_output:
        child_process = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
        return child_process
    else:
        subprocess.run(cmd, shell=True, check=True)


def check_dependency():
    logger = logging.getLogger(__name__)
    for program_name, cmd in VERSION.items():
        child_process = run(cmd, catch_output=True)
        if child_process.returncode:
            logger.error(msg=f"Could not determine version of {program_name}")
            sys.exit(0)
        else:
            # full_path = shutil.which(program_name)
            version = child_process.stdout.strip()
            logger.info(msg=f"Using {program_name:27} | {version}")


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = subprocess.getoutput(
            f"kmc -sm -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
        )
    return int(output.strip().split()[-1])


def parse_genome_size(pattern):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMG])', pattern)
    if result is None:
        sys.exit(f"Couldn't parse {pattern}")
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def run_polca(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    logger = logging.getLogger(__name__)
    logger.info("Running polca")
    
    logfile = logger.handlers[0].baseFilename
    assembly = os.path.abspath(assembly)
    os.chdir(output_dir)
    cmd = f"polca.sh -a {assembly} -r '{short_reads_1} {short_reads_2}' -t {num_threads} 2>&1 | " \
          f"tee /dev/tty | sed -r 's/\\x1b\\[[0-9;]*m//g' >> {logfile}"
    run(cmd)


def run_polypolish(assembly, alignments_1, alignments_2, output_dir):
    logger = logging.getLogger(__name__)
    logfile = logger.handlers[0].baseFilename

    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    polished_assembly = os.path.join(output_dir, 'polypolish.fasta')
    cmd = f"polypolish_insert_filter.py --in1 {alignments_1} --in2 {alignments_2} --out1 {filtered_1} --out2 {filtered_2} "\
          f"2>&1 | tee /dev/tty | sed -r 's/\\x1b\\[[0-9;]*m//g' >> {logfile}"
    run(cmd)
    cmd = f"polypolish {assembly} {filtered_1} {filtered_2} 2>&1 1> {polished_assembly} | " \
          f"tee /dev/tty | sed -r 's/\\x1b\\[[0-9;]*m//g' >> {logfile}"
    run(cmd)


def read_alignments(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    alignments_1 = os.path.join(output_dir, 'alignments_1.sam')
    alignments_2 = os.path.join(output_dir, 'alignments_2.sam')
    run(f"bwa-mem2 index {assembly}")
    run(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads_1} > {alignments_1}")
    run(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads_2} > {alignments_2}")
    return alignments_1, alignments_2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-1', '--short-reads-1', required=True)
    parser.add_argument('-2', '--short-reads-2', required=True)
    parser.add_argument('-l', '--long-reads', required=True)
    parser.add_argument('-o', '--output-dir', required=True)
    parser.add_argument('-t', '--num-threads', default=1, type=int)
    parser.add_argument('-g', '--genome-size', default=None)

    long_option = parser.add_argument_group("Long reads options")
    long_option.add_argument('-x', '--depth', default=100, type=int, help="default is 100")
    args = parser.parse_args()

    polca_dirname = os.path.join(args.output_dir, 'polish', 'polca')
    polypolish_dirname = os.path.join(args.output_dir, 'polish', 'polypolish')
    unicycler_assembly = os.path.join(args.output_dir, 'assembly.fasta')
    polca_assembly = os.path.join(polca_dirname, 'polypolish.fasta.PolcaCorrected.fa')
    polypolish_assembly = os.path.join(polypolish_dirname, 'polypolish.fasta')
    logfile = os.path.join(args.output_dir, 'autocycler.log')
    for dirname in (polca_dirname, polypolish_dirname):
        os.makedirs(dirname, exist_ok=True)

    # log setting
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    formatter = logging.Formatter(FMT, DATEFMT)
    
    file_handler = logging.FileHandler(logfile, mode="a", encoding=None, delay=False)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    check_dependency()

    logger.info(f"""Input data
    
    long reads          : {args.long_reads}
    forward short reads : {args.short_reads_1}
    reverse short reads : {args.short_reads_2}
""")

    if args.genome_size:
        genome_size = parse_genome_size(args.genome_size)
        logger.info(f"Using genome size was {genome_size}bp.")
    else:
        genome_size = estimate_genome_size(args.short_reads_1, args.num_threads)
        logger.info(f"Estimated genome size was {genome_size}bp.")

    output = subprocess.getoutput(f"zcat -f {args.short_reads_1} {args.short_reads_2} | nanoq -f -s")
    total_bases = int(output.split()[1])
    shot_depth = int(total_bases / genome_size)
    logger.info(f"Estimated short sequencing depth: {shot_depth}x.")
    
    output = subprocess.getoutput(f"nanoq -f -s -i {args.long_reads}")
    total_bases = int(output.split()[1])
    long_depth = int(total_bases / genome_size)
    logger.info(f"Estimated long sequencing depth: {long_depth}x.")

    if long_depth > args.depth:
        logger.info(f"Subsampling long reads from {long_depth}x to {args.depth}x")
        long_reads = os.path.join(args.output_dir, 'READS.sub.fq')
        run(f"rasusa --coverage {args.depth} --genome-size {genome_size} --input {args.long_reads} "
            f"--output {long_reads} 2>&1 | tee -a {logfile}")
    else:
        long_reads = args.long_reads

    logger.info("Running unicycler")
    run(f"unicycler -1 {args.short_reads_1} -2 {args.short_reads_2} -l {long_reads} "
        f"-o {args.output_dir} -t {args.num_threads} --spades_tmp_dir /tmp --no_correct --no_pilon")

    alignments_1, alignments_2 = read_alignments(
        unicycler_assembly, args.short_reads_1, args.short_reads_2, polypolish_dirname, args.num_threads
    )

    logger.info("Running polypolish")
    run_polypolish(
        unicycler_assembly,
        alignments_1,
        alignments_2,
        polypolish_dirname,
    )
    run_polca(
        polypolish_assembly,
        args.short_reads_1,
        args.short_reads_2,
        polca_dirname,
        args.num_threads,
    )
    final_assembly = os.path.join(args.output_dir, 'contigs.fasta')
    shutil.copyfile(polca_assembly, final_assembly)
    logger.info(f'Saving {final_assembly}')

    for filename in ('READS.sub.fq', 'READS.fit.fq', 'assembly.fasta.fai', 'polish'):
        file = os.path.join(args.output_dir, filename)
        if os.path.exists(file):
            if os.path.isdir(file):
                shutil.rmtree(file)
            else:
                os.remove(file)

    logger.info('Done')


if __name__ == '__main__':
    main()
