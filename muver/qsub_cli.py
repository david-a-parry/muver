# -*- coding: utf-8 -*-

"""Console script for muver."""

import click
import os

from qsub_pipeline import run_pipeline as _run_qsub_pipeline


@click.group()
def main(args=None):
    pass

@main.command()
@click.option('--processes', '-p', default=1, type=int,
              help='Number of processes to use.')
@click.option('--excluded_regions', default=None, type=click.Path(exists=True),
              help='Regions to exclude in mutation calling (BED format).')
@click.option('--fwer', default=0.01, type=float,
              help='Familywise error rate. Default = 0.01.')
@click.option('--max_records', default=1000000, type=int,
              help='Maximum number of reads to store in memory when \
              sorting BAM files. Modify when any FASTQ pair \
              contains >1 billion reads. See manual for guidance \
              in setting this parameter. Default = 1000000.')
@click.option('--runtime', default=96, type=int,
              help='Runtime for qsub scripts in hours. Default=96.')
@click.option('--dummy_run', is_flag=True,
              help='Create qsub scripts but do not submit.')
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('fastq_list', type=click.Path(exists=True))
@click.argument('control_sample_name', type=str)
@click.argument('experiment_directory', type=str)
def run_qsub_pipeline(reference_assembly, fastq_list, control_sample_name,
                      experiment_directory, processes, excluded_regions,
                      fwer, max_records, qsub=False, runtime=96,
                      dummy_run=False):
    '''
    Run MuVer pipeline, starting with FASTQ files.

    FASTQ files and other parameters for each sample are specified in the
    FASTQ_LIST, a TXT file in tab-delimited format. For file specification,
    see tool documentation.

    Accepted fields for the FASTQ_LIST include:

    \b
    "Sample Name"
    "Mate 1 FASTQ"
    "Mate 2 FASTQ"
    "Ploidy"
    "CNV BedGraph"
    '''
    _run_qsub_pipeline(
        reference_assembly,
        fastq_list,
        control_sample_name,
        experiment_directory,
        p=processes,
        excluded_regions=excluded_regions,
        fwer=fwer,
        max_records=max_records,
        runtime=runtime,
        dummy_run=dummy_run,
    )

if __name__ == "__main__":
    main()
