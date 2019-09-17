# -*- coding: utf-8 -*-

"""Console script for muver."""

import click
import os

from allelic_fraction import get_allelic_fractions
from bias_distribution import calculate_bias_distribution_bam
from call_mutations import call_mutations as _call_mutations
from depth_correction import write_corrected_bedgraph
from depth_distribution import (calculate_depth_distribution_bedgraph,
                                filter_regions_by_depth_bedgraph)
from depth_ratios import calculate_depth_ratios as _calculate_depth_ratios
from pipeline import run_pipeline as _run_pipeline
from qsub_pipeline import run_pipeline as _run_qsub_pipeline
from reference import create_reference_indices, read_chrom_sizes
from repeat_indels import fit_repeat_indel_rates as _fit_repeat_indel_rates
from repeats import create_repeat_file as _create_repeat_file
from repeats import extract_repeat_file_sample as _extract_repeat_file_sample
from utils import read_repeats
from wrappers.samtools import get_mpileup_depths


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
@click.option('--qsub', is_flag=True,
              help='Create multistage commands via qsub job handler.')
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('fastq_list', type=click.Path(exists=True))
@click.argument('control_sample_name', type=str)
@click.argument('experiment_directory', type=str)

def run_pipeline(reference_assembly, fastq_list, control_sample_name,
                 experiment_directory, processes, excluded_regions,
                 fwer, max_records, qsub=False):
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
    if qsub:
        _run_qsub_pipeline(
            reference_assembly,
            fastq_list,
            control_sample_name,
            experiment_directory,
            p=processes,
            excluded_regions=excluded_regions,
            fwer=fwer,
            max_records=max_records,
        )
    else:
        _run_pipeline(
            reference_assembly,
            fastq_list,
            control_sample_name,
            experiment_directory,
            p=processes,
            excluded_regions=excluded_regions,
            fwer=fwer,
            max_records=max_records,
        )

@main.command()
@click.option('--excluded_regions', default=None, type=click.Path(exists=True),
              help='Regions to exclude in mutation calling (BED format).')
@click.option('--fwer', default=0.01, type=float,
              help='Familywise error rate.')
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('control_sample_name', type=str)
@click.argument('sample_list', type=click.Path(exists=True))
@click.argument('input_vcf', type=click.Path(exists=True))
@click.argument('output_header', type=str)
def call_mutations(reference_assembly, control_sample_name, sample_list,
                  input_vcf, output_header, excluded_regions, fwer):
    '''
    Call mutations from a HaplotypeCaller VCF file.

    Parameters for each sample are specified in the SAMPLE_LIST, a TXT file in
    tab-delimited format. For file specification, see tool documentation.

    Accepted fields for the SAMPLE_LIST include:

    \b
    "Sample Name"
    "Ploidy"
    "CNV BedGraph"
    "Strand Bias STD"
    "Filtered Sites"
    "Repeat Indel Fits"
    '''
    _call_mutations(
        reference_assembly,
        control_sample_name,
        sample_list,
        input_vcf,
        output_header,
        excluded_regions=excluded_regions,
        fwer=fwer,
    )

@main.command()
@click.argument('reference_assembly', type=click.Path(exists=True))
def index_reference(reference_assembly):
    '''
    Index the reference using Bowtie2, picard, and samtools.
    '''
    create_reference_indices(
        reference_assembly,
    )

@main.command()
@click.option('--mean', default=None, type=float,
              help='Manually specify the mean.')
@click.option('--ploidy', default=2, type=int,
              help='Sample ploidy.')
@click.argument('bedgraph_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_file', type=str)
def calculate_depth_ratios(bedgraph_file, reference_assembly, output_file,
                           mean, ploidy):
    '''
    Considering distance from chromosome ends, find the median ratio of depth
    to the average over 500-bp bins.
    '''
    _calculate_depth_ratios(
        bedgraph_file,
        reference_assembly,
        output_file,
        mean=mean,
        ploidy=ploidy,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_file', type=str)
def plot_allelic_fraction(bam_file, reference_assembly, output_file):
    '''
    Plot allelic fractions over a BAM file.

    The OUTPUT_FILE contains a histogram of allelic fractions (tab-delimited
    TXT format).
    '''
    get_allelic_fractions(
        bam_file,
        reference_assembly,
        output_file,
    )

@main.command()
@click.option('--output_repeat_file', type=str,
    help='Specify the output file; otherwise generates the file in place.')
@click.argument('fasta_file', type=click.Path(exists=True))
def create_repeat_file(fasta_file, output_repeat_file):
    '''
    Create repeat file for the FASTA_FILE sequence.
    '''
    if output_repeat_file:
        out = output_repeat_file
    else:
        out = '{}.repeats'.format(os.path.splitext(fasta_file)[0])
    _create_repeat_file(
        fasta_file,
        out,
    )

@main.command()
@click.argument('repeat_file', type=click.Path(exists=True))
@click.argument('sample_size', type=int)
def extract_repeat_file_sample(repeat_file, sample_size):
    '''
    Extract a random sample of repeats.
    '''
    sample_file = repeat_file + '.sample'
    _extract_repeat_file_sample(
        repeat_file,
        sample_file,
        sample_size,
    )

@main.command()
@click.argument('y_int', type=float)
@click.argument('scalar', type=float)
@click.argument('mean_log', type=float)
@click.argument('sd_log', type=float)
@click.argument('slope', type=float)
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('input_bedgraph', type=click.Path(exists=True))
@click.argument('output_bedgraph', type=str)
def correct_depths(y_int, scalar, mean_log, sd_log, slope, reference_assembly,
                     input_bedgraph, output_bedgraph):
    '''
    Correct values in a depth bedGraph file.

    Correction is performed using the sum of a log-normal cumulative
    distribution function and linear function.
    '''
    chrom_sizes = read_chrom_sizes(reference_assembly)
    write_corrected_bedgraph(
        input_bedgraph,
        chrom_sizes,
        output_bedgraph,
        y_int,
        scalar,
        mean_log,
        sd_log,
        slope,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_bedgraph', type=str)
def calculate_read_depths(bam_file, reference_assembly, output_bedgraph):
    '''
    Give depths in a BAM_FILE in bedGraph format.

    Depth values are calculated using the mpileup function of samtools.
    '''
    get_mpileup_depths(
        bam_file,
        reference_assembly,
        output_bedgraph,
    )

@main.command()
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_bias_distribution', type=str)
def calculate_bias_distribution(bam_file, reference_assembly,
                                output_bias_distribution):
    '''
    Calculate distribution of strand biases.
    '''
    calculate_bias_distribution_bam(
        bam_file,
        reference_assembly,
        output_bias_distribution,
    )

@main.command()
@click.option('--output_filtered_regions', type=str, default=None,
              help='If OUTPUT_FILTERED_REGIONS is specified, regions to be '
              'filtered based on abnormal depth will be written to this file.')
@click.option('--ploidy', type=int, default=2, help='Global ploidy')
@click.option('--cnv_bedgraph_file', type=str, default=None,
              help='bedGraph file describing CNV regions')
@click.option('--p_threshold', type=float, default=0.0001,
              help='p-value threshold for abnormal depth')
@click.option('--merge_window', type=int, default=1000,
              help='maximum distance of adjacent abnormal sites \
              for creation of filtered regions')
@click.argument('bedgraph_file', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_depth_distribution', type=str)
def calculate_depth_distribution(bedgraph_file, reference_assembly,
                                 output_depth_distribution,
                                 output_filtered_regions, ploidy,
                                 cnv_bedgraph_file, p_threshold,
                                 merge_window):
    '''
    Calculate distribution of depths in a bedGraph file.
    '''
    mu, sigma = calculate_depth_distribution_bedgraph(
        bedgraph_file,
        output_depth_distribution,
        ploidy,
        cnv_bedgraph_file,
    )
    chrom_sizes = read_chrom_sizes(reference_assembly)
    if output_filtered_regions:
        filter_regions_by_depth_bedgraph(
            bedgraph_file,
            chrom_sizes,
            mu,
            sigma,
            output_filtered_regions,
            ploidy,
            cnv_bedgraph_file,
            p_threshold,
            merge_window,
        )

@main.command()
@click.option('--output_plot_header', type=str, default=None,
              help='If OUTPUT_PLOT_HEADER is specified, PNGs are written '
              'that show the fit relative to the observed repeat indel rates.')
@click.argument('bam_file', type=click.Path(exists=True))
@click.argument('repeats_file', type=click.Path(exists=True))
@click.argument('output_fits_file', type=str)
def fit_repeat_indel_rates(bam_file, repeats_file, output_fits_file,
                           output_plot_header):
    '''
    Fit a logistic function to log-transformed repeat indel rates.

    The OUTPUT_FITS_FILE gives the parameters for the derived fits.
    '''
    repeats = read_repeats(repeats_file)
    _fit_repeat_indel_rates(
        repeats,
        bam_file,
        output_fits_file,
        output_plot_header,
    )

if __name__ == "__main__":
    main()
