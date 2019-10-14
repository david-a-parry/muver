# -*- coding: utf-8 -*-

"""Console script for muver."""

import click
import os
import sys

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
@click.option('--human_autosomes', is_flag=True,
               help='Only analyze data in human autosomes (chromosomes 1-22).')
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
        human_autosomes=human_autosomes,
        dummy_run=dummy_run,
    )

@main.command()
@click.option('--ploidy',  default=2, type=int,
              help='Sample ploidy.')
@click.option('--subsample_size', default=1e8, type=float,
              help='Amount of datapoints to subsample for strand bias and \
              depth analyses.')
@click.option('--human_autosomes', is_flag=True,
               help='Only analyze data in human autosomes (chromosomes 1-22).')
@click.option('--cnv_bedgraph_file', type=str, default=None,
              help='bedGraph file describing CNV regions')
@click.option('--p_threshold', type=float, default=0.0001,
              help='p-value threshold for abnormal depth')
@click.option('--merge_window', type=int, default=1000,
              help='maximum distance of adjacent abnormal sites \
              for creation of filtered regions')
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('reference_assembly', type=click.Path(exists=True))
@click.argument('output_prefix', type=str)

def depth_and_strand_bias_ratios(input_bam, reference_assembly, output_prefix,
                                 ploidy, subsample_size, human_autosomes,
                                 cnv_bedgraph_file, p_threshold, merge_window):
    ''' Calculate depth and strand biases from bam. Produce a list of regions
        to filter.'''
    from depth_and_strand_bias import depth_and_strand_ratios_from_bam
    chrom_whitelist = set()
    if human_autosomes:
        chrom_whitelist = [str(x) for x in range(1, 23)] #GRCh37 style
        chrom_whitelist += ["chr" + x for x in chrom_whitelist] #hg19/38 style
        chrom_whitelist = set(chrom_whitelist)
    strand_output = output_prefix + ".strand_bias.txt.gz"
    depth_output = output_prefix + ".depth.txt.gz"
    strand_stats = output_prefix + ".strand_bias_stats.txt"
    depth_stats = output_prefix + ".depth_stats.txt"
    filtered_regions_output = output_prefix + ".filtered_regions.bed"
    depth_and_strand_ratios_from_bam(input_bam, reference_assembly,
                                     strand_output, strand_stats, depth_output,
                                     depth_stats, filtered_regions_output,
                                     ploidy, cnv_bedgraph=cnv_bedgraph_file,
                                     chrom_whitelist=chrom_whitelist,
                                     sample_size=subsample_size,
                                     p_threshold=p_threshold,
                                     merge_window=merge_window)

@main.command()
@click.option('--human_autosomes', is_flag=True,
               help='Only analyze data in human autosomes (chromosomes 1-22).')
@click.argument('input_bam', type=click.Path(exists=True))
@click.argument('repeats_file', type=click.Path(exists=True))
@click.argument('output_prefix', type=str)

def repeat_rates(input_bam, repeats_file, output_prefix, human_autosomes):
    repeat_bgz = repeats_file + '.bgz'
    repeat_tbi = repeats_bgz + '.tbi'
    if not os.path.exists(repeat_file):
        sys.stderr.write('Repeats not found for reference assembly. Run '
            '"muver create_repeat_file".\n')
        exit()
    if not os.path.exists(repeat_bgz) or not os.path.exists(repeats_tbi):
        sys.stderr.write('{} or {} not found.'.format(repeat_bgz, repeat_tbi) +
                         ' Run muver_qsub compress_and_index_repeats".\n')
        exit()
    from repeat_indels import fit_repeat_indel_rates_low_mem
    chrom_whitelist = None
    if human_autosomes:
        chrom_whitelist = [str(x) for x in range(1, 23)] #GRCh37 style
        chrom_whitelist += ["chr" + x for x in chrom_whitelist] #hg19/38 style
        chrom_whitelist = set(chrom_whitelist)
    if os.path.isfile(repeats_file + '.sample'):
        repeats_file = repeats_file + '.sample'
    fit_repeat_indel_rates_low_mem(repeats_file, bam_file,
                                   output_prefix + '.repeat_indel_fits.txt',
                                   output_prefix + '.repeat_indel_header.txt',
                                   chromosome_whitelist=chrom_whitelist)

@main.command()
@click.argument('repeats_file', type=click.Path(exists=True))

def compress_and_index_repeats(repeats_file):
    try:
        import pysam
    except ImportError:
        raise RuntimeError("Can not import pysam. Please install pysam in " +
                           "order to write and index compressed repeats file.")
    output = repeats_file + '.bgz'
    sys.stderr.write("Compressing {} to {}\n".format(repeats_file, output))
    pysam.tabix_compress(repeats_file, output, force=True)
    sys.stderr.write("Indexing {}\n".format(output))
    pysam.tabix_index(output, force=True, seq_col=0, start_col=4, end_col=5)
    sys.stderr.write("Finished compressing and indexing {}".format(
                     repeats_file))

if __name__ == "__main__":
    main()
