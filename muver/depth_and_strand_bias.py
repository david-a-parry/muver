import sys
import math
import numpy
import re
import gzip
import random
import pysam
from itertools import islice
from Bio import bgzf
from collections import defaultdict
from scipy.stats import norm
from scipy.optimize import curve_fit
from fitting import gaussian
from wrappers import samtools
from depth_distribution import calculate_depth_distribution
from depth_distribution import process_chromosome_values
import reference

def read_cnv_bedgraph(cnv_bedgraph, ploidy):
        '''
        Read in a bedGraph file to find ploidy at individual positions,
        used later to overwrite sample-wide ploidy.
        '''
        cnv_regions = dict()

        with open(cnv_bedgraph) as f:

            for line in f:

                chromosome, start, end, ploidy = line.strip().split()
                start = int(start) + 1
                end = int(end)
                ploidy = int(ploidy)

                for i in range(start, end + 1):
                    if ploidy != ploidy:
                        cnv_regions[(chromosome, i)] = ploidy

        return cnv_regions


def depth_and_strand_ratios_from_bam(input_bam, ref_fn, strand_output,
                                     strand_stats, depth_output, depth_stats,
                                     filtered_regions_output, ploidy=2,
                                     cnv_bedgraph=None, chrom_whitelist=set(),
                                     sample_size=1e8, p_threshold=0.0001,
                                     merge_window=1000):
    cnv_regions = dict()
    if cnv_bedgraph:
        cnv_regions = read_bedgraph(cnv_bedgraph, ploidy)
    chrom_sizes = reference.read_chrom_sizes(ref_fn)
    piter = samtools.mpileup_iter(input_bam, ref_fn)
    with gzip.open(strand_output, 'wt') as sout , \
            bgzf.BgzfWriter(depth_output) as dout:
        for line in piter:
            present_alleles = set()

            plus_tally = defaultdict(int)
            minus_tally = defaultdict(int)

            line_split = line.strip().split()

            chromosome, position, ref, coverage = line_split[:4]
            if chrom_whitelist and chromosome not in chrom_whitelist:
                continue
            position = int(position)
            coverage = float(coverage)
            if int(coverage) > 0:
                bases = line_split[4]
            else:
                bases = ''

            i = 0
            while i < len(bases):

                if bases[i] == '.':
                    present_alleles.add(ref)
                    plus_tally[ref] += 1
                elif bases[i] == ',':
                    present_alleles.add(ref)
                    minus_tally[ref] += 1
                elif re.match('[ACGT]', bases[i]):
                    present_alleles.add(bases[i])
                    plus_tally[bases[i]] += 1
                elif re.match('[acgt]', bases[i]):
                    present_alleles.add(bases[i].upper())
                    minus_tally[bases[i]] += 1

                elif re.match('[+-]', bases[i]):
                    indel_type = bases[i]
                    i += 1

                    indel_length = int(bases[i])
                    i += 1

                    indel = indel_type + bases[i:i+indel_length].upper()
                    present_alleles.add(indel)
                    if re.match('[ACGT]', bases[i:i+indel_length]):
                        plus_tally[indel] += 1
                    elif re.match('[acgt]', bases[i:i+indel_length]):
                        minus_tally[indel] += 1
                    i += indel_length - 1

                elif bases[i] == '^':
                    i += 1
                elif bases[i] == '*':
                    coverage += 1
                i += 1

            for allele in present_alleles:
                if allele in plus_tally and allele in minus_tally:
                    ratio = float(plus_tally[allele]) / minus_tally[allele]
                    sout.write("{}\t{}\t{}\t{}\t{}".format(
                        line_split[0], line_split[1], line_split[2], allele,
                        math.log(ratio)) + "\n")

            if coverage > 0:
                if (chromosome, position) in cnv_regions:
                    dp = int(coverage / cnv_regions[(chromosome, position)])
                else:
                    dp = int(coverage / ploidy)
                dout.write("{}\t{}\t{}\n".format(chromosome, position, dp))
    pysam.tabix_index(depth_output, seq_col=0, start_col=1, end_col=1,
                      force=True)
    sample_size = int(sample_size)
    sys.stderr.write("Subsampling from {}\n".format(strand_output))
    strand_sample = subsample_from_file(fh=gzip.open(strand_output),
                                        sample_size=sample_size,
                                        dtype=float,
                                        line_lambda=lambda x: x.split("\t")[4])
    if len(strand_sample) < sample_size:
        sys.stderr.write("WARNING: retrieved {:,}".format(
            len(strand_sample)) + " records from {}".format(strand_output)+
            " which is fewer than the requested subsample " +
            "({:,})".format(sample_size) + "\n")
    sys.stderr.write("Calculating strand bias stats\n")
    strand_mu, strand_sigma = calculate_strand_bias(strand_sample,
                                                    strand_stats)
    del strand_sample
    sys.stderr.write("Subsampling from {}\n".format(depth_output))
    dp_sample = subsample_from_file(fh=gzip.open(depth_output),
                                    sample_size=sample_size, dtype=int,
                                    line_lambda=lambda x: x.split("\t")[2])
    if len(dp_sample) < sample_size:
        sys.stderr.write("WARNING: retrieved {:,}".format(
            len(dp_sample)) + " records from {}".format(depth_output)+
            "which is fewer than the requested subsample " +
            "({:,})".format(sample_size) + "\n")
    sys.stderr.write("Calculating depth bias stats\n")
    dp_mu, dp_sigma = calculate_depth_distribution(dp_sample, depth_stats)
    del dp_sample
    sys.stderr.write("Filtering regions on depth\n")
    with open(filtered_regions_output, 'wt') as region_outfile:
        filter_regions_by_depth(depth_output, chrom_sizes, dp_mu, dp_sigma,
                                region_outfile, p_threshold, merge_window)

def parse_chrom_depth(tbx, chrom, chrom_sizes, mu, sigma, out,
                      p_threshold=0.0001, merge_window=1000, window=51):
    #collect depth at each nucleotide of a chromosome from tabix indexed depths
    depths = list(numpy.zeros(chrom_sizes[chrom], dtype=numpy.int32))
    for row in tbx.fetch(chrom):
        try:
            pos, dp = row.split("\t")[1:]
            depths[int(pos) - 1] = int(dp)
        except TypeError:
            sys.stderr.write("WARN: Error processing depth row: {}\n".format(
                row))
    depths = numpy.array(depths, dtype=numpy.int32)
    process_chromosome_values(chromosome=chrom, chromosome_values=depths,
                              mu=mu, sigma=sigma, OUT=out,
                              p_threshold=p_threshold,
                              merge_window=merge_window, window=window)


def filter_regions_by_depth(depth_output, chrom_sizes, mu, sigma,
                            out, p_threshold=0.0001, merge_window=1000):
    tbx = pysam.TabixFile(depth_output)
    for chrom in tbx.contigs:
        parse_chrom_depth(tbx, chrom, chrom_sizes, mu, sigma, out,
                          p_threshold=0.0001, merge_window=1000, window=51)


def calculate_strand_bias(log_ratios, output):
    p0_mu, p0_sigma = norm.fit(log_ratios)
    if p0_sigma == 0:
        p0_sigma = 0.01

    hist = numpy.histogram(log_ratios,
        bins=[float(x)/10 for x in range(-50, 51)], density=True)
    popt, pcov = curve_fit(gaussian, hist[1][:-1], hist[0],
        p0=[p0_mu, p0_sigma], maxfev=100000)
    mu, sigma = popt
    sigma = abs(sigma)
    with open(output, 'w') as OUT:

        OUT.write('Average log ratio: {}\n'.format(str(mu)))
        OUT.write(
            'Standard deviation of log ratios: {}\n\n'.format(str(sigma)))

        OUT.write('Bias distribution:\n\n')
        OUT.write('\t'.join(['Strand log ratio', 'Frequency', 'Fit value']) +
            '\n')

        for hist_value, _bin in zip(hist[0], hist[1]):
            OUT.write('\t'.join((
                str(round(_bin, 1)),
                str(hist_value),
                str(norm.pdf(_bin, mu, sigma)),
            )) + '\n')

    return mu, sigma


def subsample_from_file(fh, sample_size=1e8, dtype=float,
                        line_lambda=lambda x: x):
    sample_size = int(sample_size)
    reservoir = list(line_lambda(x) for x in islice(fh, sample_size))
    i = 0
    for i, line in enumerate(fh, sample_size):
        r = random.randint(0, i)
        if r < sample_size:
            reservoir[r] = line_lambda(line)
        if i % sample_size == 0:
            sys.stderr.write("Read {:,} lines for subsampling\n".format(i))
    sys.stderr.write("Finished reading {:,} lines for subsampling\n".format(i))
    return [dtype(x) for x in reservoir]

