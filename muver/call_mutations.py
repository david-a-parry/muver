import os
import sys

import reference
import sample
import variant_list


def call_mutations(reference_assembly, control_sample, sample_list, input_vcf,
                  output_header, chrom_sizes=None, excluded_regions=None):
    '''
    Using a reference assembly, sample list, and VCF file from GATK
    HaplotypeCaller, call mutations.

    chrom_sizes -- If specified, chromosome sizes are taken from here.
    excluded_regions -- Regions to exclude from variant calling (BED format).
    '''
    if not reference.check_reference_indices(reference_assembly):
        sys.stderr.write('Reference assembly not indexed. Run "muver '
            'index_reference".\n')
        exit()

    samples = sample.read_samples_from_text(sample_list)
    control_sample = next(
        (x for x in samples if x.sample_name == control_sample),
        None,
    )

    repeat_file = '{}.repeats'.format(
        os.path.splitext(reference_assembly)[0])
    if not os.path.exists(repeat_file):
        create_repeat_file(reference_assembly, repeat_file)

    if chrom_sizes:
        chrom_sizes = reference.read_chrom_sizes_from_file(chrom_sizes)
    else:
        chrom_sizes = reference.read_chrom_sizes(reference_assembly)

    variants = variant_list.VariantList(
        input_vcf, samples, excluded_regions, repeat_file,
        control_sample, chrom_sizes)

    text_output = '{}.mutations.txt'.format(output_header)
    vcf_output = '{}.mutations.vcf'.format(output_header)

    variants.write_output_table(text_output)
    variants.write_output_vcf(vcf_output)
