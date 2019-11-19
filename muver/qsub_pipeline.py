import os
import sys
from subprocess import call
import reference
from sample import (Sample, read_samples_from_text,
                    generate_experiment_directory, write_sample_info_file)
from math import ceil
from collections import defaultdict
import ConfigParser

PATHS = dict()
config = ConfigParser.ConfigParser()
config.read(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'paths.cfg'))
for key, value in config.items('paths'):
    PATHS[key] = value

_bash_check_log = '''for LOG in $@
do
    echo $(date) - checking previous exit status from $LOG
    if [ "$(tail -n 1 $LOG)" != '0' ]
    then
        echo non-zero status from $LOG - exiting
        exit 1
    fi
    echo $(date) - status OK
done
'''

class CommandGenerator(object):

    def __init__(self, sample, fq1, i, ref, fq2=None, processes=1,
                 align_mem=8, java_mem=8, java_overhead=4, runtime=96,
                 max_records=1000000):
        self.sample = sample
        self.fq1 = fq1
        self.fq2 = fq2
        self.ref = ref
        self.max_records = max_records
        self.runtime = runtime
        self.align_mem = float(align_mem)
        self.java_mem = float(java_mem)
        self.qsub_java_mem = float(java_mem + java_overhead)
        self.processes = processes
        self.i = i

    def _java_cmd(self, script_name, cmd, parallel=False):
        if parallel:
            pe = '#$ -pe sharedmem {}'.format(self.processes)
            mem = int(ceil(self.qsub_java_mem/self.processes))
        else:
            pe = ''
            mem = int(self.qsub_java_mem)
        return '''#!/bin/bash
#$ -e {}.stderr
#$ -o {}.stdout
#$ -V
#$ -cwd
#$ -l h_rt={}:00:00
#$ -l h_vmem={}G
{}
set -euo pipefail
{}

echo $(date) Starting
{}

echo $(date) Finished
echo $?
'''.format(script_name, script_name, self.runtime, mem, pe, _bash_check_log,
           cmd)


    def align(self, script_name):
        out = self.sample._mapq_filtered_sams[self.i].name
        self.sample._mapq_filtered_sams[self.i].close()
        if self.fq2:
            fq_string = "-1 {} -2 {} --fr".format(self.fq1, self.fq2)
            write_cmd = "perl -wane 'print if /^\@/ or $F[6] eq \"=\";' | "
        else:
            fq_string = "-U {}".format(self.fq1)
            write_cmd = ''
        cmd = '''
{} {} \\
    -q --phred33 -p {} \\
    -I 0 -X 1000 --local --sensitive-local \\
    -x {} | \\
    {}{} view \\
        -S -O BAM -q 20 \\
        -o {} -
'''.format(PATHS['bowtie2'], fq_string, self.processes, self.ref,
           write_cmd, PATHS['samtools'], out)
        q_script = '''#!/bin/bash
#$ -e {}.stderr
#$ -o {}.stdout
#$ -V
#$ -cwd
#$ -l h_rt={}:00:00
#$ -l h_vmem={}G
#$ -pe sharedmem {}

set -euo pipefail

echo $(date) Doing alignment
{}

echo $(date) Finished
echo $?
'''.format(script_name, script_name, self.runtime,
           int(ceil(self.align_mem/self.processes)), self.processes, cmd)
        return q_script

    def add_read_groups(self, script_name):
        in_bam = self.sample._mapq_filtered_sams[self.i].name
        out_bam = self.sample._read_group_bams[self.i].name
        self.sample._read_group_bams[self.i].close()
        cmd = '''java -Xmx{}g -jar {} \\
    AddOrReplaceReadGroups \\
    VALIDATION_STRINGENCY=SILENT \\
    SO=coordinate \\
    RGPL=illumina \\
    RGPU={} \\
    RGSM={} \\
    RGLB={} \\
    RGID={} \\
    I={} \\
    O={} \\
    TMP_DIR={} \\
    MAX_RECORDS_IN_RAM={}
'''.format(int(self.java_mem), PATHS['picard'], self.sample.sample_name,
           self.sample.sample_name, self.sample.sample_name,
           self.sample.sample_name, in_bam, out_bam,
           self.sample.tmp_dirs[self.i], self.max_records)
        return self._java_cmd(script_name, cmd)

    def deduplicate(self, script_name):
        in_bam = self.sample._read_group_bams[self.i].name
        out_bam = self.sample._deduplicated_bams[self.i].name
        metrics_file = self.sample._deduplication_metrics[self.i].name
        self.sample._deduplicated_bams[self.i].close()
        self.sample._deduplication_metrics[self.i].close()
        cmd = '''java -Xmx{}g -jar {} \\
    MarkDuplicates \\
    VALIDATION_STRINGENCY=SILENT \\
    REMOVE_DUPLICATES=TRUE \\
    I={} \\
    O={} \\
    M={} \\
    TMP_DIR={} \\
    CREATE_INDEX=true \\
    MAX_RECORDS_IN_RAM={}
'''.format(int(self.java_mem), PATHS['picard'], in_bam, out_bam, metrics_file,
           self.sample.tmp_dirs[self.i], self.max_records)
        return self._java_cmd(script_name, cmd)

    def realigner_target_creator(self, script_name):
        in_bam = self.sample._deduplicated_bams[self.i].name
        out_intervals = self.sample._interval_files[self.i].name
        self.sample._interval_files[self.i].close()
        cmd = '''java -Xmx{}g -jar {} \\
    -R {} \\
    -T RealignerTargetCreator \\
    -I {} \\
    -o {} \\
    -nt {}
'''.format(int(self.java_mem), PATHS['gatk'], self.ref, in_bam, out_intervals,
           self.processes)
        return self._java_cmd(script_name, cmd, parallel=True)

    def indel_realigner(self, script_name):
        in_bam = self.sample._deduplicated_bams[self.i].name
        intervals = self.sample._interval_files[self.i].name
        out_bam = self.sample._realigned_bams[self.i].name
        self.sample._realigned_bams[self.i].close()
        cmd = '''java -Xmx{}g -jar {} \\
    -R {} \\
    -T IndelRealigner \\
    --maxReadsForRealignment 100000 \\
    -I {} \\
    -targetIntervals {} \\
    -o {} \\
    -log {}
'''.format(int(self.java_mem), PATHS['gatk'], self.ref, in_bam, intervals, out_bam,
           self.sample.realignment_logs[self.i])
        return self._java_cmd(script_name, cmd, parallel=False)



    def fix_mate_information(self, script_name):
        in_bam = self.sample._realigned_bams[self.i].name
        out_bam = self.sample._fixed_mates_bams[self.i].name
        self.sample._fixed_mates_bams[self.i].close()
        if len(self.sample.fastqs) > 1: #need indices for merge command
            create_index = 'CREATE_INDEX=true'
        else:
            create_index = ''
        cmd = '''java -Xmx{}g -jar {} \\
    FixMateInformation \\
    VALIDATION_STRINGENCY=SILENT \\
    SO=coordinate \\
    I={} \\
    O={} \\
    TMP_DIR={} \\
    MAX_RECORDS_IN_RAM={} {}
'''.format(int(self.java_mem), PATHS['picard'], in_bam, out_bam,
           self.sample.tmp_dirs[self.i], self.max_records, create_index)
        return self._java_cmd(script_name, cmd)


def generate_qsub_scripts(fq1, fq2, ref, sample, i, p=1, max_records=1000000,
                         runtime=96):
    command_gen = CommandGenerator(sample=sample, ref=ref, fq1=fq1, fq2=fq2,
                                   i=i, processes=p, max_records=max_records,
                                   runtime=runtime)
    directory = os.path.join(sample.exp_dir, 'qsub_scripts', '')
    if not os.path.isdir(directory):
        os.makedirs(directory)
    qsub_stages = ['align', 'add_read_groups', 'deduplicate',
                   'realigner_target_creator', 'indel_realigner',
                   'fix_mate_information', ]
    scripts = list()
    for stage in qsub_stages:
        scr = os.path.join(directory, "{}_{}_{}.sh".format(sample.sample_name,
                                                           i, stage))
        with open(scr, 'wt') as outfile:
            outfile.write(getattr(command_gen, stage)(scr))
        scripts.append(scr)
    return scripts

def get_merge_and_index_scripts(sample, runtime=96, mem=4, processes=1):
    pe = ''
    if len(sample.fastqs) > 1:
        bam = sample.merged_bam.name
        sample.merged_bam.close()
        cmd = '{} merge -f '.format(PATHS['samtools'])
        if processes > 1:
            cmd += '-@ {} '.format(processes -1)
            pe = '#$ -pe sharedmem {}'.format(processes)
        cmd += '{}'.format( bam) + " ".join(sample._fixed_mates_bams)
    else:
        bam = sample._fixed_mates_bams[0].name
        cmd = ''
    directory = os.path.join(sample.exp_dir, 'qsub_scripts', '')
    scr = os.path.join(directory, "{}_merge.sh".format(sample.sample_name))
    with open(scr, 'wt') as outfile:
        outfile.write('''#!/bin/bash
#$ -e {}.stderr
#$ -o {}.stdout
#$ -V
#$ -cwd
#$ -l h_rt={}:00:00
#$ -l h_vmem={}G
{}
set -euo pipefail
{}

echo $(date) starting
{}
echo $(date) indexing
{} index {}
echo $(date) done
echo $?
'''.format(scr, scr, runtime, mem, pe, _bash_check_log, cmd, PATHS['samtools'],
           bam))
    return scr, bam

def submit_sample_scripts(sample2scripts, sample2merge, dummy=False):
    final_jobs = []
    for samp, script_lists in sample2scripts.items():
        for l in script_lists:
            prev_script = None
            for script in l:
                if prev_script:
                    cmd = ["qsub", "-hold_jid", os.path.basename(prev_script),
                           script, prev_script + '.stdout']
                else:
                    cmd = ["qsub", script]
                prev_script = script
                if dummy:
                    print(" ".join(cmd))
                else:
                    call(cmd)
        merge_script = sample2merge[samp]
        cmd = ["qsub", "-hold_jid", os.path.basename(prev_script),
               merge_script, prev_script + '.stdout']
        if dummy:
            print(" ".join(cmd))
        else:
            call(cmd)
        final_jobs.append(merge_script)
    return final_jobs

def submit_haplotype_caller(hc_script, bams, reference_assembly, vcf, log,
                            nct, holding_jobs, runtime=96, java_mem=8,
                            java_overhead=4, dummy=False):
    bam_args = "-I " + " \\\n    -I ".join(bams)
    hold_args = ",".join(os.path.basename(x) for x in holding_jobs)
    stdouts = " ".join(x + '.stdout' for x in holding_jobs)
    cmd = '''java -Xmx{}g -jar {} \\
    -T HaplotypeCaller \\
    -o {} \\
    -A StrandAlleleCountsBySample \\
    -A DepthPerSampleHC \\
    -R {} \\
    -nct {} \\
    -mmq 5 \\
    -log {} \\
    --minPruning 0 \\
    --minDanglingBranchLength 0 \\
    --pcr_indel_model NONE \\
    {}
'''.format(java_mem, PATHS['gatk'], vcf, reference_assembly, nct, log,
           bam_args)
    with open(hc_script, 'wt') as outfile:
        outfile.write('''#!/bin/bash
#$ -e {}.stderr
#$ -o {}.stdout
#$ -V
#$ -cwd
#$ -l h_rt={}:00:00
#$ -l h_vmem={}G
#$ -pe sharedmem {}

set -euo pipefail
{}

echo $(date) Calling variants
{}

echo $(date) Finished
echo $?
'''.format(hc_script, hc_script, runtime,
           int(ceil(float(java_mem + java_overhead)/nct)), nct,
           _bash_check_log, cmd))
    call_args = ['qsub', '-hold_jid', hold_args, hc_script, stdouts]
    if dummy:
        print(" ".join(call_args))
    else:
        call(call_args)

def submit_dp_strand_scripts(script_prefix, bams, reference_assembly,
                             holding_jobs, human_autosomes=False, runtime=96,
                             mem=16, dummy=False):
    h_flag = ""
    jobs = []
    if human_autosomes:
        h_flag = "--human_autosomes"
    hold_args = ",".join(os.path.basename(x) for x in holding_jobs)
    stdouts = " ".join(x + '.stdout' for x in holding_jobs)
    for bam in bams:
        prefix = os.path.splitext(bam)[0]
        cmd = '''muver_qsub depth_and_strand_bias_ratios {} \\
    {} \\
    {} \\
    {}'''.format(h_flag, bam, reference_assembly, prefix)
        script = script_prefix + "_" + os.path.basename(prefix) + ".sh"
        with open(script, 'wt') as outfile:
                outfile.write('''#!/bin/bash
#$ -e {}.stderr
#$ -o {}.stdout
#$ -V
#$ -cwd
#$ -l h_rt={}:00:00
#$ -l h_vmem={}G

set -euo pipefail
{}

echo $(date) Getting depth and strand bias stats from {}
{}

echo $(date) Finished
echo $?
'''.format(script, script, runtime, mem, _bash_check_log, bam, cmd))
        call_args = ['qsub', '-hold_jid', hold_args, script, stdouts]
        if dummy:
            print(" ".join(call_args))
        else:
            call(call_args)
        jobs.append(script)
    return jobs

def run_pipeline(reference_assembly, fastq_list, control_sample,
                 experiment_directory, p=1, excluded_regions=None,
                 fwer=0.01, max_records=1000000, runtime=96,
                 human_autosomes=False, dummy_run=False):
    '''
    Run the MuVer pipeline considering input FASTQ files.  All files written
    to the experiment directory.
    '''
    repeat_file = '{}.repeats'.format(os.path.splitext(reference_assembly)[0])
    repeat_bgz = repeat_file + '.bgz'
    repeat_tbi = repeat_bgz + '.tbi'
    if not reference.check_reference_indices(reference_assembly):
        sys.stderr.write('Reference assembly not indexed. Run '
            '"muver index_reference".\n')
        exit()
    if not os.path.exists(repeat_file):
        sys.stderr.write('Repeats not found for reference assembly. Run '
            '"muver create_repeat_file".\n')
        exit()
    if not os.path.exists(repeat_bgz) or not os.path.exists(repeats_tbi):
        sys.stderr.write('{} or {} not found.'.format(repeat_bgz, repeat_tbi) +
                         ' Run muver_qsub compress_and_index_repeats".\n')
        exit()
    generate_experiment_directory(experiment_directory)
    samples = read_samples_from_text(
        fastq_list, exp_dir=experiment_directory)
    control_sample = next(
        (x for x in samples if x.sample_name == control_sample),
        None,
    )

    for sample in samples:
        sample.generate_intermediate_files()
    # Align
    sample2bam = dict()
    sample2scripts = defaultdict(list)
    sample2merge = dict()
    for sample in samples:
        for i, fastqs in enumerate(sample.fastqs):
            if len(fastqs) == 2:
                f1, f2 = fastqs
            else:
                f1 = fastqs[0]
                f2 = None
            scripts = generate_qsub_scripts(fq1=f1, fq2=f2,
                                            ref=reference_assembly,
                                            sample=sample, i=i, p=p,
                                            max_records=max_records,
                                            runtime=runtime)
            sample2scripts[sample.sample_name].append(scripts)
        merge_script, final_bam = get_merge_and_index_scripts(sample)
        sample2merge[sample.sample_name] = merge_script
        sample2bam[sample.sample_name] = final_bam

    holding_jobs = submit_sample_scripts(sample2scripts, sample2merge,
                                         dummy=dummy_run)

    haplotype_caller_log = os.path.join(
        experiment_directory,
        'logs',
        'haplotype_caller.log'
    )
    haplotype_caller_vcf = os.path.join(
        experiment_directory,
        'gatk_output',
        'haplotype_caller_output.vcf.gz'
    )
    bams = list(sample2bam.values())
    hc_script = os.path.join(sample.exp_dir, 'qsub_scripts',
                             'haplotype_caller.sh')
    submit_haplotype_caller(
        hc_script,
        bams,
        reference_assembly,
        haplotype_caller_vcf,
        haplotype_caller_log,
        nct=p,
        holding_jobs=holding_jobs,
        runtime=runtime,
        dummy=dummy_run,
    )
    dp_strand_script = os.path.join(sample.exp_dir, 'qsub_scripts',
                                   'strand_and_depth_biases')
    bias_jobs = submit_dp_strand_scripts(dp_strand_script, bams,
                                         reference_assembly, holding_jobs,
                                         human_autosomes)

#    chrom_sizes = reference.read_chrom_sizes(reference_assembly)
#
#    strand_bias_std_values = pool.map(analyze_depth_distribution, zip(
#        range(len(samples)),
#        [s.get_intermediate_file_names() for s in samples],
#        repeat(reference_assembly),
#        repeat(chrom_sizes),
#        [s.ploidy for s in samples],
#        [s.cnv_regions for s in samples],
#    ))
#    for index, strand_bias_std in strand_bias_std_values:
#        samples[index].strand_bias_std = strand_bias_std
#
#    # Characterize repeats
#
#    if os.path.isfile(repeat_file + '.sample'):
#        repeats = read_repeats(repeat_file + '.sample')
#    else:
#        repeats = read_repeats(repeat_file)
#
#    pool.map(characterize_repeat_indel_rates, zip(
#        [s.get_intermediate_file_names() for s in samples],
#        repeat(repeats),
#        [s.repeat_indel_header for s in samples],
#    ))
#    for sample in samples:
#        sample.repeat_indel_fits_dict = read_fits(sample.repeat_indel_fits)
#
#    variants = VariantList(
#        haplotype_caller_vcf, samples, excluded_regions, repeat_file,
#        control_sample, chrom_sizes, fwer)
#
#    text_output = os.path.join(
#        experiment_directory,
#        'output',
#        'mutations.txt'
#    )
#    vcf_output = os.path.join(
#        experiment_directory,
#        'output',
#        'mutations.vcf'
#    )
#    variants.write_output_table(text_output)
#    variants.write_output_vcf(vcf_output)
#
#    for sample in samples:
#        sample.clear_temp_file_indices()
#
#    sample_info_file = os.path.join(
#        experiment_directory, 'sample_info.txt')
#    write_sample_info_file(samples, sample_info_file)
