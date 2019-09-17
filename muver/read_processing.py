import pysam

def remove_diff_chr_pairs(in_bam, out_bam):
    '''
    Go through a SAM file and remove entries where the other read is on a
    different chromosome.
    '''
    bamfile = pysam.AlignmentFile(in_bam, 'rb')
    good_pairs = pysam.AlignmentFile(out_bam, "wb", template=bamfile)
    for read in bamfile.fetch(until_eof=True):
        #if read.is_proper_pair: #TODO - check if this is author's original intention and intended handling of single-end reads
        if read.mrnm == read.rname:
            good_pairs.write(read)
    bamfile.close()
    good_pairs.close()
