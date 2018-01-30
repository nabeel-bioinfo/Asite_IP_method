from __future__ import division
from optparse import OptionParser
from Bio import SeqIO
import sys
import subprocess
import os


def samparser_genome(sfile, frag_min, frag_max, three_prime=False):
    # Parse the SAM file and quantify mapped reads to chromosome positions.
    dict_count = {}  # This dictionary will contain the count values and other attributes for every genomic position
    dict_mul_count = {}  # This dictionary will contain the count values of multiple mapped reads
    dict_total = {}
    # Multiple mapped reads are counted separately and if they are greater than certain % of total reads mapped on a gene, that gene will not be considered for analysis
    log_file = open('sam_parser_logfile.log', 'w')
    frag_range = frag_max - frag_min + 1
    read_count = 0
    mul_read_count = 0
    if three_prime:
        print 'Mapping reads by 3\' end...'
    else:
        print 'Mapping reads by 5\' end...'
    counter = 0
    # Open the SAM file and parse the mapped reads according to chromosome positions, read length and strand
    with open(sfile) as sam_file:
        for line in sam_file:
            # Ignore header files
            if line[0] != '@':
                fields = line.split('\t')
                if len(fields) > 19:
                    multi_flag = fields[19].strip()  # Sometimes there are "NH:i:1\n"
                    # If the read has more than one alignment then report it as multiple mapping
                    if multi_flag != 'NH:i:1':
                        multi_align = True
                    else:
                        multi_align = False
                else:
                    multi_align = False
                sam_flag = int(fields[1])
                chr_no = fields[2]
                read_length = len(fields[9])

                if sam_flag in (0, 256):
                    # If mapping by 3' end
                    if three_prime:
                        position = int(fields[3]) + read_length - 1
                    # If mapping by 5' end
                    else:
                        position = int(fields[3])

                elif sam_flag in (16, 272):
                    # If mapping by 3' end
                    if three_prime:
                        position = int(fields[3])
                    # If mapping by 5' end
                    else:
                        position = int(fields[3]) + read_length - 1

                # Initialize the dict if not done already
                if chr_no not in dict_count:
                    dict_count[chr_no] = {}
                    dict_total[chr_no] = 0
                # Initialize an array of zeroes corresponding to reads of different lengths. Two values for each length representing the positive and negative strands
                # The last two zeroes will correspond to multiple aligned reads according to strand
                if position not in dict_count[chr_no]:
                    dict_count[chr_no][position] = [0] * (frag_range * 2 + 2)

                # If read is multiple mapped, then initialize the multiple count dict as well
                if multi_align and chr_no not in dict_mul_count:
                    dict_mul_count[chr_no] = {}
                if multi_align and position not in dict_mul_count[chr_no]:
                    dict_mul_count[chr_no][position] = [0] * (frag_range * 2 + 2)

                pos_in_value = (read_length - frag_min) * 2
                # Count the read according to its length and strand
                if frag_min <= read_length <= frag_max:
                    try:
                        if sam_flag == 0:  # Primary alignment on forward strand
                            if not multi_align:
                                dict_count[chr_no][position][pos_in_value] += 1
                                dict_total[chr_no] += 1
                                read_count += 1
                            else:
                                # Multiple mapped reads are counted separately
                                dict_count[chr_no][position][-2] += 1
                                dict_mul_count[chr_no][position][pos_in_value] += 1
                                mul_read_count += 1
                        if sam_flag == 16:  # Primary alignment on reverse strand
                            if not multi_align:
                                dict_count[chr_no][position][pos_in_value + 1] += 1
                                dict_total[chr_no] += 1
                                read_count += 1
                            else:
                                # Multiple mapped reads are counted separately. Last two columns are initialized for mul mapped reads for +ve and -ve strands respectively
                                dict_count[chr_no][position][-1] += 1
                                dict_mul_count[chr_no][position][pos_in_value + 1] += 1
                                mul_read_count += 1
                        # Not primary alignment. It will counted under multiple aligned reads
                        if sam_flag == 256:
                            position = int(fields[3])
                            dict_count[chr_no][position][-2] += 1
                            dict_mul_count[chr_no][position][pos_in_value] += 1

                        if sam_flag == 272:  # Not primary alignment and on reverse strand
                            position = int(fields[3]) + read_length - 1
                            dict_count[chr_no][position][-1] += 1
                            dict_mul_count[chr_no][position][pos_in_value + 1] += 1
                    except KeyError:
                        log_file.write("The KeyError line is \n " + line + "\n")
            counter += 1
            sys.stdout.write("Line of SAM file currently being parsed: {0}.\t\r".format(counter))
            sys.stdout.flush()

    print '\nSAM file parsed for total ' + str(read_count) + ' reads mapping onto ' + str(len(dict_count)) + ' chromosomes.'
    print str(mul_read_count) + ' reads are multiple aligned mapped to ' + str(len(dict_mul_count)) + ' chromosomes.'

    return dict_count, dict_mul_count


def samparser_transcriptome(sfile, frag_min, frag_max, three_prime=False):
    # Parse the SAM file and quantify mapped reads to chromosome positions.
    dict_count = {}  # This dictionary will contain the count values and other attributes for every genomic position
    dict_mul_count = {}  # This dictionary will contain the count values of multiple mapped reads
    total_count = {}
    # Multiple mapped reads are counted separately and if they are greater than 0.1% of total reads mapped on a gene, that gene will not be considered for analysis
    log_file = open('sam_parser_logfile.log', 'w')
    frag_range = frag_max - frag_min + 1
    read_count = 0
    mul_read_count = 0
    discarded_reads = 0
    if three_prime:
        print 'Mapping reads by 3\' end...'
    else:
        print 'Mapping reads by 5\' end...'
    counter = 0
    with open(sfile) as sam_file:
        for line in sam_file:
            # Ignore header files
            if line[0] != '@':
                fields = line.split('\t')
                sam_flag = int(fields[1])
                gene = fields[2]
                # if gene not in filtered_genes:
                #     continue
                if len(fields) > 19:
                    multi_flag = fields[19].strip()  # Sometimes there are "NH:i:1\n"
                    # If the read has more than one alignment then report it as multiple mapping
                    if multi_flag != 'NH:i:1':
                        multi_align = True
                    else:
                        multi_align = False
                else:
                    multi_align = False

                read_length = len(fields[9])
                pos_in_value = (read_length - frag_min)
                if gene not in dict_count:
                    dict_count[gene] = {}
                    total_count[gene] = 0
                # If mapping by 3' end
                if three_prime:
                    position = int(fields[3]) + read_length - 1
                # If mapping by 5' end
                else:
                    position = int(fields[3])
                if position not in dict_count[gene]:
                    dict_count[gene][position] = [0] * (frag_range + 1)
                if multi_align:
                    if gene not in dict_mul_count:
                        dict_mul_count[gene] = {}
                    if position not in dict_mul_count[gene]:
                        dict_mul_count[gene][position] = [0] * (frag_range + 1)

                if frag_min <= read_length <= frag_max:
                    try:
                        if sam_flag == 0:  # Primary alignment on forward strand
                            if not multi_align:
                                dict_count[gene][position][pos_in_value] += 1
                                read_count += 1
                                total_count[gene] += 1
                            else:
                                # Multiple mapped reads are counted separately
                                dict_count[gene][position][-1] += 1
                                dict_mul_count[gene][position][pos_in_value] += 1
                                mul_read_count += 1
                        # Not primary alignment. It will counted under multiple aligned reads
                        elif sam_flag == 256:
                            dict_count[gene][position][-1] += 1
                            dict_mul_count[gene][position][pos_in_value] += 1
                        elif sam_flag == 16:
                            discarded_reads += 1

                    except KeyError:
                        log_file.write("The KeyError line is \n " + line + "\n")
            counter += 1
            sys.stdout.write("Line of SAM file currently being parsed: {0}.\t\r".format(counter))
            sys.stdout.flush()

    print '\nSAM file parsed for total ' + str(read_count) + ' reads mapping onto ' + str(len(dict_count)) + ' genes.'
    print str(mul_read_count) + ' reads are multiple aligned mapped to ' + str(len(dict_mul_count)) + ' genes.'
    print str(discarded_reads) + ' reads were mapped spuriously.'
    return dict_count, dict_mul_count, total_count


def create_cds_table_transcriptome(idx_file, seq_file, output, sam_count_dict, dict_mul_count, total_count, frag_min, frag_max, three_prime=False, filter_file='', fast_mode=True):
    outfile = open(output+'CDS_count_table', 'w')
    logfile = open(output+'makecdstable.log', 'w')
    frag_range = frag_max - frag_min + 1
    mul_gene_list = []
    dict_mul_genes = {}
    dict_len = {}
    idx_dict = {}
    dict_count_len = {}
    dict_mul_count_len = {}
    # These dicts will contain dicts of genes and their read counts (as lists) according to fragment length
    for fsize in range(frag_min, frag_max + 1):
        dict_count_len[fsize] = {}

    with open(idx_file) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            # fields = [gene name, CDS start position, CDS length, Start codon, End codon]
            idx_dict[fields[0]] = [int(fields[1]), int(fields[2]), fields[3], fields[4]]

    seq_dict = {}

    for seq_record in SeqIO.parse(seq_file, "fasta"):
        gene_id = str(seq_record.id).split(' ')[0]
        seq_dict[gene_id] = seq_record.seq

    if filter_file:
        filter_genes = True
        select_genes = []
        with open('/gpfs/group/epo2/default/nxa176/reference/Mus_musculus/single_isoform/mm10_start_index_single_isoforms.tab') as f:
            for lines in f:
                fields = lines.strip().split('\t')
                gene = fields[0]
                select_genes.append(gene)
        canonical_genes = []
        with open('/gpfs/group/epo2/default/nxa176/reference/Mus_musculus/Canonical_orfs_tables3_ingolia_mm10_included_missing_identifiers.txt') as f:
            for lines in f:
                canonical_genes.append(lines.strip())
    else:
        filter_genes = False

    dict_start = {}
    # gene_count_file = open('Gene_counts_and_length.tab', 'w')
    print 'Starting to make the CDS table'
    for gene in sam_count_dict:
        try:
            start_pos, cds_len = idx_dict[gene][:2]
        except KeyError:
            logfile.write('No index available for gene ' + gene + '\n')
            continue
        # gene_count_file.write(gene + '\t' + str(total_count[gene]) + '\t' + str(cds_len) + '\n')

        if start_pos in [0, -1]:
            logfile.write('No UTR region for gene ' + gene + '\n')
            continue
        # If a gene has very sparse reads, it is better to leave it out as it will not meet the filtering criteria.
        # below criteria states that avg is less than 1. This is done for faster processing.
        if fast_mode and total_count[gene] < cds_len:
            logfile.write('SPARSELY POPULATED - Being filtered:\t' + gene + '\t' + str(total_count[gene]) + '\t' + str(cds_len) + '\n')
            continue
        dict_len[gene] = cds_len
        multi_genes = 'N'
        start_utr = -start_pos
        cds_pos = start_utr
        for pos in range(1, len(seq_dict[gene]) + 1):
            # for pos in sorted(sam_count_dict[gene]):
            nuc = seq_dict[gene][pos - 1]
            # We are not interested in positions beyond 50 nt around the CDS region. Hence we skip these positions
            if cds_pos < -50 or cds_pos > cds_len + 50:
                cds_pos += 1
                continue
            if cds_pos == 0:
                cds_pos += 1
            try:
                count_list = sam_count_dict[gene][pos][0:frag_range]
            except KeyError:
                count_list = [0] * frag_range

            try:
                multi_aligned = sam_count_dict[gene][pos][-1]
            except KeyError:
                # print 'KeyError in multi aligned for gene ' + gene + ' at position ' + str(pos)
                multi_aligned = 0

            try:
                if multi_aligned > 0:
                    mul_count_list = dict_mul_count[gene][pos][0:frag_range]
            except KeyError:
                mul_count_list = [0] * frag_range

            outfile.write(gene + '\t' + str(cds_pos) + '\t' + '\t'.join(map(str, ['X', cds_len, nuc, 'transcriptome', pos, multi_aligned, multi_genes, sum(count_list)] +
                                                                            count_list)) + '\n')
            # We will create dictionaries with fragment length as keys and dict of gene with list of read counts at every position as our values
            for fsize in range(frag_min, frag_max + 1):
                if gene not in dict_count_len[fsize]:
                    dict_count_len[fsize][gene] = []
                    dict_start[gene] = cds_pos
                dict_count_len[fsize][gene].append(count_list[fsize - frag_min])

            if multi_aligned > 0 and gene not in mul_gene_list:
                mul_gene_list.append(gene)
                if gene not in dict_mul_genes:
                    dict_mul_genes[gene] = {}
            if multi_aligned > 0:
                dict_mul_genes[gene][cds_pos] = ['X', cds_len, nuc, 'transcriptome', pos, multi_aligned, multi_genes] + mul_count_list
                # As done for actual counts, multiple mapped reads are counted according to fragment length
                for fsize in range(frag_min, frag_max + 1):
                    if gene not in dict_mul_count_len:
                        dict_mul_count_len[gene] = {}
                    if fsize not in dict_mul_count_len[gene]:
                        dict_mul_count_len[gene][fsize] = {0: 0, 1: 0, 2: 0}
                    if cds_pos < 0:
                        frame = cds_pos % 3
                    else:
                        frame = cds_pos % 3 - 1
                        if frame == -1:
                            frame = 2
                    if three_prime:
                        pos_range = range(1, cds_len + frag_max)
                    else:
                        pos_range = range(-frag_max, cds_len)
                    if cds_pos in pos_range:
                        dict_mul_count_len[gene][fsize][frame] += mul_count_list[fsize - frag_min]

            cds_pos += 1

    # Write out the number of mutliple mapped reads for each gene according to fragment size and frame. For each fragment size, write the reads for frame 0, 1 and 2.
    mul_out_file = open(output+'Multiple_mapped_gene_read_counts_' + str(frag_min) + '_' + str(frag_max) + '.tab', 'w')

    for gene in dict_mul_count_len:
        mul_out_file.write(str(gene))
        for fsize in range(frag_min, frag_max + 1):
            for frame in xrange(3):
                mul_out_file.write('\t' + str(dict_mul_count_len[gene][fsize][frame]))
        mul_out_file.write('\n')
    mul_out_file.close()

    for fsize in range(frag_min, frag_max + 1):
        count_file = open(output+'Read_counts_' + str(fsize) + '.tab', 'w')
        for gene, reads_list in dict_count_len[fsize].iteritems():
            start_pos, cds_len = idx_dict[gene][:2]
            start_idx = dict_start[gene]
            if filter_genes:
                if gene in select_genes and gene in canonical_genes:
                    count_file.write(gene + '\t' + str(start_idx) + '\t' + str(cds_len) + '\t' + ','.join(map(str, reads_list)) + '\n')
            else:
                count_file.write(gene + '\t' + str(start_idx) + '\t' + str(cds_len) + '\t' + ','.join(map(str, reads_list)) + '\n')
        count_file.close()


def make_cds_table(annotation_file, genome, output, sam_parsed_count_dict, dict_mul_count, frag_min, frag_max, remove_overlapped_genes=True, three_prime=False, extra_overlap=0):
    # This module will take the count dictionary parsed from the sam file containing the read counts at each genomic position and map it to gene positions
    outfile = open(output+'CDS_count_table', 'w')

    complimentarydict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    dict_gene, dict_cds_count, dict_cds_info, genome_dict, overlap_genes, dict_len = cdsparser(annotation_file, genome, extra_overlap=extra_overlap)

    dict_count_len = {}
    dict_mul_count_len = {}
    # These dicts will contain dicts of genes and their read counts (as lists) according to fragment length
    for fsize in range(frag_min, frag_max + 1):
        dict_count_len[fsize] = {}

    counter = 0

    # For every gene in the dataset
    for gene_name in dict_cds_count:
        chr_num, strand = dict_gene[gene_name]
        print gene_name, chr_num, strand
        # Do not include genes which have overlapping annotations as the source of the reads will be ambiguous
        if gene_name in overlap_genes:
            multi_genes = 'Y'
        else:
            multi_genes = 'N'
        frag_range = frag_max - frag_min + 1
        try:
            gene_length = dict_len[gene_name]
        except KeyError:
            print dict_cds_info[gene_name]
            print dict_len[gene_name]
            print 'KeyError in multicdslen for dict_cds_info in gene ' + str(gene_name)

        # dict_cds_info contains lists of CDS start and end positions as a list. For example if there are two exons for a gene X
        # dict_cds_info['X'] = [[111234, 111345],[112122,112543]]
        # We get reads at 50 nt around the first cds position and the last cds position and update these lists
        dict_cds_info[gene_name][0][0] -= 50
        dict_cds_info[gene_name][-1][1] += 50
        # Gene position index
        gene_position_start = -50
        # For -ve strand, this will be start of position as they go in opposite direction
        gene_position_end = gene_length + 50
        for k in dict_cds_info[gene_name]:
            for pos in range(k[0], k[1] + 1):
                if gene_position_start == 0:
                    gene_position_start += 1
                if gene_position_end == 0:
                    gene_position_end += -1
                try:
                    nuc = genome_dict[chr_num][pos - 1]
                except KeyError:
                    print 'KEYERROR in genome dict for chromosome ' + str(chr_num)
                    print 'Genome dict is ' + str(genome_dict.keys())
                except IndexError:
                    # A gene maybe present at the extreme end of chromosome (in E.coli) and hence 50 positions downstream will not be present
                    nuc = 'N'

                try:
                    counts_list = sam_parsed_count_dict[chr_num][pos][0:frag_range * 2]
                except KeyError:
                    # sam_parsed_count_dict is initialized only for genes with read counts.
                    # Positions with zero reads will lead to KeyErrors and we assign zero value to this position
                    counts_list = [0] * (frag_range * 2)
                # check if this position has multiple aligned reads
                if chr_num not in dict_mul_count:
                    pass
                elif pos in dict_mul_count[chr_num]:
                    mul_count_list = dict_mul_count[chr_num][pos][0:frag_range * 2]

                if strand == '+':
                    try:
                        multi_mapped = sam_parsed_count_dict[chr_num][pos][-2]
                    except KeyError:
                        multi_mapped = 0
                    # Get the read counts on the positive strand
                    strand_counts = countstrand(counts_list, '+')
                    if multi_mapped > 0:
                        # If there are multiple mapped reads at this position, get the read counts for those multiple mapped reads
                        mul_strand_counts = countstrand(mul_count_list, '+')
                    # For positive strand, the current position in "for" loop is from gene start side
                    current_pos = gene_position_start

                if strand == '-':
                    try:
                        multi_mapped = sam_parsed_count_dict[chr_num][pos][-1]
                    except KeyError:
                        multi_mapped = 0

                    # Get the read counts on the negative strand
                    strand_counts = countstrand(counts_list, '-')
                    if multi_mapped > 0:
                        mul_strand_counts = countstrand(mul_count_list, '-')
                    # Since the strand is negative, we take the complement of the current nucleotide
                    nuc = complimentarydict[nuc]
                    current_pos = gene_position_end

                # We will create dictionaries with fragment length as keys and dict of gene with list of read counts at every position as our values
                for fsize in range(frag_min, frag_max + 1):
                    if gene_name not in dict_count_len[fsize]:
                        dict_count_len[fsize][gene_name] = []
                    dict_count_len[fsize][gene_name].append(strand_counts[fsize - frag_min])

                outfile.write(gene_name + '\t' + str(current_pos) + '\t' + '\t'.join(map(str, [strand, gene_length, nuc, chr_num, pos,
                                                                                               multi_mapped, multi_genes] + strand_counts)) + '\n')

                if multi_mapped > 0:
                    # As done for actual counts, multiple mapped reads are counted according to fragment length
                    for fsize in range(frag_min, frag_max + 1):
                        if gene_name not in dict_mul_count_len:
                            dict_mul_count_len[gene_name] = {}
                        if fsize not in dict_mul_count_len[gene_name]:
                            dict_mul_count_len[gene_name][fsize] = {0: 0, 1: 0, 2: 0}
                        if current_pos < 0:
                            frame = current_pos % 3
                        else:
                            frame = current_pos % 3 - 1
                            if frame == -1:
                                frame = 2
                        if three_prime:
                            pos_range = range(1, gene_length + frag_max)
                        else:
                            pos_range = range(-frag_max, gene_length)
                        if current_pos in pos_range:
                            dict_mul_count_len[gene_name][fsize][frame] += mul_strand_counts[fsize - frag_min]
                if strand == '-':
                    gene_position_end += -1
                if strand == '+':
                    gene_position_start += 1
        counter += 1
        sys.stdout.write("Out of " + str(len(dict_cds_count)) + " transcripts, currently processing transcript {0}.\t\r".format(counter))
        sys.stdout.flush()

    outfile.close()

    # Determine the percentage of reads which are multiple mapped to a gene and discard it if it is greater than the threshold set for multiple map filter.
    # This is done specific to each fragment size
    mul_out_file = open(output+'Multiple_mapped_gene_read_counts_' + str(frag_min) + '_' + str(frag_max) + '.tab', 'w')

    for gene in dict_mul_count_len:
        mul_out_file.write(gene)
        for fsize in range(frag_min, frag_max + 1):
            for frame in xrange(3):
                mul_out_file.write('\t' + str(dict_mul_count_len[gene][fsize][frame]))
        mul_out_file.write('\n')
    mul_out_file.close()

    for fsize in range(frag_min, frag_max + 1):
        count_file = open(output+'Read_counts_' + str(fsize) + '.tab', 'w')
        for gene, reads_list in dict_count_len[fsize].iteritems():
            chr_num, strand = dict_gene[gene]
            length = dict_len[gene]
            if remove_overlapped_genes and gene in overlap_genes:
                continue
            if strand == '-':
                # For negative strand, read count were appended from the opposite end
                reads_list = reversed(reads_list)
            count_file.write(gene + '\t-50\t' + str(length) + '\t' + ','.join(map(str, reads_list)) + '\n')
        count_file.close()


# Calculate total length of multiCDS region
def multicdslen(cds_positions):
    length_total = 0
    for each_pos in cds_positions:
        each_len = each_pos[1] - each_pos[0] + 1
        length_total += each_len
    return length_total


# Count the reads on 5' end with same strand
def countstrand(counts_list, strand):
    # print 'Length of counts_list/2 is ' + str(len(counts_list)/2)
    frag_range = int(len(counts_list) / 2)
    strand_counts = []
    # initialize the strand count list
    for q in range(0, frag_range):
        strand_counts.append(0)

    # Start filling the strand count list with actual values from the count_list
    for r in range(0, frag_range):
        if strand == '+':
            strand_counts[r] = counts_list[r * 2]
        if strand == '-':
            strand_counts[r] = counts_list[r * 2 + 1]
    return strand_counts


# Parse the genome reference file
def read_fasta(genome_file):
    chrom, sequence = None, []
    for line in genome_file:
        line = line.rstrip()  # Removing any trailing whitespace characters
        if line.startswith(">"):
            if chrom:
                yield (chrom, ''.join(sequence))
            chrom, sequence = line, []
        else:
            sequence.append(line)  # Appends each sequence line to list. Finally joins in at the end of chromosome
    if chrom:
        yield (chrom, ''.join(sequence))


# Parse the GFF file to get the positions and no. of CDS regions in a gene.
# Important data structures: dictcdscount and dictcdsinfo
def process_gff(gff):
    dictgene = {}  # Contains gene names as keys and their start and end positions as a list of values
    dictcdscount = {}  # Contains gene names as keys and no. of cds regions(exons) for that particular gene
    dictcdsinfo = {}  # Contains the start and end positions of the cds regions as list of lists
    dict_len = {}

    with open(gff) as gff_file:
        for line in gff_file:
            line_list = line.strip().split('\t')
            if line.startswith('chr'):
                left_pos = int(line_list[3])
                right_pos = int(line_list[4])
                chrnum = line_list[0]
                strand = line_list[6]
                if line_list[2] == 'CDS':
                    cds_annotation = line_list[8]
                    # This is for sacCer3 genome annotation file
                    cds_name = cds_annotation.split(';')[0].split('=')[1].split('_')[0]
                    # cds_name = cds_annotation.split(';')[0].split('=')[1]
                    if cds_name not in dictcdscount:
                        dictcdscount[cds_name] = 1
                        dictcdsinfo[cds_name] = []
                    else:
                        dictcdscount[cds_name] += 1
                    dictgene[cds_name] = [chrnum, strand]
                    dictcdsinfo[cds_name].append([left_pos, right_pos])
    # intron_genes_file = open("list_genes_introns.tab", 'w')
    cds_file = open("cds_info.tab", "w")

    for gene in dictcdsinfo:
        chrnum, strand = dictgene[gene]
        dict_len[gene] = multicdslen(dictcdsinfo[gene])
        #  Gene name  Chr number Strand  Length of gene (CDS regions only) No of CDS regions
        cds_file.write(gene + '\t' + chrnum + '\t' + strand + '\t' + str(dict_len[gene]) + '\t' + str(dictcdscount[gene]))
        # Start and end of each CDS region
        for cds in dictcdsinfo[gene]:
            cds_file.write('\t' + str(cds[0]) + '\t' + str(cds[1]))
        cds_file.write('\n')

    # for gene in dictcdscount:
    #     if dictcdscount[gene] > 1:
    #         intron_genes_file.write(gene + '\n')
    gff_file.close()
    return cds_file


def process_gff_ecoli(gff):
    dictgene = {}  # Contains gene names as keys and their start and end positions as a list of values
    dictcdsinfo = {}  # Contains the start and end positions of the cds regions as list of lists
    dict_psuedo = {}
    dict_len = {}
    problem_genes = []
    chrom = ''
    with open(gff) as gff_file:
        for line in gff_file:
            line_list = line.strip().split('\t')
            if not line.startswith('##'):
                try:
                    left_pos = int(line_list[3])
                    right_pos = int(line_list[4])
                    strand = line_list[6]
                    if left_pos == 3698003:
                        print left_pos, right_pos, strand
                    chrom = line_list[0]
                    if line_list[2] == 'gene':
                        gene_anno = line_list[8]
                        fields = gene_anno.split(';')
                        gene_info = {}
                        for field in fields:
                            fd = field.split('=')
                            gene_info[fd[0]] = fd[1]
                        gene_name = gene_info['ID']
                        if 'pseudo' in gene_info:
                            if gene_info['pseudo'] == 'true':
                                # print 'Got pseudo gene '+gene_name
                                dict_psuedo[gene_name] = [left_pos, right_pos, strand]
                                continue
                            else:
                                dictgene[gene_name] = [left_pos, right_pos, gene_info]
                        else:
                            dictgene[gene_name] = [left_pos, right_pos, gene_info]
                    if line_list[2] == 'CDS':
                        cds_annotation = line_list[8]
                        fields = cds_annotation.split(';')
                        cds_info = {}
                        for field in fields:
                            fd = field.split('=')
                            cds_info[fd[0]] = fd[1]
                        gene_alias = cds_info['Parent']
                        if gene_alias in dict_psuedo or 'pseudo' in cds_info:
                            continue
                        if gene_alias in dictgene:
                            cds_name = dictgene[gene_alias][2]['locus_tag']

                        if gene_alias in dictgene:
                            if left_pos != dictgene[gene_alias][0]:
                                print 'CDS for gene ' + cds_name + ' does not start at gene start. The cds start is at ' + str(left_pos) + ' and gene start is at ' + str(
                                    dictgene[gene_alias][0])
                                print 'CDS for gene ' + cds_name + ' does not stop at gene stop. The cds stop is at ' + str(right_pos) + ' and gene stop is at ' + str(
                                    dictgene[gene_alias][1])
                                problem_genes.append(cds_name)
                        else:
                            print 'Gene with only CDS annotation is ', cds_name
                        if cds_name not in dictcdsinfo:
                            dictcdsinfo[cds_name] = [left_pos, right_pos, strand]
                        else:
                            print 'Second CDS present for gene ' + cds_name
                            print 'First CDS'
                            print dictcdsinfo[cds_name]
                            print 'Second CDS'
                            print left_pos, right_pos, strand
                            dictcdsinfo[cds_name] = [left_pos, right_pos, strand]
                        gene_length = right_pos - left_pos + 1
                        dict_len[cds_name] = gene_length
                except IndexError:
                    print 'Weird line:', line
    outfile = open("../data_files/ecoli/CDS_info.tab", 'w')
    for gene in problem_genes:
        del dictcdsinfo[gene]
    for gene in dictcdsinfo:
        left_pos, right_pos, strand = dictcdsinfo[gene]
        outfile.write(gene + '\t' + str(chrom) + '\t' + str(strand) + '\t' + str(dict_len[gene]) + '\t1\t' + str(left_pos) + '\t' + str(right_pos) + '\n')

    return outfile


def cdsparser(annot_file, genome, extra_overlap=0):
    dict_cds_count = {}
    dict_cds_info = {}
    dict_len = {}
    dict_gene_info = {}
    # print 'annotation file type ' + str(type(annot_file))
    with open(annot_file) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            gene, chrnum, strand, length, count = fields[:5]
            dict_cds_count[gene] = int(count)
            dict_len[gene] = int(length)
            dict_gene_info[gene] = [chrnum, strand]
            dict_cds_info[gene] = []
            for i in range(5, len(fields), 2):
                dict_cds_info[gene].append([int(fields[i]), int(fields[i + 1])])

    genomedict = {}  # This dictionary contains chromosomes (chrI, chrII, etc) as keys and the sequence as its value
    with open(genome) as fp:
        for name, seq in read_fasta(fp):
            genomedict[name[1::]] = seq
    # Get the list of genes which have overlapping cds regions
    overlap_genes = find_overlapping_genes(dict_cds_info, dict_gene_info, extra_len=extra_overlap)
    # print 'dictcdscount length ' + str(len(dict_cds_count))
    return dict_gene_info, dict_cds_count, dict_cds_info, genomedict, overlap_genes, dict_len


# Find overlapping genes from a GFF annotation file
def find_overlapping_genes(dict_cds_info, dict_gene, extra_len=0):
    dict_search_overlap = {}
    dict_overlap = {}
    for gene_name in dict_cds_info:
        chrnum, strand = dict_gene[gene_name]
        # Left most position
        leftpos = int(dict_cds_info[gene_name][0][0])
        # Right most position
        rightpos = int(dict_cds_info[gene_name][-1][1])
        try:
            regions = dict_search_overlap[chrnum]
            for gene in regions:
                if strand == gene[3]:
                    if leftpos - extra_len < gene[2] and rightpos + extra_len > gene[1]:
                        dict_overlap[gene[0]] = [chrnum, gene[1], gene[2], strand, gene_name, leftpos, rightpos, strand]
                        dict_overlap[gene_name] = [chrnum, leftpos, rightpos, strand, gene[0], gene[1], gene[2], strand]
            dict_search_overlap[chrnum].append([gene_name, leftpos, rightpos, strand])
        except KeyError:
            dict_search_overlap[chrnum] = [[gene_name, leftpos, rightpos, strand]]

    # overlap_file = open('Overlap_genes', 'w')
    overlap_list = []
    for gene, info in dict_overlap.iteritems():
        # Write the overlapping genes to a file
        # overlap_file.write(gene + '\t' + '\t'.join(map(str, info)) + '\n')
        overlap_list.append(gene)
    overlap_list.sort()

    # overlap_file.close()
    return overlap_list


def sort_cds_table(output, cds_table):
    dict_gene = {}
    dict_len = {}
    with open(cds_table) as f:
        for lines in f:
            fields = lines.strip().split('\t')
            gene = fields[0]
            if gene not in dict_gene:
                dict_gene[gene] = {}
            pos = int(fields[1])
            length = int(fields[3])
            vals = fields[2:]
            dict_gene[gene][pos] = vals
            if gene not in dict_len:
                dict_len[gene] = length

    out_file = open(output+'CDS_count_table.tab', 'w')
    for gene in sorted(dict_gene):
        for pos in sorted(dict_gene[gene]):
            out_file.write(gene + '\t' + str(pos) + '\t')
            try:
                out_file.write('\t'.join(map(str, dict_gene[gene][pos])) + '\n')
            except TypeError:
                print dict_gene[gene][pos]
    subprocess.call('rm ' + cds_table, shell=True)


def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-i", "--input",
                      dest="inputPath",
                      help="The path of the input SAM/BAM file.")
    parser.add_option("-a", "--alignment",
                      dest="alignment",
                      help="If aligned to a genome, enter genome(Default). If aligned to a transcriptome, enter 'transcriptome'.")
    parser.add_option("-g", "--annotation",
                      dest="annotation_file",
                      help="The path of the annotation file. If the alignment is done to a genome, the annotation file should be of following format."
                           "Gene\tChromosome number\tStrand\tNo. of CDS regions\tCDS_start1\tCDS_end1\tCDS_start2\tCDS_end2\t..."
                           "If the alignment is done to a transcriptome, the annotation file should be of following format"
                           "Gene\tStart Index\tLength of CDS")
    parser.add_option("-e", "--sequence",
                      dest="sequence",
                      help="The path of the sequence file. If alignment is done to a genome, the genome fasta file. "
                           "If alignment is done to a transcriptome, the fasta files of all the gene sequences aligned to.")
    parser.add_option("-m", "--min",
                      dest="min",
                      help="The minimum fragment length for which mapped reads are counted.")
    parser.add_option("-x", "--max",
                      dest="max",
                      help="The max fragment length for which mapped reads are counted.")
    parser.add_option("-3", "--three_prime",
                      dest="thr_prime",
                      help=" Yes if reads are to be quantified by 3' end instead of 5' end(Default).")
    parser.add_option("-o", '--output',
                      dest="output",
                      help="Output folder where the resulting files will be stored")
    parser.add_option("-v", "--overlap",
                      dest="overlap",
                      help="Overlapping genes are removed to avoid ambiguity in assignment of reads. This parameter is number of nucleotides beyond CDS region on both sides "
                           "to be avoid overlap. Default = 0")
    (options, arg) = parser.parse_args()
    if not options.inputPath:
        parser.error("Requires an input file, the path of this script. Run with --help for help.")
        quit(1)
    return options


# ===Program Start===
if __name__ == "__main__":
    # Parse arguments from command line
    print 'Starting to parse arguments'
    arguments = parse_arguments()

    if arguments.min:
        min_frag = int(arguments.min)
    else:
        min_frag = 20  # If no minimum fragment length is given, the default is 20
    if arguments.max:
        max_frag = int(arguments.max)
    else:
        max_frag = 35

    if arguments.output:
        output_dir = arguments.output
    else:
        output_dir = os.getcwd()
    if output_dir[-1] != '/':
        output_dir += '/'
    if not os.path.exists(output_dir+'output/'):
        os.makedirs(output_dir+'output/')
        out = output_dir+'output/'
    SAM = arguments.inputPath
    # If the alignment format is BAM, then we first convert to SAM format. Requires samtools to be installed.
    if SAM.split('.')[-1] == 'bam':
        subprocess.call('samtools view ' + SAM + ' > ' + SAM.split('.')[0] + '.sam', shell=True)
        SAM = SAM.split('.')[0] + '.sam'
    elif SAM.split('.')[-1] == "sam":
        pass
    else:
        print 'Unrecognized file format of alignment input file. Exiting...'
        sys.exit()
    if arguments.alignment:
        if arguments.alignment == 'transcriptome':
            transcriptome = True
        elif arguments.alignment == 'genome':
            transcriptome = False
        else:
            print 'Invalid alignment value entered. Exiting...'
            sys.exit()
    else:
        transcriptome = False
    if arguments.sequence:
        if transcriptome:
            transcripts = arguments.sequence
        else:
            genome_loc = arguments.sequence
    else:
        # Default alignment to an yeast sacCer3 genome
        genome_loc = '../data_files/sacCer3/sacCer3_R64-2-1_genome.fa'

    if arguments.annotation_file:
        annotations = arguments.annotation_file
    else:
        # Default annotations of yeast sacCer3
        annotations = '../data_files/sacCer3/sacCer3_R64-2-1_20150113.gff'
    print 'Parsed all arguments'
    if transcriptome:
        print 'Entering transcriptome mode'
        print 'Transcript sequence location: ' + transcripts
    else:
        print 'Entering genome mode '
        print 'Genome file location: ' + genome_loc

    if annotations.split('.')[-1] == 'gff' or annotations.split('.')[-1] == 'gff3':
        print '[Warning]: GFF file entered as input for gene annotations. Exercise caution. Parsing of GFF file is optimized only for sacCer3 and E.coli. Check the format'
        if 'sacCer3' in annotations or 'saccer3' in annotations:
            parsed_gff = process_gff(annotations)
        elif 'ecoli' in annotations:
            parsed_gff = process_gff_ecoli(annotations)
        print 'Got annotation file from processing the GFF file'
        annotations = parsed_gff.name
        print 'Annotation file: ' + annotations
    else:
        print 'Annotation file: ' + annotations

    if arguments.thr_prime:
        if arguments.thr_prime == 'Yes':
            thr_prime = True
        else:
            thr_prime = False
    else:
        thr_prime = False

    if arguments.overlap:
        overlap = int(arguments.overlap)
    else:
        overlap = 0

    if transcriptome:
        count_dict, mul_count_dict, total_dict = samparser_transcriptome(SAM, min_frag, max_frag)
        print 'Parsed the SAM file aligned to the transcriptome. Starting to make CDS table'
        create_cds_table_transcriptome(annotations, transcripts, out, count_dict, mul_count_dict, total_dict, min_frag, max_frag, filter_file='filter')
        print 'Created CDS table. Sorting...'
        sort_cds_table(out, out+'CDS_count_table')
        print 'Done'
    else:
        count_dict, mul_count_dict = samparser_genome(SAM, min_frag, max_frag, three_prime=thr_prime)
        print 'Parsed the SAM file. Starting to make CDS table'
        make_cds_table(annotations, genome_loc, out, count_dict, mul_count_dict, min_frag, max_frag, three_prime=thr_prime, extra_overlap=overlap)
        print 'Created CDS table. Sorting...'
        sort_cds_table(out, out+'CDS_count_table')
        print 'Done'
