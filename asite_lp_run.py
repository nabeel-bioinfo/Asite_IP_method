from __future__ import division
import cPickle as Pickle
import os
import sys
from collections import OrderedDict
from optparse import OptionParser
import numpy as np
import scipy.stats as stat

CODON_TYPES = ['UUU', 'UUC', 'UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG', 'AUU', 'AUC', 'AUA', 'AUG', 'GUU', 'GUC', 'GUA',
               'GUG', 'UCU', 'UCC', 'UCA', 'UCG', 'CCU', 'CCC', 'CCA', 'CCG', 'ACU', 'ACC', 'ACA', 'ACG', 'GCU', 'GCC',
               'GCA', 'GCG', 'UAU', 'UAC', 'CAU', 'CAC', 'CAA', 'CAG', 'AAU', 'AAC', 'AAA', 'AAG', 'GAU', 'GAC', 'GAA',
               'GAG', 'UGU', 'UGC', 'UGG', 'CGU', 'CGC', 'CGA', 'CGG', 'AGU', 'AGC', 'AGA', 'AGG', 'GGU', 'GGC', 'GGA',
               'GGG', 'UAA', 'UAG', 'UGA']

genetic_code = {'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',
                'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', 'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W',
                'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}

AMINO_ACIDS = ['A', 'R', 'D', 'N', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '*']


def parse_read_counts(count_directory, frag_min, frag_max, three_prime=False):
    if count_directory[-1] != '/':
        count_directory += '/'

    dict_reads = {}
    for fsize in range(frag_min, frag_max + 1):
        infile = open(count_directory + 'Read_counts_' + str(fsize) + '.tab', 'r')
        if fsize not in dict_reads:
            dict_reads[fsize] = {}
        for line in infile:
            fields = line.strip().split('\t')
            gene = fields[0]
            start_index = int(fields[1])
            cds_length = int(fields[2])
            reads = map(int, fields[3].split(','))
            if not three_prime:
                dict_reads[fsize][gene] = reads[abs(start_index + fsize):abs(start_index + fsize) + cds_length]
            else:
                dict_reads[fsize][gene] = reads[abs(start_index):abs(start_index) + cds_length + fsize]


# Using cutoffs for all possible frag sizes
def select_high_cov_genes(cds_folder, mul_map_file, frag_min, frag_max, output, threshold=1, three_prime=False):
    if cds_folder[-1] != '/':
        cds_folder += '/'

    # Parsing the cds file or loading it from the stored pickle dict
    dict_length = {}
    # get the read dict for every gene according to fragment size and at each gene position
    reads_dict = {}
    # A dict containing the total number of reads mapped to a gene
    total_reads = {}
    for fsize in range(frag_min, frag_max + 1):
        if fsize not in reads_dict:
            reads_dict[fsize] = {}
        with open(cds_folder + 'Read_counts_' + str(fsize) + '.tab') as f:
            for lines in f:
                line_list = lines.strip().split('\t')
                gene_name = line_list[0]
                if gene_name not in reads_dict[fsize]:
                    reads_dict[fsize][gene_name] = {}
                start_index = int(line_list[1])
                gene_len = int(line_list[2])
                if gene_name not in dict_length:
                    dict_length[gene_name] = gene_len
                count_list = map(int, line_list[3].split(','))
                if gene_name not in total_reads:
                    total_reads[gene_name] = 0
                idx = start_index
                for val in count_list:
                    if three_prime:
                        # For 3' end alignments, we include mapped reads upto fragment size length beyond the stop codon
                        if 1 <= idx <= gene_len + fsize:
                            reads_dict[fsize][gene_name][idx] = val
                            total_reads[gene_name] += val
                    else:
                        # We will count reads from fsize positions before the start position upto the last nt position of the gene
                        if -fsize <= idx <= gene_len:
                            reads_dict[fsize][gene_name][idx] = val
                            total_reads[gene_name] += val
                    idx += 1
                    if idx == 0:
                        idx += 1

        sys.stdout.write("Parsed read counts for fragment size {0}.\t\r".format(fsize))
        sys.stdout.flush()

    print '\nParsed read counts for all fragment sizes '
    # Get the number of mul mapped reads to decide whether to delete the gene or not. If a gene has more than 0.1% of reads multiple mapped, we delete it
    mul_map_dict = {}
    mul_map_gene_reads = {}
    with open(mul_map_file) as f:
        for lines in f:
            line_list = lines.strip().split('\t')
            gene_name = line_list[0]
            read_counts = map(int, line_list[1:])
            mul_map_gene_reads[gene_name] = sum(read_counts)
            idx = 0
            for fsize in range(frag_min, frag_max + 1):
                if fsize not in mul_map_dict:
                    mul_map_dict[fsize] = {}
                mul_map_dict[fsize][gene_name] = {0: read_counts[idx], 1: read_counts[idx + 1], 2: read_counts[idx + 2]}
                idx += 3

        mulfile = open(output + 'Genes_multiple_mapped_reads_agg.tab', 'w')
        # List of all genes which have > 1% mul mapped reads and hence will not be considered
        mul_map_genes = []
        for gene in mul_map_gene_reads:
            if mul_map_gene_reads[gene] > 0:
                if gene not in total_reads:
                    print gene + ' not in total reads. It must be overlapping gene'
                    continue
                else:
                    perc_mul_map = float(mul_map_gene_reads[gene]) * 100 / float(mul_map_gene_reads[gene] + total_reads[gene])
                    mulfile.write(gene + '\t' + str(perc_mul_map) + '%\t' + str(mul_map_gene_reads[gene]) + '\t' + str(total_reads[gene]) + '\n')
                    if perc_mul_map > 1:
                        mul_map_genes.append(gene)
    print '\nParsed the multiple mapped read counts\n'
    dict_gene = {}
    good_genes = {}
    good_genes_mul_map = {}

    log_file = open(output + 'select_genes.log', 'w')
    print 'Starting to select genes for each fragment size and frame'
    for fsize in range(frag_min, frag_max + 1):
        good_genes[fsize] = {}
        good_genes_mul_map[fsize] = {}
        for frame in xrange(3):
            good_genes[fsize][frame] = []
            good_genes_mul_map[fsize][frame] = []
        for gene_name, dict_reads in reads_dict[fsize].iteritems():
            if gene_name not in dict_gene:
                dict_gene[gene_name] = {}
                for size in range(frag_min, frag_max + 1):
                    dict_gene[gene_name][size] = {}
            reads = []

            try:
                short_utr = False
                # 5' end
                if not three_prime:
                    start = -fsize
                    end = dict_length[gene_name] + 1
                else:
                    start = 1
                    end = dict_length[gene_name] + fsize + 1
                for k in range(start, end):
                    # Ignoring since there is no zero position
                    if k == 0:
                        continue
                    try:
                        reads.append(dict_reads[int(k)])

                    except KeyError:
                        log_file.write(gene_name + '\t' + str(k) + '\tline 257\t' + str(dict_length[gene_name]))
                        short_utr = True
            except KeyError:
                print 'Length not available for gene ' + gene_name
                short_utr = False
            for frame in xrange(3):
                if three_prime:
                    ref = reads[frame::3]

                else:
                    # The extra number of nucleotides on one side of CDS
                    if short_utr:
                        extra_s = min(dict_reads)
                        if extra_s < 12:
                            continue
                    else:
                        extra_s = fsize
                    # and then choose the reads in positions of multiples of 3 according to the frame
                    ref = reads[extra_s % 3 + frame::3]
                    # select the reads of the given 5' end frame  # last_off is -fsize. That many positions will be left out from the right side of the list (after gene ending)
                if gene_name in mul_map_genes:
                    mul_map_ref = mul_map_dict[fsize][gene_name][frame]

                avg_reads = np.mean(ref)

                dict_gene[gene_name][fsize][frame] = (sum(ref), len(ref), avg_reads)
                if avg_reads > threshold:
                    if gene_name in mul_map_genes:
                        try:
                            perc_mul_map_fs_fr = float(mul_map_ref) * 100 / float(mul_map_ref + sum(ref))
                        except ZeroDivisionError:
                            perc_mul_map_fs_fr = 0
                        if perc_mul_map_fs_fr > 1:
                            good_genes_mul_map[fsize][frame].append(gene_name)
                        else:
                            good_genes[fsize][frame].append(gene_name)
                    else:
                        good_genes[fsize][frame].append(gene_name)

    # for fsize in sorted(good_genes):
    #     for frame in xrange(3):
    #         print 'Good genes for frag size %d and frame %d is %d with avg reads per nt greater than %f reads' % (fsize, frame, len(good_genes[fsize][frame]), threshold)

    outfile = open(output + 'Sorted_genes_by_Avg_reads_frag_size.tab', 'w')
    for name in dict_gene:
        outfile.write(name)
        for fsize in range(frag_min, frag_max + 1):
            try:
                for frame in xrange(3):
                    outfile.write('\t' + str(dict_gene[name][fsize][frame][0]) + '\t' + str(dict_gene[name][fsize][frame][1]) + '\t' + str(dict_gene[name][fsize][frame][2]))
            except KeyError:
                continue
        outfile.write('\n')

    # final_list = open("Final_selected_genes.tab", "w")
    filtered_cds_dict = {}
    for i in range(frag_min, frag_max + 1):
        filtered_cds_dict[i] = {}
        for frame in xrange(3):
            filtered_cds_dict[i][frame] = {}
    log_file = open(output + 'Gene_filter_statistics.tab', 'w')
    for fsize in filtered_cds_dict:
        for frame in xrange(3):
            for name in good_genes[fsize][frame]:
                # if name in ['YFR032C-A', 'YFR032C-B', 'YGL188C-A', 'YGL189C']:
                #    continue
                filtered_cds_dict[fsize][frame][name] = reads_dict[fsize][name]

            log_file.write(str(fsize) + '\t' + str(frame) + '\t' + str(len(filtered_cds_dict[fsize][frame])) + '\t' + str(len(good_genes_mul_map[fsize][frame])) + '\n')
    return filtered_cds_dict, dict_length


def asite_algorithm_improved_second_offset_correction(reads_dict, dict_len, frag_min, frag_max, output, offset_threshold=70, off_correction_threshold=5, get_coverage_stats=True,
                                                      advanced=True, conserve_frame=True, bootstrap=False, three_prime=False,
                                                      cov_range=(1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 25, 30, 35,
                                                                 40, 45, 50)):
    print 'Got gene lengths for ' + str(len(dict_len)) + ' genes'

    # Used for debugging
    offset_dic = {}
    correction_dic = {}
    details = {}
    debug_details = {}

    for frame in xrange(3):
        details[frame] = {}
        debug_details[frame] = {}
    log_file = open(output + 'asite_lp.log', 'w')

    # The following dict will contain meta gene for every gene in every fsize and frame. The meta-data include no. of zeros, perc zeros, avg, avg at start and avg at end
    sum_total = {}
    dict_cov_info = {}
    # Initialize files for all fragment size and frame combinations
    for fsize in range(frag_min, frag_max + 1):
        dict_cov_info[fsize] = {}
        sum_total[fsize] = {}
        for frame in xrange(3):
            dict_cov_info[fsize][frame] = {}
            sum_total[fsize][frame] = 0

    # MAIN ANALYSIS STARTS HERE
    for fsize in range(frag_min, frag_max + 1):
        print 'Starting analysis for fragment size ' + str(fsize)

        # The following will be used as an index in the read dictionary for each gene
        last_off = -fsize

        for frame in xrange(3):
            if frame not in offset_dic:
                offset_dic[frame] = {}
                correction_dic[frame] = {}
            if fsize not in offset_dic[frame]:
                offset_dic[frame][fsize] = {}
            '''
            FOLLOWING CODE FOR GENES IN EVERY FSIZE AND FRAME
            '''
            for gene in reads_dict[fsize][frame]:
                skip_gene = False
                short_utr = False

                """
                ***  PARSING THE READ VALUES   ***
                """
                # The following will be reads dictionary with nucleotide position as key and number of reads as value
                dict_reads = reads_dict[fsize][frame][gene]

                reads = []
                try:
                    # 5' end
                    if not three_prime:
                        start = -fsize
                        end = dict_len[gene] + 1
                    # For 3' end, we get the reads at the CDS nucleotide positions as well as fsize length after the end of CDS
                    else:
                        start = 1
                        end = dict_len[gene] + fsize + 1

                    for k in range(start, end):
                        # Ignoring since there is no zero position
                        if k == 0:
                            continue
                        try:
                            reads.append(dict_reads[k])
                        except KeyError:
                            log_file.write('ERROR: --Short UTR--KeyError in reads dictionary for ' + gene + ' at position ' + str(k) + 'with gene length of ' +
                                           str(dict_len[gene]) + '\n')
                            short_utr = True
                            if k > 0:
                                skip_gene = True
                                log_file.write('ERROR: KeyError in reads dictionary for ' + gene + ' at position ' + str(k) + 'with gene length of ' + str(dict_len[gene]) + '\n')
                                # Using break instead of continue as this will break the inner for loop and continue the outer for loop
                                break
                except KeyError:
                    # Length not available for this gene
                    skip_gene = True
                    log_file.write('ERROR: KeyError in length dict for ' + gene + '\n')

                if skip_gene:
                    # If gene does not have proper read values, we will not consider it for analysis and hence we remove it and move to next gene in the for loop
                    continue
                # The extra number of nucleotides on one side of CDS
                if short_utr:
                    extra_s = min(dict_reads)
                    if extra_s < 12:
                        continue
                else:
                    extra_s = fsize

                if three_prime:
                    ref = [0] * frame
                    for r in reads[frame::3]:
                        ref += [r, 0, 0]
                else:
                    # To adjust for the len and replace zeroes for out of frame reads, we do the following. This will be followed by deletion/addition of
                    # additional zeroes to make length equal to original seq
                    ref = [0] * (extra_s % 3 + frame)
                    # Choose the reads in positions of multiples of 3 according to the frame
                    for r in reads[extra_s % 3 + frame::3]:  # select the reads of the given 5' end frame
                        ref += [r, 0, 0]
                ref = ref[:-2]  # we exclude the last [0,0] we added at the end
                if (len(reads) - len(ref)) > 0:
                    ref += (len(reads[:last_off]) - len(ref)) * [0]
                # we put it back to the original length (which might have changed when multiple of 3).
                if three_prime:
                    avg_reads = np.mean(ref[frame::3])
                else:
                    avg_reads = np.mean(ref[extra_s % 3 + frame::3])
                """
                ***    CALCULATING THE OFFSET FOR THE GENE  ***
                """
                # try:
                if conserve_frame:
                    if three_prime:
                        best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf_3_end(ref, extra_s, dict_len[gene], centre_on_p_not_on_a=False,
                                                                                                            advanced=advanced,
                                                                                                            go_three_at_a_time=True)
                    else:
                        best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf(ref, extra_s, dict_len[gene], centre_on_p_not_on_a=False,
                                                                                                      advanced=advanced, go_three_at_a_time=True)
                else:
                    # we get the offset to centre on frame 1, hence +1 and -2. (-2 is because we use the length to get to the last index by start+len,
                    # but we gave +1 to start so we need to take out an extra one to len).
                    best_start_index, best_end_index, offset, score_per_offset = put_reads_in_orf(ref, extra_s + 1, dict_len[gene] - 2, centre_on_p_not_on_a=False,
                                                                                                  advanced=advanced, go_three_at_a_time=True)

                # If secondary selection criteria is to be applied, we compare the scores of top two offset and the reads near the start codon.
                if advanced:
                    # sort the offsets based on scores.  If scores are same for two or more offsets, they will be sorted according to offset values.
                    sorted_scores = sorted(sorted(score_per_offset), key=score_per_offset.get, reverse=True)
                    # Quality check to make sure the highest offset is the same as the best offset we get from our function
                    if sorted_scores[0] != offset:
                        log_file.write('ERROR: Sorted offsets do not match the offset we get from put_reads_in_orf function for gene ' + gene + ' in fragment size' +
                                       str(fsize) + 'and frame ' + str(frame) + '. The sorted scores are ' + str(sorted_scores) + ' and the offset is ' + str(offset) + '\n')
                        print gene, fsize, frame, sorted_scores[0], offset, sorted_scores

                    # Difference of top two offsets
                    diff = score_per_offset[offset] - score_per_offset[sorted_scores[1]]

                    # If the difference in scores is less than the avg number of reads across the gene, we apply the secondary selection criteria
                    if diff < avg_reads:
                        # Offsets wit diff less than avg will be listed in here
                        list_offsets_to_compare = []
                        # Add the top two offsets to the list of offsets to compare
                        if score_per_offset[sorted_scores[0]] >= score_per_offset[sorted_scores[1]]:
                            list_offsets_to_compare.append(sorted_scores[0])
                            list_offsets_to_compare.append(sorted_scores[1])
                            # If any other offsets down the order have equal scores with second best offset, then they get added to the list as well
                            for i in range(2, len(sorted_scores)):
                                if score_per_offset[sorted_scores[i]] == score_per_offset[sorted_scores[1]]:
                                    list_offsets_to_compare.append(sorted_scores[i])
                        # The offsets for which the condition is met will be added in here
                        off_true = []
                        # The dict will contain the difference between the average of R2, R3 and R4 and the reads in first codon R1
                        diff_dict = {}

                        # Check the secondary selection criteria of the listed offsets
                        for off in list_offsets_to_compare:
                            # quality check.
                            if off > fsize:
                                log_file.write('Unusual offset ' + str(off) + ' being considered for fsize ' + str(fsize) + ' frame ' + str(frame) + ' in gene ' + gene + '\n')
                                continue
                            # quality check
                            if off % 3 != 0:
                                log_file.write('Unusual offset ' + str(off) + ' being considered for fsize ' + str(fsize) + ' frame ' + str(frame) + ' in gene ' + gene + '\n')
                            # Getting the first 4 codon values in the particular offset
                            if three_prime:
                                reads_offset = reads[off:off + 12]
                            else:
                                reads_offset = reads[extra_s - off:extra_s - off + 12]
                            if not reads_offset:
                                print 'Reads offset list is empty'
                                print extra_s, off
                            # Checking the condition whether the R1 is less than one-fifth of the average of R2, R3 and R4
                            bool_off, diff_avg = secondary_selection_conditions(reads_offset, frame, threshold=off_correction_threshold)
                            # Adding this offset to the list if the condition is met
                            if bool_off:
                                off_true.append(off)
                                diff_dict[off] = diff_avg
                        # Select the offset which meets the secondary selection criteria
                        if len(off_true) == 1:
                            offset_correction = off_true[0]
                        # If more than one offset meets the secondary selection criteria, then choose the offset with the maximum score
                        elif len(off_true) > 1:
                            diff_compare = {}
                            # check if the scores are equal or not. If the scores are also equal, add the difference of avg(R2, R3, R4) and R1 to diff compare and
                            # choose the offset with max difference
                            max_score = score_per_offset[sorted_scores[0]]
                            for i in range(0, len(off_true)):
                                if score_per_offset[off_true[i]] == max_score:
                                    diff_compare[off_true[i]] = diff_dict[off_true[i]]
                            if len(diff_compare) > 1:
                                sorted_diff = sorted(diff_compare, key=diff_compare.get, reverse=True)
                                offset_correction = sorted_diff[0]
                            else:
                                # check if anything changes
                                # offset_correction = sorted_scores[0]
                                offset_correction = off_true[0]
                        # If no offset meets the secondary selection criteria, we let the offset with the maximum score remain the best offset.
                        # For offsets with equal scores, the smallest offset is the optimal offset
                        else:
                            offset_correction = sorted_scores[0]

                        ''' 
                        CHECK FOR OFFSET CORRECTION 
                        '''
                        # If the offset after applying secondary selection criteria is not the same as offset with top score, change the optimal offset
                        if offset_correction != offset:
                            offset = offset_correction
                        # If the offset is the same as the original one, no need to change the offset
                        elif offset_correction == offset:
                            pass
                        # If the offset correction did not yield any offset (secondary conditions not met for any offset), the initial offset remains the optimal offset
                        elif offset_correction == '':
                            pass

                # OFFSET IS ASSIGNED FOR THE GENE IN THIS FRAGMENT SIZE AND FRAME
                if not skip_gene:
                    offset_dic[frame][fsize][gene] = offset
                else:
                    continue

                if three_prime:
                    sum_total[fsize][frame] += sum(ref[frame::3])
                else:
                    sum_total[fsize][frame] += sum(ref[extra_s % 3 + frame::3])
                dict_cov_info[fsize][frame][gene] = avg_reads

    """
    ***   ANLAYSE DATA FOR DIFFERENT COVERAGE THRESHOLDS ***
    """
    if get_coverage_stats:
        # This dict will contain the bootstrapped distributions from which the error bars will be calculated
        bootstrap_dict = {}
        dict_most_prob_offsets = {}
        print 'Running the coverage analysis with offset threshold ' + str(offset_threshold) + '%'
        # Get coverage statistics
        for fsize in range(frag_min, frag_max + 1):
            # Needed to create the error bars for the plot for Figure 3.
            bootstrap_dict[fsize] = {}
            dict_most_prob_offsets[fsize] = {}
            for frame in xrange(3):
                bootstrap_dict[fsize][frame] = {}
                read_avg_dict = {}
                for off in range(0, fsize, 3):
                    bootstrap_dict[fsize][frame][off] = {}
                    read_avg_dict[off] = []
                dict_most_prob_offsets[fsize][frame] = {'off': '', 'perc': ''}
                gene_count = 0
                # Append the meta gene properties to their respective dictionaries according to the offsets
                for gene, offset in offset_dic[frame][fsize].iteritems():
                    try:
                        # Get the meta data of each gene from the dict_cov_info dictionary
                        read_avg = dict_cov_info[fsize][frame][gene]
                    except KeyError:
                        print 'KeyError for dict_cov_info in ' + str(fsize) + ' and frame ' + str(frame)
                        continue
                    read_avg_dict[offset].append(read_avg)
                    gene_count += 1
                list_genes = []  # Required for bootstrapping

                for off in range(0, fsize, 3):
                    for j in read_avg_dict[off]:
                        # this is appending the offset value with the read average
                        list_genes.append((off, j))
                if len(list_genes) == 0:
                    continue

                got_offset = False
                # First list will contain the cutoff values and the second list will contain the percentage values
                trend_list = [[], []]
                for c in cov_range:
                    sum_dict, perc_dict = count_stuff_all_offsets(c, read_avg_dict, False)
                    no_of_genes = sum(sum_dict.values())
                    sorted_perc = sorted(perc_dict, key=perc_dict.get, reverse=True)
                    if not got_offset:
                        if no_of_genes >= 10:
                            if perc_dict[sorted_perc[0]] >= offset_threshold:
                                dict_most_prob_offsets[fsize][frame]['off'] = str(sorted_perc[0])
                                dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                                next_highest_perc = perc_dict[sorted_perc[1]]
                                next_best_off = sorted_perc[1]
                                got_offset = True
                            else:
                                dict_most_prob_offsets[fsize][frame]['off'] = str(sorted_perc[0]) + '/' + str(sorted_perc[1])
                                dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]]) + '/' + str(perc_dict[sorted_perc[1]])

                        elif no_of_genes < 10 and dict_most_prob_offsets[fsize][frame]['off'] == '':
                            dict_most_prob_offsets[fsize][frame]['off'] = 'NA'
                            dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                    elif perc_dict[sorted_perc[0]] >= offset_threshold and no_of_genes >= 10:
                        dict_most_prob_offsets[fsize][frame]['perc'] = str(perc_dict[sorted_perc[0]])
                    if no_of_genes >= 10:
                        trend_list[0].append(c)
                        trend_list[1].append(perc_dict[sorted_perc[0]])

                # When the threshold is crossed at an higher cutoff, we do a linear trend analysis to check whether the trend is significant
                if trend_list[1] and trend_list[1][0] < offset_threshold and '/' not in dict_most_prob_offsets[fsize][frame]['off'] and 'NA' not in \
                        dict_most_prob_offsets[fsize][frame]['off']:
                    print 'Doing Linear trend analysis for fsize ' + str(fsize) + ' and frame ' + str(frame) + ' for trends for offset ' + \
                          str(dict_most_prob_offsets[fsize][frame]['off']) + ' with trends ' + str(trend_list)
                    b1, b0, r, pval, std_err = stat.linregress(trend_list[0], trend_list[1])

                    if pval > 0.05:
                        dict_most_prob_offsets[fsize][frame]['off'] += '/' + str(next_best_off)
                        dict_most_prob_offsets[fsize][frame]['perc'] += '/' + str(next_highest_perc) + '(NS, slope = ' + str(np.round(b1, 3)) + ' pval =' + str(
                            np.round(pval, 3)) + ')'
                    else:
                        dict_most_prob_offsets[fsize][frame]['perc'] += '(Significant, slope = ' + str(np.round(b1, 3)) + ' pval =' + str(np.round(pval, 3)) + ')'

                if bootstrap:
                    mega_sum_dict, mega_perc_dict = bootstrap_gene_count(c, list_genes)
                    Pickle.dump(mega_perc_dict, open('Pickle_dicts/bootstrap_perc_dict.p', 'wb'))
                    Pickle.dump(mega_sum_dict, open('Pickle_dicts/bootstrap_count_dict.p', 'wb'))
                    # Write the properties of the genes at each offset for eack threshold value
                    for off in range(0, fsize, 3):
                        if off not in mega_perc_dict:
                            print str(off) + ' was not in mega_perc_dict so adding a empty directory for fsize ' + str(fsize) + ' and frame ' + str(frame)
                            mega_perc_dict[off] = [0]
                        try:
                            avg = np.mean(mega_perc_dict[off])
                        except TypeError:
                            avg = 'TE'
                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['avg'] = avg
                    for off in range(0, fsize, 3):
                        try:
                            se_mean = np.std(mega_perc_dict[off])
                        except TypeError:
                            se_mean = 'TE'

                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['se_mean'] = se_mean
                    for off in range(0, fsize, 3):
                        try:
                            low_ci = np.percentile(mega_perc_dict[off], 2.5)
                        except TypeError:
                            low_ci = 'TE'

                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['low_ci'] = low_ci
                    for off in range(0, fsize, 3):
                        try:
                            high_ci = np.percentile(mega_perc_dict[off], 97.5)
                        except TypeError:
                            high_ci = 'TE'
                        if c == 1:
                            bootstrap_dict[fsize][frame][off]['high_ci'] = high_ci
                    for off in range(0, fsize, 3):
                        if off not in mega_sum_dict:
                            print str(off) + ' was not in mega_sum_dict so adding a empty directory for fsize ' + str(fsize) + ' and frame ' + str(frame)
                            mega_sum_dict[off] = [0]

                        try:
                            # the average of the bootstrapped distribution should be close to the original value
                            avg = np.mean(mega_sum_dict[off])
                        except TypeError:
                            avg = 'TE'

                    for off in range(0, fsize, 3):
                        try:
                            # The standard deviation of bootstrap distribution is the Standard error of mean of the original value
                            se_mean = np.std(mega_sum_dict[off])
                        except TypeError:
                            se_mean = 'TE'

                    for off in range(0, fsize, 3):
                        try:
                            low_ci = np.percentile(mega_sum_dict[off], 2.5)
                        except TypeError:
                            low_ci = 'TE'

                    for off in range(0, fsize, 3):
                        try:
                            high_ci = np.percentile(mega_sum_dict[off], 97.5)
                        except TypeError:
                            high_ci = 'TE'

    """
    ***    WRITE THE RESULTS AND PLOT DISTRIBUTION OF OFFSETS ***
    """
    outfile = open(output + "Results_LP_algorithm.tab", "w")
    perc_file = open(output + "Perc_of_genes_for_all_offsets.tab", "w")

    outfile.write('\n\nMost probable Offsets for Fragment Size and Frame (including coverage data)\n')
    outfile.write('Frag size\tFrame_0\tFrame_1\tFrame_2\n')
    for fsize in range(frag_min, frag_max + 1):
        try:
            outfile.write(str(fsize) + '\t' + str(dict_most_prob_offsets[fsize][0]['off']) + '\t' + str(dict_most_prob_offsets[fsize][1]['off']) + '\t' + str(
                dict_most_prob_offsets[fsize][2]['off']) + '\n')
        except KeyError:
            outfile.write(str(fsize) + '\tNA\tNA\tNA\n')

    outfile.write('\n\nPerc of genes (including coverage data)\n')
    outfile.write('Frag size\tFrame_0\tFrame_1\tFrame_2\n')
    for fsize in range(frag_min, frag_max + 1):
        try:
            outfile.write(str(fsize) + '\t' + str(dict_most_prob_offsets[fsize][0]['perc']) + '\t' + str(dict_most_prob_offsets[fsize][1]['perc']) + '\t' + str(
                dict_most_prob_offsets[fsize][2]['perc']) + '\n')
        except KeyError:
            outfile.write(str(fsize) + '\tNA\tNA\tNA\n')

    outfile.write('\n\n\tNumber of genes\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        outfile.write(str(fsize))
        for frame in xrange(3):
            try:
                outfile.write('\t' + str(len(offset_dic[frame][fsize])))
            except KeyError:
                outfile.write('\tNA')
        outfile.write('\n')

    outfile.write('\n\n\tNumber of reads\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        outfile.write(str(fsize))
        for frame in xrange(3):
            outfile.write('\t' + str(sum_total[fsize][frame]))
        outfile.write('\n')

    perc_file.write('Percentage of genes\nFrag/Frame\t0\t1\t2\n')
    for fsize in range(frag_min, frag_max + 1):
        perc_file.write(str(fsize))
        for frame in xrange(3):
            offset_list = {}
            for off in range(0, fsize, 3):
                offset_list[off] = []
            try:
                for gene, val in offset_dic[frame][fsize].iteritems():
                    if val in offset_list:
                        offset_list[val].append(gene)
                for off in sorted(offset_list):
                    try:
                        perc_file.write('\t' + str(float(len(offset_list[off])) * 100 / float(len(offset_dic[frame][fsize]))))
                    except KeyError:
                        perc_file.write('\tNA')
                    except ZeroDivisionError:
                        perc_file.write('\tNA')
            except KeyError:
                pass
            perc_file.write('\t')
        perc_file.write('\n')
    if bootstrap:
        perc_file.write('Perc of gene bootstrap results (AVERAGE)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')

        for fsize in range(frag_min, frag_max + 1):
            outfile.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['avg']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

        perc_file.write('Perc of gene bootstrap results (SE)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['se_mean']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')
        perc_file.write('Perc of gene bootstrap results (LOW CI)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['low_ci']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

        perc_file.write('Perc of gene bootstrap results (HIGH CI)\nFrag\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\t0\t3\t6\t9\t12\t15\t18\t21\n')
        for fsize in range(frag_min, frag_max + 1):
            perc_file.write(str(fsize))
            for frame in xrange(3):
                try:
                    for off in [0, 3, 6, 9, 12, 15, 18, 21]:
                        perc_file.write('\t' + str(bootstrap_dict[fsize][frame][off]['high_ci']))
                except KeyError:
                    perc_file.write('\tNA')
                perc_file.write('\t')
            perc_file.write('\n')

    return dict_most_prob_offsets


def count_stuff_all_offsets(count, dict_count, less=True):
    sum_dict = {}
    perc_dict = {}

    # for offset in [0,3,6,9,12,15,18,21]:

    for offset in dict_count:

        if less:
            sum_dict[offset] = sum(i < count for i in dict_count[offset])
        else:
            sum_dict[offset] = sum(i > count for i in dict_count[offset])

    for offset in dict_count:
        try:
            perc_dict[offset] = np.round(float(sum_dict[offset]) * 100 / float(sum(sum_dict.values())), 2)

        except ZeroDivisionError:
            perc_dict[offset] = 'NA'

    return sum_dict, perc_dict


#  Improved check of secondary offset correction.  Just checks the condition, whether it is valid or not
def secondary_selection_conditions(reads, frame, threshold=5):
    try:
        avg_three_codons = float(reads[3 + frame] + reads[6 + frame] + reads[9 + frame]) / 3
    except IndexError:
        print 'IndexError in calculating average of R2,R3,R4'
        print reads
        print frame
        avg_three_codons = 0
    bool_first = False

    # The average of 2nd, 3rd and 4th codon should be greater than reads in first codon and reads in second codon must be greater in third and fourth codons
    if (avg_three_codons > threshold * reads[frame]) and (reads[3 + frame] > reads[6 + frame]):
        bool_first = True
    diff = avg_three_codons - reads[frame]

    return bool_first, diff


def put_reads_in_orf(reads, gene_start_index, gene_len, centre_on_p_not_on_a=False, advanced=True, go_three_at_a_time=False):
    """
    tries to put most of the reads in the ORF.
    Returns best_start_index, best_end_index, offset (the indices are in reads, offset is how much - in nucleotides - one should move the reads to be aligned with the A site).
     Note that if centre_on_P_not_on_A is False then best_end_index-best_start_index+1 == gene_len-3 since we exclude the first codon when we centre on the A site
     (translation initiate at with AUG in the P site).
     ASSUMES missing reads value is zero
     Assumes the offset is positive (it never moves backward.. reverse genes should be 'straightened beforehand (this happens when they are annotated with reads from parsing of
     SAM file).

    Parameters
    ----------
    go_three_at_a_time
    advanced
    centre_on_p_not_on_a
    gene_len
    gene_start_index
    reads
    """
    if not centre_on_p_not_on_a:
        # Start from second codon
        gene_start_index += 3
        # The length of gene will be counted from second codon until stop codon
        gene_len -= 3
    off = gene_start_index
    # It looks like iteration will happen in opposite direction Starting from the fragment size value+3 (since we are dealing from second codon).
    # However the offset is 0 below even though we are looking at pos with index off and off+gene_len. We are using pos indexes to get sum of the reads while offset is 0.
    reads_in_cds = sum(reads[off: off + gene_len])  # this is for the initial guessed start_index equal to gene_start_index
    score = reads_in_cds
    best_offset = off
    score_per_offset = None
    if advanced:
        score_per_offset = OrderedDict()
        # gene_start_index-off is our actual offset and now we create a dict of scores with offsets as keys
        score_per_offset[gene_start_index - off] = reads_in_cds
    # We are increasing the offset by 3 by shifting the pos indexes by -3 positions to calculate sum of reads
    if go_three_at_a_time:
        off -= 3
    else:
        off -= 1
    while off > 3:
        if go_three_at_a_time:
            reads_in_cds = reads_in_cds + sum(reads[off:off + 3]) - sum(reads[off + gene_len:off + 3 + gene_len])
        else:
            reads_in_cds = reads_in_cds + reads[off] - reads[off + gene_len]
        if advanced:
            score_per_offset[gene_start_index - off] = reads_in_cds
        if reads_in_cds > score:
            # If the read sum is greater than the best score, make it the best score and the offset
            score = reads_in_cds
            best_offset = off
        if go_three_at_a_time:
            off -= 3
        else:
            off -= 1

    return best_offset, best_offset + gene_len - 1, gene_start_index - best_offset, score_per_offset


def put_reads_in_orf_3_end(reads, fsize, gene_len, centre_on_p_not_on_a=False, advanced=True, go_three_at_a_time=False):
    """
    tries to put most of the reads in the ORF.
    Returns best_start_index, best_end_index, offset (the indices are in reads, offset is how much - in nucleotides - one should move the reads to be aligned with the A site).
     Note that if centre_on_P_not_on_A is False then best_end_index-best_start_index+1 == gene_len-3 since we exclude the first codon when we centre on the A site
     (translation initiate at with AUG in the P site).
     ASSUMES missing reads value is zero
     Assumes the offset is positive (it never moves backward.. reverse genes should be 'straightened beforehand (this happens when they are annotated with reads from parsing of
     SAM file).

    Parameters
    ----------
    go_three_at_a_time
    advanced
    centre_on_p_not_on_a
    gene_len
    fsize
    reads
    """
    if not centre_on_p_not_on_a:
        # Start from second codon
        gene_start_index = 3
        # The length of gene will be counted from second codon until stop codon
        gene_len -= 3
    off = gene_start_index
    # Iteration will happen in opposite direction. Starting from the fragment size value+3 (since we are dealing from second codon).
    # However the offset is 0 below even though we are looking at pos with index off and off+gene_len. We are using pos indexes to get sum of the reads while offset is 0.
    reads_in_cds = sum(reads[off: off + gene_len])  # this is for the initial guessed start_index equal to gene_start_index
    score = reads_in_cds
    best_offset = off
    score_per_offset = None
    if advanced:
        score_per_offset = OrderedDict()
        # gene_start_index-off is our actual offset and now we create a dict of scores with offsets as keys
        score_per_offset[gene_start_index - off] = reads_in_cds
    # We are increasing the offset by 3 by shifting the pos indexes by -3 positions to calculate sum of reads
    # if go_three_at_a_time:
    #     off += 3
    # else:
    #     off += 1
    while off < fsize:
        if go_three_at_a_time:
            reads_in_cds = reads_in_cds - sum(reads[off:off + 3]) + sum(reads[off + gene_len:off + 3 + gene_len])
        else:
            reads_in_cds = reads_in_cds + reads[off] - reads[off + gene_len]
        if go_three_at_a_time:
            off += 3
        else:
            off += 1
        if advanced:
            score_per_offset[off - gene_start_index] = reads_in_cds
        if reads_in_cds > score:
            # If the read sum is greater than the best score, make it the best score and the offset
            score = reads_in_cds
            best_offset = off

    return best_offset, best_offset + gene_len - 1, best_offset - gene_start_index, score_per_offset


def bootstrap_gene_count(cutoff, list_genes):
    mega_sum_dict = {}
    mega_perc_dict = {}
    for itr in range(1, 10001):
        # For each iteration, declare a random dict
        rand_dict = {}
        # Select randomly with replacement the indexes of the actual dict
        list_rand_idx = list(np.random.choice(range(0, len(list_genes)), size=len(list_genes)))
        # Populate the random dict with the actual value of offset and the read avg
        for l in list_rand_idx:
            offset, value = list_genes[l]
            if offset not in rand_dict:
                rand_dict[offset] = []
            # Append the gene value
            rand_dict[offset].append(value)
        # Calculate the no of genes for each offset along with the percentage of genes according to offsets for the total genes at this threshold
        sum_dict, perc_dict = count_stuff_all_offsets(cutoff, rand_dict, less=False)
        for off in [0, 3, 6, 9, 12, 15, 18, 21]:
            if off not in mega_sum_dict:
                mega_sum_dict[off] = []
                mega_perc_dict[off] = []
            try:
                mega_sum_dict[off].append(sum_dict[off])
                mega_perc_dict[off].append(perc_dict[off])
            except KeyError:
                mega_sum_dict[off].append(0)
                mega_perc_dict[off].append(0)

    return mega_sum_dict, mega_perc_dict


def parse_arguments():
    parser = OptionParser("Usage: %prog -i <inputpath>")
    parser.add_option("-i", "--input",
                      dest="inputPath",
                      help="The full absolute path of the CDS file.")
    parser.add_option("-j", "--input2",
                      dest="inputPath2",
                      help="The full absolute path of the Mul mapped file.")
    parser.add_option("-m", "--min",
                      dest="min",
                      help="The minimum fragment length for which mapped reads are counted.")
    parser.add_option("-x", "--max",
                      dest="max",
                      help="The max fragment length for which mapped reads are counted.")
    parser.add_option("-o", "--output",
                      dest="output",
                      help="Output folder.")
    parser.add_option("-n", "--number",
                      dest="number",
                      help="Threshold of avg reads per codon for filtering genes")
    parser.add_option("-a", "--method",
                      dest="method",
                      help="The A-site method to be compared with LP method.")
    parser.add_option("-d", "--dataset",
                      dest="dataset",
                      help="The A-site method to be compared with LP method.")
    parser.add_option("-f", "--frac",
                      dest="frac_threshold",
                      help="Threshold for filtering genes with less than frac of positions with zero reads")
    parser.add_option("-p", "--xpp",
                      dest="xpp",
                      help="Boolean. For False (default), carry out analysis for PPX. For True, carry out analysis for XPP")
    parser.add_option("-t", "--threshold",
                      dest="threshold",
                      help="Threshold for assigning a unique offset for a combination of fragment size and frame")
    parser.add_option("-s", "--second_threshold",
                      dest="second_threshold",
                      help="Threshold for secondary offset correction. No. of times the reads in start codon should be less than times the avg of 2nd, 3rd and 4th codons")
    parser.add_option("-3", "--three-prime",
                      dest="three_prime",
                      help="Align reads by 3' end. Default(5' end)")
    (options, arg) = parser.parse_args()
    # if not options.inputPath:
    #     parser.error("Requires an input file, the path of this script. Run with --help for help.")
    #     quit(1)
    return options


# ===Program Start===
if __name__ == "__main__":
    # Parse arguments from commandline
    print 'Starting to parse arguments'
    arguments = parse_arguments()

    if arguments.output:
        output_dir = arguments.output
    else:
        output_dir = os.getcwd()
    if output_dir[-1] != '/':
        output_dir += '/'
    if not os.path.exists(output_dir + 'output/'):
        os.makedirs(output_dir + 'output/')
    out = output_dir + 'output/'

    if arguments.inputPath:
        input_file = arguments.inputPath
    else:
        input_file = out
    if arguments.min:
        min_frag = int(arguments.min)
    else:
        min_frag = 20  # If no minimum fragment length is given, the default is 20
    if arguments.max:
        max_frag = int(arguments.max)
    else:
        max_frag = 35

    if not os.path.exists(out+'Pickle_dicts/'):
        os.makedirs(out+'Pickle_dicts/')
    if arguments.number:
        threshold_avg_reads = float(arguments.number)
    else:
        threshold_avg_reads = 1
    if arguments.threshold:
        threshold_cutoff = int(arguments.threshold)
    else:
        threshold_cutoff = 70
    if arguments.second_threshold:
        second_threshold = int(arguments.second_threshold)
    else:
        second_threshold = 5
    input_file2 = arguments.inputPath2
    if arguments.three_prime == 'Yes':
        three_prime_end = True
        print 'Doing analysis for 3\' aligned reads'
    else:
        three_prime_end = False
    # Filter genes which have greater than threshold (default=1) reads per codon on average
    if not os.path.isfile(out+"Pickle_dicts/filtered_dict.p"):
        filtered_genes, dataset_gene_len = select_high_cov_genes(input_file, input_file2, min_frag, max_frag, out, threshold_avg_reads, three_prime=three_prime_end)
        print 'Parsed the CDS file'
        Pickle.dump(filtered_genes, open(out+"Pickle_dicts/filtered_dict.p", "wb"))
        Pickle.dump(dataset_gene_len, open(out+"Pickle_dicts/gene_len.p", "wb"))
    else:
        filtered_genes = Pickle.load(open(out+"Pickle_dicts/filtered_dict.p", "rb"))
        print ' Need not select high coverage genes. Loaded the Pickle dict for filtered genes'
        dataset_gene_len = Pickle.load(open(out+"Pickle_dicts/gene_len.p", "rb"))
    offset_dict = asite_algorithm_improved_second_offset_correction(filtered_genes, dataset_gene_len, min_frag, max_frag, out, off_correction_threshold=second_threshold,
                                                                    bootstrap=False, offset_threshold=threshold_cutoff, three_prime=three_prime_end)
