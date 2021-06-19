import argparse
import sys
import operator
import distance

from datetime import datetime
from pysam import VariantFile, FastaFile
from collections import defaultdict

"""
This script generates variants based on the depth and other small variants we observe in the telomeric regions of the genome. This was the script we ran to polish
T2T v1.0:

telomere_depths.bed : bed file with coverage generated with samtools depth
chm13.draft_v1.0.telomere.bed : telomere annotated regions in v1.0
chm13.draft_v1.0.fasta : v1.0 assembly
CHM13_v1.0_telomere_edits_kmer_sensitive.vcf: Output file
CHM13v1_ONT_PEPPER_HP_CANDIDATES.vcf.gz: Input file


python3 generate_telomere_edits.py \
--telomere_depth_bed telomere_depths.bed \
--telomere_annotation chm13.draft_v1.0.telomere.bed \
--fasta chm13.draft_v1.0.fasta \
--output_vcf CHM13_v1.0_telomere_edits_kmer_sensitive.vcf \
--small_variant_vcf CHM13v1_ONT_PEPPER_HP_CANDIDATES.vcf.gz \
--min_depth 5 \
--min_gq 2 \
--min_vaf 0.5
"""


def is_restoring_canonical_kmer(contig, start_pos, end_pos, reference_allele, alternate_allele, canonical_kmer, assembly_fasta_file, contig_length):
    left_start = max(1, start_pos - 10)
    right_end = min(contig_length, end_pos + 10)

    contig_sequence = assembly_fasta_file.fetch(reference=contig, start=left_start, end=right_end)

    canonical_repeat = canonical_kmer * int(len(contig_sequence) / len(canonical_kmer))
    old_distance_with_canonical = distance.levenshtein(canonical_repeat, contig_sequence)

    replaced_sequence = assembly_fasta_file.fetch(reference=contig, start=left_start, end=start_pos) + alternate_allele + assembly_fasta_file.fetch(reference=contig, start=end_pos, end=right_end)
    new_distance_with_canonical = distance.levenshtein(canonical_repeat, replaced_sequence)

    # if contig == "chr19" and start_pos < 200:
    #     print(left_start, start_pos)
    #     print(assembly_fasta_file.fetch(reference=contig, start=left_start, end=start_pos))
    #     print(alternate_allele)
    #     print(assembly_fasta_file.fetch(reference=contig, start=end_pos, end=right_end))
    #     print(left_start)
    #     print(right_end)
    #     print(contig_sequence)
    #     print(replaced_sequence)
    #     print(old_distance_with_canonical)
    #     print(new_distance_with_canonical)

    if new_distance_with_canonical < old_distance_with_canonical:
        return 1
    elif new_distance_with_canonical > old_distance_with_canonical:
        return -1
    else:
        return 0



def telomere_pruning(telomere_depth_bed, telomere_annotation, fasta, small_variant_vcf, output_vcf, min_depth, min_gq, min_vaf):
    """
    Find regions to delete in the telomere
    :param telomere_depth_bed:
    :param telomere_annotation:
    :param fasta:
    :param output_vcf:
    :return:
    """
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: READING DEPTH BED FILE " + "\n")
    sys.stderr.flush()

    assembly_fasta_file = FastaFile(fasta)

    # directionary to keep track of depth at each position of the telomere
    position_wise_depth = defaultdict()

    # outputs regions that are going to be edited
    telomere_edit_regions = open("CHM13_v1_telomere_edit_regions.bed", "w")

    # read the depth file
    depth_bed_file = open(telomere_depth_bed, "r")

    small_variant_vcf = VariantFile(small_variant_vcf)
    output_vcf_file = VariantFile(output_vcf, 'w', header=small_variant_vcf.header)

    # populate the position dictionary
    for bed_record in depth_bed_file:
        contig, position, depth = bed_record.rstrip().split("\t")
        # the bedfile has an offset of 1
        position_wise_depth[(contig, int(position) - 1)] = int(depth)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: DEPTH BED LOADED. " + "\n")
    sys.stderr.flush()

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: READING ANNOTATION BED. " + "\n")
    sys.stderr.flush()

    regions_of_deletion = defaultdict(lambda: list)
    telomere_regions = defaultdict(lambda: list)
    contig_length_dict = defaultdict()
    telomere_annotation_bed_file = open(telomere_annotation, "r")

    all_vcf_records = []

    for bed_record in telomere_annotation_bed_file:
        contig, start_pos, end_pos, contig_length = bed_record.rstrip().split("\t")
        contig_length_dict[contig] = contig_length
        if contig not in regions_of_deletion.keys():
            regions_of_deletion[contig] = []

        if contig not in telomere_regions.keys():
            telomere_regions[contig] = []

        telomere_regions[contig].append((int(start_pos), int(end_pos)))

        start_pos = int(start_pos)
        end_pos = int(end_pos)
        contig_length = int(contig_length)

        if start_pos == 0:
            # this is the left side of the telomere, so scan left to right
            if (contig, 0) in position_wise_depth.keys() and position_wise_depth[(contig, 0)] >= min_depth:
                # it has full coverage, so simply do nothing.
                continue
            # otherwise scan to the point we hit min_depth
            record_start_pos = 0
            current_position = 1

            while True:
                current_depth = 0
                if (contig, current_position) in position_wise_depth.keys():
                    current_depth = position_wise_depth[(contig, current_position)]

                if current_depth >= min_depth or current_position == end_pos:
                    break
                current_position += 1
            record_end_position = current_position
            length_of_record = record_end_position - record_start_pos + 1
            # pad the reference allele by one place
            reference_allele = assembly_fasta_file.fetch(reference=contig, start=record_start_pos, end=record_end_position + 1)
            # alternate allele is the last base of the
            alternate_allele = reference_allele[-1]

            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PRUNING: " + contig + " " + str(record_start_pos) + " " + str(record_end_position) + " " + str(length_of_record) + "\n")
            sys.stderr.flush()

            regions_of_deletion[contig].append((record_start_pos, record_end_position))
            telomere_edit_regions.write(contig + "\t" + str(record_start_pos) + "\t" + str(record_end_position) + "\n")

            # write this deletion to the VCF file
            alleles = [reference_allele, alternate_allele]
            vcf_record = output_vcf_file.new_record(contig=contig, start=record_start_pos,
                                                    stop=record_end_position + 1, id='.', qual=60,
                                                    filter='PASS', alleles=alleles, GT=[1, 1],
                                                    GQ=60, VAF=[1.0])
            all_vcf_records.append((vcf_record.contig, vcf_record.start, vcf_record.stop, vcf_record))
        elif end_pos == contig_length:
            # this is the right side of the telomere, so scan right to left
            if (contig, end_pos) in position_wise_depth.keys() and position_wise_depth[(contig, end_pos)] >= min_depth:
                # it has full coverage, so simply do nothing.
                continue

            record_end_position = end_pos
            current_position = end_pos - 1
            while True:
                current_depth = 0
                if (contig, current_position) in position_wise_depth.keys():
                    current_depth = position_wise_depth[(contig, current_position)]

                if current_depth >= min_depth or current_position == start_pos:
                    break
                current_position -= 1
            record_start_pos = current_position
            length_of_record = record_end_position - record_start_pos + 1

            # pad the reference allele by one place
            reference_allele = assembly_fasta_file.fetch(reference=contig, start=record_start_pos - 1, end=record_end_position)
            # alternate allele is the last base of the
            alternate_allele = reference_allele[0]
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PRUNING: " + contig + " " + str(record_start_pos) + " " + str(record_end_position) + " " + str(length_of_record) + "\n")
            sys.stderr.flush()

            regions_of_deletion[contig].append((record_start_pos, record_end_position))
            telomere_edit_regions.write(contig + "\t" + str(record_start_pos) + "\t" + str(record_end_position) + "\n")

            # write this deletion to the VCF file
            alleles = [reference_allele, alternate_allele]
            vcf_record = output_vcf_file.new_record(contig=contig, start=record_start_pos - 1,
                                                    stop=record_end_position, id='.', qual=60,
                                                    filter='PASS', alleles=alleles, GT=[1, 1],
                                                    GQ=60, VAF=[1.0])
            all_vcf_records.append((vcf_record.contig, vcf_record.start, vcf_record.stop, vcf_record))

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: READING SMALL VARIANT VCF. " + "\n")
    sys.stderr.flush()

    # true_positive_positions = defaultdict(list)
    #
    # # filter the file
    for rec in small_variant_vcf.fetch():
        contig_deletion_regions = list(regions_of_deletion[rec.chrom])
        record_overlaps = False
        for region_start, region_end in contig_deletion_regions:
            if region_start <= rec.pos <= region_end:
                record_overlaps = True
                break

        # this small variant overalps with a region we are going to delete.
        if record_overlaps:
            continue

        telomere_region = list(telomere_regions[rec.chrom])
        record_overlaps = False
        for region_start, region_end in telomere_region:
            if region_start <= rec.pos <= region_end:
                record_overlaps = True
                break

        # small variant is outside telomere region
        if record_overlaps is False:
            continue

        sample_vafs = []
        for sample in rec.samples:
            sample_vafs = rec.samples[sample]['VAF']

        selected_alleles = [rec.alleles[0]]
        selected_allele_vaf = []

        for i in range(0, len(rec.alts)):



            if rec.pos < 10000:
                restoring_canonical = is_restoring_canonical_kmer(rec.contig, rec.start, rec.stop, rec.alleles[0], rec.alts[i], "CCCTAA", assembly_fasta_file, int(contig_length_dict[rec.contig]))
            else:
                restoring_canonical = is_restoring_canonical_kmer(rec.contig, rec.start, rec.stop, rec.alleles[0], rec.alts[i], "GGGTTA", assembly_fasta_file, int(contig_length_dict[rec.contig]))

            # if rec.contig == "chr19" and rec.pos < 200:
            #     print(rec, end='')
            #     print(restoring_canonical)

            # restoring canonical: 0 is no change, 1 is positive change, -1 means it's moving away from canonical
            if restoring_canonical == 1:
                found_positive_change = True
                selected_alleles.append(rec.alts[i])
                selected_allele_vaf.append(sample_vafs[i])
            elif restoring_canonical == 0:
                # meaning this allele has no affect on canonical k-mer restoration, so we simply fall back to set thresholds.
                if sample_vafs[i] >= min_vaf and rec.qual >= min_gq:
                    selected_alleles.append(rec.alts[i])
                    selected_allele_vaf.append(sample_vafs[i])

        # no allele passed the thresholds or is restoring canonical k-mer
        if len(selected_alleles) == 1:
            continue

        vcf_record = output_vcf_file.new_record(contig=rec.contig, start=rec.start,
                                                stop=rec.stop, id='.', qual=rec.qual,
                                                filter='PASS', alleles=selected_alleles, GT=[1, 1],
                                                GQ=rec.qual, VAF=selected_allele_vaf)

        all_vcf_records.append((vcf_record.contig, vcf_record.start, vcf_record.stop, vcf_record))

        telomere_edit_regions.write(rec.contig + "\t" + str(rec.start) + "\t" + str(rec.stop) + "\n")

    all_vcf_records = sorted(all_vcf_records, key=lambda x: (x[0], x[1], x[2]))
    for contig, start, stop, record in all_vcf_records:
        output_vcf_file.write(record)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--telomere_depth_bed",
        type=str,
        required=True,
        help="Path to a bed file describing the depth at each location in the telomere position."
    )
    parser.add_argument(
        "--telomere_annotation",
        "-v2",
        type=str,
        required=True,
        help="Path to telomere annotation bed file."
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Path to the assembly file."
    )
    parser.add_argument(
        "--small_variant_vcf",
        type=str,
        required=True,
        help="Path to VCF with small variants to apply."
    )
    parser.add_argument(
        "--min_gq",
        type=int,
        required=True,
        help="Min genotype quality for small variant."
    )
    parser.add_argument(
        "--min_vaf",
        type=float,
        required=True,
        help="Min variant allele frequency for small variant."
    )
    parser.add_argument(
        "--min_depth",
        type=int,
        required=True,
        help="Min depth allowed in telomere."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    telomere_pruning(FLAGS.telomere_depth_bed, FLAGS.telomere_annotation, FLAGS.fasta, FLAGS.small_variant_vcf, FLAGS.output_vcf, FLAGS.min_depth, FLAGS.min_gq, FLAGS.min_vaf)

