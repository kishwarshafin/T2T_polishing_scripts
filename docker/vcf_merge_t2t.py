import argparse
import sys
import operator

from datetime import datetime
from pysam import VariantFile
from collections import defaultdict

"""
This script takes two VCF files and output from hap.py to output variants that are union of the two sets.
Way to run it:
1. Run hap.py with the two VCFs.

time docker run -it -v /data:/data \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
/data/Variant_calls_1.vcf.gz \
/data/Variant_calls_2.vcf.gz \
-r /data/CHM13_T2T.fa \
-o /data/output/variant_1_vs_variant1_chr20 \
--pass-only \
--engine=vcfeval \
--threads=32

This will generate variant_1_vs_variant_2_chr20.vcf.gz in the output, this will have TP annotation for variants which will indicate which variants are present in both sets.

Then run this script:
python3 vcf_merge_t2t.py \
-v1 /data/Variant_calls_1.vcf.gz  \
-v2 /data/Variant_calls_2.vcf.gz \
-hv /data/output/variant_1_vs_variant1_chr20.vcf.gz \
-o /data/VARIANT_12_MERGED.vcf
"""


def vcf_merge_vcfs(in_vcf1, in_vcf2, happy_vcf, output_vcf):
    """
    Merge two vcf files
    :param in_vcf1: Input VCF file 1
    :param in_vcf2: Input VCF file 2
    :param happy_vcf: Hap.py input file
    :param output_vcf: Output VCF file
    :return:
    """
    happy_vcf_file = VariantFile(happy_vcf)

    # counter to keep track to of true positive cases
    true_positive_positions = defaultdict(list)

    # filter the file
    for rec in happy_vcf_file.fetch():
        for sample in rec.samples:
            sample_bd = rec.samples[sample]['BD']
            if sample_bd == 'TP':
                # record a true positive case
                true_positive_positions[rec.contig].append(rec.pos)

    # read the two inpt files
    vcf1_vcf_file = VariantFile(in_vcf1)
    vcf2_vcf_file = VariantFile(in_vcf2)

    # for VCF1 we add all of the records.
    merged_records = []
    position_dict = set()
    for rec in vcf1_vcf_file.fetch():
        position_dict.add((rec.contig, rec.pos))
        merged_records.append((rec.contig, rec.pos, rec))

    # for VCF2 we add records that are not true positives
    for rec in vcf2_vcf_file.fetch():
        if rec.pos not in true_positive_positions[rec.contig]:
            if (rec.contig, rec.pos) not in position_dict:
                merged_records.append((rec.contig, rec.pos, rec))

    # sort the records
    merged_records.sort(key=operator.itemgetter(0, 1))

    # output the file
    vcf_out = VariantFile(output_vcf, 'w', header=vcf1_vcf_file.header)

    # write the VCF
    for cotig, pos, rec in merged_records:
        vcf_out.write(rec)

    # process completed
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED " + "\n")
    sys.stderr.flush()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_vcf1",
        "-v1",
        type=str,
        required=True,
        help="Path to input vcf1."
    )
    parser.add_argument(
        "--input_vcf2",
        "-v2",
        type=str,
        required=True,
        help="Path to input vcf2."
    )
    parser.add_argument(
        "--happy_vcf",
        "-hv",
        type=str,
        required=True,
        help="Path to happy vcf."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_merge_vcfs(FLAGS.input_vcf1, FLAGS.input_vcf2, FLAGS.happy_vcf, FLAGS.output_vcf)

