import vcf
import pysam

# read VCF file
def read_vcf_file(vcf_path):
    vcf_reader = vcf.Reader(open(vcf_path, "r"))
    return list(vcf_reader)


# get all mutant sites
def get_variants_positions(vcf_file):
    positions = []
    vcf_records = read_vcf_file(vcf_file)
    for record in vcf_records:
        chrom_pos = f"{record.CHROM}:{record.POS}"
        positions.append(chrom_pos)
    return positions


# write AD.tsv header
def write_ad_header(ad_file, positions):
    ad_file.write("cell_id,{}\n".format(",".join(positions)))


# write AD.tsv data
def write_ad_data(ad_file, vcf_file, positions, counter):
    vcf_records = read_vcf_file(vcf_file)
    # init row
    row = ["0"] * len(positions)
    #  row
    for record in vcf_records:
        chrom_pos = f"{record.CHROM}:{record.POS}"
        index = positions.index(chrom_pos)
        # get variant depth
        sample_record = record.samples[counter]
        sample_name = sample_record.sample
        if sample_record["AD"]:
            depth = sample_record["AD"][1]
            row[index] = str(depth)
    # write row
    print(sample_name + "   AD       write!")
    ad_file.write("{},{}\n".format(sample_name, ",".join(row)))


# write DP.tsv header
def write_dp_header(dp_file, positions):
    dp_file.write("cell_id,{}\n".format(",".join(positions)))


# write DP.tsv data
def write_dp_data(dp_file, vcf_file, positions, counter):
    vcf_records = read_vcf_file(vcf_file)
    # init row
    row = ["0"] * len(positions)
    # write row
    for record in vcf_records:
        chrom_pos = f"{record.CHROM}:{record.POS}"
        index = positions.index(chrom_pos)
        # get variant depth
        sample_record = record.samples[counter]
        sample_name = sample_record.sample
        if sample_record["DP"]:
            depth = sample_record["DP"]
            row[index] = str(depth)
    # write row
    print(sample_name + "  DP   write!")
    dp_file.write("{},{}\n".format(sample_name, ",".join(row)))


# main
def main(vcf_file, ad_file, dp_file):
    # all variant position
    positions = get_variants_positions(vcf_file)
    # write AD.tsv and DP.tsv header
    write_ad_header(ad_file, positions)
    write_dp_header(dp_file, positions)

    vcf_records = read_vcf_file(vcf_file)
    for i in range(0, len(vcf_records[0].samples)):
        write_ad_data(ad_file, vcf_file, positions, i)
        write_dp_data(dp_file, vcf_file, positions, i)


# call main function
vcf_file = "monovar.vcf"
with open("AD.csv", "w") as ad_file, open("DP.csv", "w") as dp_file:
    main(vcf_file, ad_file, dp_file)
