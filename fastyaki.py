import csv
import os
#import pandas as pd

#from matplotlib import pyplot as plt
#from matplotlib import cm
from ont_fast5_api.fast5_interface import get_fast5_file


def print_all_raw_data():
    fast5_filepath = 'test/mapped_reads.hdf5'  # This can be a single- or multi-read file
    with get_fast5_file(fast5_filepath, mode='r') as f5:
        for read_id in f5.get_read_ids():
            read = f5.get_read(read_id)
            raw_data = read.get_raw_data()
            print(read_id, raw_data)


def get_unique_bacteria_parts():
    umi_ref_link = 'test/umi_ref_link.csv'
    bacteria = {}
    with open(umi_ref_link) as file:
        reader = csv.reader(file)
        for row in reader:
            bac = row[1]
            if not bac in bacteria:
                bacteria[bac] = 0
            bacteria[bac] += 1
    return bacteria


def get_unique_bacteria():
    bacteria_parts = get_unique_bacteria_parts()
    bacteria = {}
    for part in bacteria_parts.keys():
        bac = part.split('_')[0]
        if not bac in bacteria:
            bacteria[bac] = 0
        bacteria[bac] += bacteria_parts[part]
    return bacteria


def print_bacteria_part_umibins(bacteria_part):
    umi_ref_link = 'test/umi_ref_link.csv'
    umi_bins = []
    with open(umi_ref_link) as file:
        reader = csv.reader(file)
        for row in reader:
            if row[1] == bacteria_part:
                umi_bins.append(row[0])
    print('.fast5 '.join(umi_bins))


def print_bacteria_umibins(bacteria):
    all_parts = get_unique_bacteria_parts()
    bac_parts = [part for part, value in all_parts.items() if bacteria in part]
    for bac_part in bac_parts:
        print_bacteria_part_umibins(bac_part)


def n_bases(reference_path):
    res = 0
    with open(reference_path, 'r') as ref_file:
        for line in ref_file:
            if line.startswith('<'):
                continue
            res += len(line.strip())
    return res


def get_bacteria_bins(umi_ref_file, bacterias, only_part=False):
    res = {}
    with open(umi_ref_file) as file:
        reader = csv.reader(file)
        for row in reader:

            bin, bac_part = row
            if not only_part:
                bac = bac_part.split('_')[0]
            else:
                bac = bac_part

            if bac in bacterias:
                if res.get(bac_part) is None:
                    res[bac_part] = []
                res[bac_part].append(bin)

    return res


def reads_from_bins(bins, bins_dir):
    reads = []
    for bin in bins:
        with get_fast5_file(bins_dir + bin + '.fast5', mode='r') as f5:
            reads.extend(f5.get_read_ids())
    return reads


def get_bacteria_reads(umi_ref_file, bins_dir, bacterias, only_part=False):
    reads = {}
    bacteria_bins = get_bacteria_bins(umi_ref_file, bacterias, only_part)
    for bac, bins in bacteria_bins.items():
        reads[bac] = reads_from_bins(bins, bins_dir)
    return reads


def split_bac_reads_evenly(bac_part_reads, distribution=0.8, max_reads=None):
    part1 = []
    part2 = []

    if max_reads:
        max_reads_per_part = int(max_reads / len(bac_part_reads))
        for bac_part in bac_part_reads:
            bac_part_reads[bac_part] = bac_part_reads[bac_part][:max_reads_per_part]

    for reads in bac_part_reads.values():
        split_at = int(len(reads) * distribution)
        part1.extend(reads[:split_at])
        part2.extend(reads[split_at + 1:])
    return part1, part2


def reads_to_tsv(reads, filename):
    with open(filename, 'w+') as file:
        file.write('read_id\n')
        file.write('\n'.join(reads))


def make_read_id_to_fasq_file(fastq_dir):
    output = ''
    files = [f for f in os.listdir(fastq_dir) if f.endswith('.fastq')]
    for file in files:
        with open(os.path.join(fastq_dir, file)) as fastq:
            headers = [l for l in fastq if l.startswith('@')]
            for header in headers:
                header = header[1:]  # remove @
                header = header[:header.find(' ')]  # remove everything after readid
                output += f'{header}\t{file}\n'
    with open(os.path.join(fastq_dir, 'read_id_to_fastq.tsv'), 'w+') as file:
        file.write(output)


def avg_read_identity_from_tsv(tsv_path):
    sum = 0
    count = 0
    with open(tsv_path) as file:
        next(file)  # skip header
        for line in file:
            identity = float(line.split()[2])
            if identity == 0.0:
                continue
            sum += identity
            count += 1
    return sum / count


#cmaps = OrderedDict()
#print_all_raw_data()
#bacteria = get_unique_bacteria()
#print(bacteria)
#get_bacteria_bins('test/umi_ref_link.csv', 'Salmonella')
#print_bacteria_umibins('Salmonella')
reads = get_bacteria_reads('test/umi_ref_link.csv', '/home/mac/data/100GBraw/single/', 'Salmonella_5', only_part=True)
#train, validation = split_bac_reads_evenly(reads, distribution=1)
#reads_to_tsv(reads, 'test/Salmonella.tsv')
#reads_to_tsv(validation, 'test/Escherichia_validation.tsv')
print(reads)
#print_bacteria_umibins('Escherichia_1')
#print(n_bases('test/zymo-ref-uniq_2019-03-15.fa'))
#print(avg_read_identity_from_tsv('/home/mac/winlenbasecalls/winlen9_identity.tsv'))
#print(avg_read_identity_from_tsv('/home/mac/winlenbasecalls/winlen19_identity.tsv'))
#print(avg_read_identity_from_tsv('/home/mac/winlenbasecalls/winlen28_identity.tsv'))
#print(avg_read_identity_from_tsv('/home/mac/winlenbasecalls/guppy_standard_identity.tsv'))
#make_read_id_to_fasq_file('/home/mac/winlenbasecalls/winlen9basecalls')