import csv
import os

from ont_fast5_api.fast5_interface import get_fast5_file


def get_bacteria_bins(umi_ref_file, in_parts=True,):
    res = {}
    with open(umi_ref_file) as file:
        reader = csv.reader(file)
        for row in reader:

            bin, bac_part = row
            if not in_parts:
                bac = bac_part.split('_')[0]
            else:
                bac = bac_part

            if res.get(bac_part) is None:
                res[bac_part] = []
            res[bac_part].append(bin)

    return res


def reads_from_bins(bins, bins_dir):
    if isinstance(bins, str):
        bins = [bins]
    reads = []
    for bin in bins:
        with get_fast5_file(os.path.join(bins_dir, bin + '.fast5'), mode='r') as f5:
            reads.extend(f5.get_read_ids())
    return reads


def read2bac_file_to_dict(read2bac_file):
    read2bac = {}
    with open(read2bac_file) as file:
        reader = csv.reader(file)
        for row in reader:
            read2bac[row[0]] = row[1:]
    return read2bac


def copy_h5_groups(source, dest, group_names):
    n_not_found = 0
    for group in group_names:
        try:
            source.copy(group, dest)
        except KeyError:
            n_not_found += 1
    return n_not_found