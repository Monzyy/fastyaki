import argparse
import h5py
import helpers
import os
import random
import shutil


def bac2read(args):
    bac_to_read = {}
    bac_bins = helpers.get_bacteria_bins(args.umi_ref_link)

    for bac, bins in bac_bins.items():
        bac_to_read[bac] = helpers.reads_from_bins(bins, args.bins_dir)

    if args.hdf5source:
        hdf5 = h5py.File(args.hdf5source, 'r')
        h5groups = [val.name.split('/')[-1] for val in hdf5['Reads'].values()]

        for bac, read_ids in bac_to_read.items():
            intersect = set(h5groups).intersection(read_ids)
            bac_to_read[bac] = list(intersect)

    with open(args.output, 'w+') as outfile:
        for bac, read_ids in bac_to_read.items():
            outfile.write(f'{bac},{",".join(read_ids)}\n')


def split_bac8(args):
    bac2read_dict = helpers.read2bac_file_to_dict(args.bac2read)
    hdf5 = h5py.File(args.source, 'r')
    source = hdf5['Reads']

    test_bacs = args.test_bacs.split(',')

    n_not_found = 0
    os.mkdir(args.dest)
    train_set = helpers.h5_create_copy_without_reads(os.path.join(args.dest, args.train_output), hdf5)
    dest_train_read_group = train_set['Reads']
    test_set = helpers.h5_create_copy_without_reads(os.path.join(args.dest, args.test_output), hdf5)
    dest_test_read_group = test_set['Reads']

    # Make training set
    for bac_part, reads_ids in bac2read_dict.items():
        if any((bac_part.startswith(bac) for bac in test_bacs)):
            continue
        subset = int(len(reads_ids) * args.percentage)
        train_groups = random.sample(reads_ids, subset)

        n_not_found += helpers.copy_h5_groups(source, dest_train_read_group, train_groups)

    # Make test set
    os.mkdir(os.path.join(args.dest, 'fast5_testdata'))
    bac_bins = helpers.get_bacteria_bins(args.umi_ref_link)
    for bac_part, reads_ids in bac2read_dict.items():
        if any((not bac_part.startswith(bac) for bac in test_bacs)):
            continue
        subset = int(len(reads_ids) * args.percentage)
        test_groups = random.sample(reads_ids, subset)

        n_not_found += helpers.copy_h5_groups(source, dest_test_read_group, test_groups)

        # Copy fast5 files
        for bin in bac_bins[bac_part]:
            shutil.copyfile(os.path.join(args.fast5_dir, f'{bin}.fast5'),
                            os.path.join(args.dest, f'fast5_testdata/{bin}.fast5'))


    hdf5.close()
    train_set.close()
    test_set.close()

    if n_not_found:
        print(f'{n_not_found} read_ids were not found in source file.'
              f' Reads may have been filtered by prepare_mapped_reads.py')
    print('Done')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Split a dataset.')
    subparsers = parser.add_subparsers(help='Sub-command help')

    # Subparsers
    bac8parser = subparsers.add_parser(
        'bac8', help='Process the 8 PCR sequences bacteria from Department of Chemistry and Bioscience')

    # Bac8
    bac8subparser = bac8parser.add_subparsers(help='Bac8 sub-command help')

    # Bac8 read2bac
    bac8readtobacparser = bac8subparser.add_parser('read2bac', help='Create a bac8 read2bac file')
    bac8readtobacparser.set_defaults(func=bac2read)
    bac8readtobacparser.add_argument('umi_ref_link', help='Path to umi_ref_link.csv file')
    bac8readtobacparser.add_argument('bins_dir', help='Path to umi bin directory')
    bac8readtobacparser.add_argument('-o', '--output', default='bac8_to_read.csv', help='Output file')
    bac8readtobacparser.add_argument('-5', '--hdf5source',
                                     help='Prune reads that does not exist in hdf5source file.')

    # Bac8 split data
    bac8splitdata = bac8subparser.add_parser('split', help='Create a subset of the bac8 dataset.')
    bac8splitdata.add_argument('bac2read', help='Bacteria to read file.')
    bac8splitdata.add_argument('source', help='Source HDF5 file containing mapped reads.')
    bac8splitdata.add_argument('test_bacs', type=str, help='Comma separated list of bacterias to use in the test '
                                                           'set. Use the rest in the train set.')
    bac8splitdata.add_argument('fast5_dir', help='Path to fast5 data directory')
    bac8splitdata.add_argument('umi_ref_link', help='Path to umi_ref_link.csv file')
    bac8splitdata.add_argument('-p', '--percentage', type=float, default=1.0,
                               help='Use a percentage of dataset, distributed evenly over bacteria')
    bac8splitdata.add_argument('--dest', default='output', help='Directory to save new data sets')
    bac8splitdata.add_argument('--train_output', default='train.hdf5',
                               help='File name for HDF5 training set.')
    bac8splitdata.add_argument('--test_output', default='test.hdf5',
                               help='File name for HDF5 test set.')
    bac8splitdata.set_defaults(func=split_bac8)

    args = parser.parse_args()
    args.func(args)
