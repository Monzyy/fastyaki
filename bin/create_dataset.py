import argparse
import h5py
import helpers


def bac2read(args):
    bac_to_read = {}
    bac_bins = helpers.get_bacteria_bins(args.umi_ref_link)

    for bac, bins in bac_bins.items():
        bac_to_read[bac] = helpers.reads_from_bins(bins, args.bins_dir)

    if args.hdf5source:
        hdf5 = h5py.File(args.hdf5source)
        h = hdf5['Reads']
        for bac, read_ids in bac_to_read.items():
            for idx in read_ids:
                if h.get(idx) is None:
                    bac_to_read[bac].remove(idx)

    with open(args.output, 'w+') as outfile:
        for bac, read_ids in bac_to_read.items():
            outfile.write(f'{bac},{",".join(read_ids)}\n')


def split_bac8(args):
    bac2read_dict = helpers.read2bac_file_to_dict(args.bac2read)
    hdf5 = h5py.File(args.source)
    source = hdf5['Reads']

    output = helpers.h5_create_copy_without_reads(args.output, hdf5)
    dest_read_group = output['Reads']

    n_not_found = 0
    if args.test:
        test_set = helpers.h5_create_copy_without_reads(args.test_output, hdf5)
        dest_test_read_group = test_set['Reads']

        for reads_ids in bac2read_dict.values():
            subset = int(len(reads_ids) * args.percentage)
            split_index = int(subset * args.test)
            test_groups = reads_ids[:split_index]
            train_groups = reads_ids[split_index:subset]

            n_not_found += helpers.copy_h5_groups(source, dest_read_group, train_groups)
            n_not_found += helpers.copy_h5_groups(source, dest_test_read_group, test_groups)

        test_set.close()

    else:
        for reads_ids in bac2read_dict.values():
            subset = int(len(reads_ids) * args.percentage)
            groups = reads_ids[:subset]
            n_not_found += helpers.copy_h5_groups(source, dest_read_group, groups)

    hdf5.close()
    output.close()

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
    bac8splitdata.add_argument('-o', '--output', default='output.hdf5',
                               help='Output HDF5 file containing a subset of mapped reads.')
    bac8splitdata.add_argument('-p', '--percentage', type=float, default=1.0,
                               help='Use a percentage of dataset, distributed evenly over bacteria')
    bac8splitdata.add_argument('-t', '--test', type=float,
                               help='Split the dataset into a training and test set, '
                                    'and choose size of test set in percentage')
    bac8splitdata.add_argument('--test_output', default='test.hdf5',
                               help='Output file for test set, if option -t is used. Default is test_out.hdf5.')
    bac8splitdata.set_defaults(func=split_bac8)

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()
