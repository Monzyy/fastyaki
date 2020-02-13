import editdistance
import h5py
import pandas as pd
from multiprocessing import cpu_count, Pool
import numpy as np
from difflib import SequenceMatcher


def h5_to_pickle(in_path, out_path):
    data = h5py.File(in_path, 'r')['Reads']
    data_slice_ids = [id for i, id in enumerate(data)]
    data_slice = [data[id]['Reference'] for id in data_slice_ids]
    reads = np.full([len(data_slice), len(max(data_slice, key=lambda x: len(x)))], -1)
    for i, j in enumerate(data_slice):
        reads[i][0:len(j)] = j

    df = pd.DataFrame(reads, data_slice_ids)
    df.to_pickle(out_path)


#h5_to_pickle('test/mapped_reads2.hdf5', 'test/escheria_1_reads.pkl')

#df = pd.read_pickle('test/escheria_1_reads.pkl')
#print(df)

with open('test/flipped_escheria.txt', 'r') as flipped_file:
    # Load read ids of flipped reads as set
    flipped_reads = set(flipped_file.read().replace(',', '').split())

# Dict that specifies which how to tranlate a base to its complementary base.
# E.g. base 3 should be converted to base 0
_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'X': 'X', 'N': 'N',
               'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'x': 'x', 'n': 'n',
               '-': '-'}


def complement(seq, compdict=_COMPLEMENT):
    return ''.join(compdict[b] for b in seq)


def reverse_complement(seq, compdict=_COMPLEMENT):
    return complement(seq, compdict)[::-1]


joined_bases = []


def int_list_to_str(int_list):
    return ('{} '*len(int_list)).format(*int_list)


def align_strands(ref):
    if ref.name in flipped_reads:
        return reverse_complement(str(ref.iloc[0]))
    return ref


def match_substring(m, n):
    seq_match = SequenceMatcher(None, list(m), list(n))
    return seq_match.find_longest_match(0, len(m), 0, len(n))


def load_read_refs_to_df():
    references = []
    read_ids = []
    with open('test/read_references.fasta') as f:
        for line in f:
            if line.startswith('>'):
                read_ids.append(line.replace('>', '').strip())
            else:
                references.append(line.strip())
    return pd.DataFrame(references, read_ids, dtype=str)


def get_reference():
    reference = ''
    with open('test/reference.fa') as f:
        for line in f:
            if line.startswith('>'):
                continue
            reference += line.strip()
    return reference


def lcss_on_dataframe(df, first=None):
    lcss = ''
    rest = df
    if first is None:
        first = str(df.iloc[0])
        rest = df.iloc[1:]
    for i in range(len(first)):
        for j in range(len(first) - i + 1):
            if j > len(lcss):
                cur_substring = first[i:i + j]
                is_css = rest.str.contains(cur_substring)
                if is_css.all():
                    lcss = cur_substring
                else:
                    break
    return lcss


def parallelize_dataframe(df, func, n_cores=1):
    df_split = np.array_split(df, n_cores, axis=0)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


def parallelize_dataframe_with_second_param(params, func, n_cores=1):
    df, second_param = params
    df_split = np.array_split(df, n_cores, axis=0)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, (second_param, df_split)))
    pool.close()
    pool.join()
    return df


def lccs_compare(df):
    return df.apply(lambda x: lcss_on_dataframe(x))


def lcss_individually_to_pkl(filter_ids=None):
    df = load_read_refs_to_df()
    #df = df.apply(lambda x: align_strands(x), axis=1)
    if filter_ids is not None:
        df = df.loc[filter_ids]
    #df['strand'] = [True if x in flipped_reads else False for x in df.index.values]

    print('lcss_individually')
    lcss_individually = parallelize_dataframe(df, lccs_compare)
    lcss_individually = pd.DataFrame(lcss_individually.values, index=lcss_individually.index, columns=['lcss', ])
    lcss_individually['length'] = lcss_individually['lcss'].str.len()
    lcss_individually = lcss_individually.sort_values(by=['length'])
    pd.to_pickle(lcss_individually, 'test/lcss_individually1.pkl')
    print(lcss_individually)


def load_lcss_individually_from_pkl():
    df = pd.read_pickle('test/lcss_individually.pkl')
    df = df[df['length'] == 10]
    return df


def df_names_to_file(df):
    with open('test/only_10_lcss1.txt', 'w') as f:
        f.write(' '.join(list(df.index)))

#df = load_lcss_individually_from_pkl()
#df_names_to_file(df)


with open('test/only_10_lcss1.txt', 'r') as f:
    only_10_lcss = f.read().split()


#lcss_individually_to_pkl(only_10_lcss)

def find_lcss(x, ys):
    lcss = ''
    if isinstance(ys, str):
        ys = [ys]
    if isinstance(ys, list):
        ys = pd.Series(ys)
    if isinstance(x, pd.Series):
        x = str(x.iloc[0])
    for i in range(len(x)):
        for j in range(len(x) - i + 1):
            if j > len(lcss):
                cur_substring = x[i:i + j]
                is_css = ys.str.contains(cur_substring)
                if is_css.all():
                    lcss = cur_substring
                else:
                    break
    return lcss


def match_by_flipping(x, y):
    def condition():
        return len(lcss) > 0.7 * len(x)
    lcss = find_lcss(x, y)
    if condition():
        return x
    x_rev_comp = reverse_complement(x)
    lcss = find_lcss(x_rev_comp, y)
    if condition():
        return x_rev_comp
    x_comp = complement(x)
    lcss = find_lcss(x_comp, y)
    if condition():
        return x_comp
    x_rev = complement(x_rev_comp)
    return x_rev


def proper_align_to_ref(reference, read_refs_df):
    return read_refs_df.apply(lambda x: match_by_flipping(x, reference), axis=1)


def proper_align_read_refs_and_pkl():
    read_refs = load_read_refs_to_df()
    reference = get_reference()
    proper_aligned_refs = parallelize_dataframe_with_second_param(proper_align_to_ref, (reference, read_refs))
    # sanity check:
    lcss = find_lcss(reference, proper_aligned_refs)
    print(len(lcss), lcss)
    proper_aligned_refs.to_pickle('test/proper_aligned_refs.pkl')


proper_align_read_refs_and_pkl()
