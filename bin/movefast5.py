import os
import re
import argparse

from shutil import copyfile

parser = argparse.ArgumentParser(description='Split a dataset.')
parser.add_argument('source', help='source directory')
parser.add_argument('dest', help='destination directory')
args = parser.parse_args()

for subdir, dirs, files in os.walk(args.source):
    for file in files:
        if re.search('.txt', file):
            continue
        new_filename = subdir[subdir.rfind('/') + 1:] + '.fast5'
        copyfile(os.path.join(subdir, file), os.path.join(args.dest, new_filename))
