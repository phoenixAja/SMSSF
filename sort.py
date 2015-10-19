#!/usr/bin/python

import glob
import os
import sys

size = int(sys.argv[1])
def mk_tree(path):
    #files to sort
    files = glob.glob(os.path.join(path, "*.mol2"))
    #creating directories to be filled with certain number of decoy .mol2 files
    chunks = [files[chunk:chunk+10] for chunk in range(0, len(files), size)]
    for i, chunk in enumerate(chunks):
        new_dir = os.path.join(path, "decoydum%05d" % i)
        os.mkdir(new_dir)
        for fn in chunk:
            os.rename(fn, os.path.join(new_dir, os.path.basename(fn)))

def main():
    mk_tree("decoys")

if __name__ == '__main__':
    main()
