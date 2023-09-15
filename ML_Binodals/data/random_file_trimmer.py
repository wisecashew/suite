import numpy as np
import re
import os

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

if __name__=="__main__":

    # Get the current directory
    current_directory = os.getcwd()

    # List all files in the current directory
    all_files = os.listdir(current_directory)

    # Filter files with the ".binodal" extension
    binodal_files = [file for file in all_files if (file.startswith("vs") and file.endswith(".binodal"))]
    # binodal_files = ["vs_10.0-vc_1.0-vp_1.0-chisc_-1.0-chips_-2.0-chipc_-12.0.binodal"] # vs_8.0-vc_1.0-vp_1.0-chisc_-1.0-chips_-2.0-chipc_-11.0.binodal"]

    lines_to_keep = 40001
    for bfile in binodal_files:
        print(f"Going through {bfile}...")
        with open(bfile, 'rb') as fp:
            c_generator = _count_generator(fp.raw.read)
            # count each \n
            lcount = sum(buffer.count(b'\n') for buffer in c_generator)
        print(f"line_count = {lcount}...")

        if lcount < lines_to_keep:
            continue

        relevant_idx = np.arange(1, lcount, dtype=int)
        np.random.shuffle(relevant_idx)
        random_index_generator = relevant_idx[:40000]
        random_index_generator = np.array(random_index_generator, dtype=int)

        print (f"len(random_index_generator) = {len(random_index_generator)}...")

        inp = open(bfile, 'r')
        oup = open("trim."+bfile, 'w')
        for idx, line in enumerate(inp):
            if idx == 0:
                oup.write(line)
            elif idx in random_index_generator:
                oup.write(line)
            else:
                continue

        inp.close()
        oup.close()

        print (f"Done with {bfile}!")







