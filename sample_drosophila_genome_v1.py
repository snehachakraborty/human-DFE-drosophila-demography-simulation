import os
import stdpopsim
import numpy as np
import pandas as pd

os.chdir("/u/home/s/sneha_c/project-klohmuel/meixi_simulations/rescaled_1Mb_no_recom_map")


def check_intervals_overlap(starts, ends):
    # intervals should be sorted already
    for i in range(len(starts) - 1):
        if starts[i+1] < ends[i]:
            print(i)
            raise ValueError("Overlapping intervals")


# bedtools subtract -a $statebed -b $statequi
# trivial and obvious (not)
# this function DOES NOT allow overlapping intervals in a
def sub_intervals(xs, ys):
    if len(xs) == 0:
        return []

    i = j = 0
    rs = []
    cur_lo = xs[0][0]

    while i < len(xs) and j < len(ys):
        x_lo, x_hi = xs[i]
        cur_lo = max(cur_lo, x_lo)
        y_lo, y_hi = ys[j]

        if cur_lo < y_lo:
            rs.append((cur_lo, min(y_lo, x_hi)))
            if x_hi < y_lo:
                i += 1
            else:
                cur_lo = y_lo
        else:
            if y_hi < x_hi:
                cur_lo = y_hi
                j += 1
            else:
                i += 1

    if i < len(xs):
        rs.append((max(cur_lo, xs[i][0]), xs[i][1]))
    rs.extend(xs[i+1:])
    return rs


def merge_intervals(intervals, sort=True):
    if len(intervals) == 0:
        return []

    if not sort:
        intervals = sorted(intervals, key=lambda x: x[0])

    rs = []
    cur_lo, cur_hi = intervals[0]
    for lo, hi in intervals[1:]:
        if lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            rs.append((cur_lo, cur_hi))
            cur_lo, cur_hi = lo, hi

    rs.append((cur_lo, cur_hi))
    return rs


def fill_in_annots(annots, fill_name, fill_col, end_val):
    """
    >>> df = pd.DataFrame({'start': [1, 5, 10, 15], 'end': [3, 9, 15, 21], 'type': ['exon']*4})
    >>> fill_in_annots(df, 'fill', 'type', 25)
       start  end   type
    0      0    1   fill
    1      1    3   exon
    2      3    5   fill
    3      5    9   exon
    4      9   10   fill
    5     10   15   exon
    6     15   21   exon
    7     21   25   fill
    """
    new_annots = []
    if len(annots) == 0:
        return pd.DataFrame([(0, end_val, fill_name)], columns=['start', 'end', fill_col])
    if annots.iloc[0, 0] != 0:
        new_annots.append((0, annots.iloc[0, 0], fill_name))
    for i in range(len(annots) - 1):
        s, e, f = annots.iloc[i]
        new_annots.append((s, e, f))
        next_s = annots.iloc[i + 1, 0]
        if e < next_s:
            new_annots.append((e, next_s, fill_name))
    new_annots.append((annots.iloc[-1, 0], annots.iloc[-1, 1], annots.iloc[-1, 2]))
    if annots.iloc[-1, 1] != end_val:
        new_annots.append((annots.iloc[-1, 1], end_val, fill_name))
    return pd.DataFrame(new_annots, columns=['start', 'end', fill_col])


# run this once
def make_chrom_annots():
    if all([os.path.exists(f"./sim_files/{chrom}_annotations.csv") for chrom in chr_names]):
        print("Annotations already exist")
        return

    os.makedirs("./sim_files", exist_ok=True)

    # CDS annots from stdpopsim
    cds = species.get_annotations('FlyBase_BDGP6.32.51_CDS')

    # for each chromosome, save all types of annotations to file
    for chrom in chr_names:
        # get annotations for each type for this chromosome
        cds_chrom = cds.get_chromosome_annotations(chrom)

        # check for overlaps within each type of annotation
        check_intervals_overlap(cds_chrom[:, 0], cds_chrom[:, 1])

        # subtraction
        # cds_chrom = sub_intervals(cds_chrom, centromeres_chrom.values)
        # cds_chrom = sub_intervals(cds_chrom, gaps_chrom.values)

        cds_chrom = pd.DataFrame(cds_chrom, columns=['start', 'end'])

        # add column to each dataframe for each type of annotation
        cds_chrom['type'] = 'exon'
        # centromeres_chrom['type'] = 'exclude'
        # gaps_chrom['type'] = 'exclude'

        # combine all annotations for this chromosome and sort by start and end
        chrom_annots = pd.concat([cds_chrom])
        chrom_annots = chrom_annots.sort_values(by=['start', 'end'])

        # check for overlaps
        if check_intervals_overlap(chrom_annots['start'].values, chrom_annots['end'].values):
            print(f"Overlapping intervals in {chrom}")
            continue

        # fill in gaps between annotations
        chrom_annots = fill_in_annots(chrom_annots, 'bkgd', 'type', chr_lengths[chrom])

        # save to file
        chrom_annots.to_csv(f"./sim_files/{chrom}_annotations.csv", index=False)


def subset_annotations(chrom_annots, start, end):
    """
    >>> df = pd.DataFrame({'start': [1, 3, 5, 9, 10, 15], 'end': [3, 5, 9, 10, 15, 21]})
    >>> subset_annotations(df, 2, 9)
       start  end
    0      2    3
    1      3    5
    2      5    9
    >>> subset_annotations(df, 2, 11)
       start  end
    0      2    3
    1      3    5
    2      5    9
    3      9   10
    4     10   11
    >>> subset_annotations(df, 10, 21)
       start  end
    2     10   15
    3     15   21
    """
    # get the left bin index for start and end
    start_idx = np.digitize(start, chrom_annots['start'].values) - 1
    start_idx = np.max([0, start_idx])
    end_idx = np.digitize(end, chrom_annots['end'].values)
    if chrom_annots.iloc[end_idx - 1, 1] == end:
        end_idx -= 1
    subset_annots = chrom_annots.iloc[start_idx:end_idx + 1].copy()
    subset_annots.iloc[0, 0] = start
    subset_annots.iloc[-1, 1] = end
    return subset_annots


def sample_chrom_segment(chr_lengths, chr_names, chr_probs, length=1_000_000):
    chrom = np.random.choice(chr_names, p=chr_probs)
    chrom_length = chr_lengths[chrom]
    start = np.random.randint(1, chrom_length - length)
    end = start + length
    return chrom, start, end


#%%
# chromosome details from stdpopsim
# chromosome details from stdpopsim
species = stdpopsim.get_species("DroMel")
chr_lengths = {f'chr{i.id}': i.length for i in species.genome.chromosomes if i.id not in ["4", "X", "Y", "mitochondrion_genome"]}
chr_probs = np.array(list(chr_lengths.values())) / sum(chr_lengths.values())
chr_names = list(chr_lengths.keys())

make_chrom_annots()

#%%
# sample 1 Mb segments for 108 simulations
np.random.seed(3079104)
segment_details = []
for i in range(1, 109):
    while True:
        # sample chromosome segment
        c, s, e = sample_chrom_segment(chr_lengths, chr_names, chr_probs)

        # read annots of segment
        chrom_annots = pd.read_csv(f"./sim_files/{c}_annotations.csv")
        segment_annots = subset_annotations(chrom_annots, s, e)

        # make sure there are no "exclude"s in segment annots
        if "exclude" in segment_annots['type'].values:
            continue
        else:
            break

    # make the annots start at 0 and end at length
    segment_annots['start'] = segment_annots['start'] - s
    segment_annots['end'] = segment_annots['end'] - s - 1  # end is inclusive in slim
    segment_annots.to_csv(f'./sim_files/sim_{i}_annots.csv', index=False)

    # save segment details
    segment_details.append((c, s, e))

# save segment details to file
segment_details_df = pd.DataFrame(segment_details, columns=['chrom', 'start', 'end'])
segment_details_df.to_csv(f'./sim_files/sim_segment_details.csv', index=False)
