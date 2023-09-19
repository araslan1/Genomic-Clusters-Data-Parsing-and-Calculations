import pandas as pd
import os
import numpy as np
import math
import csv


header = ['sample', '(0,0)', '(0,1)', '(0,2)', '(0,3)', '(0,4)', '(1,0)', '(1,1)', '(1,2)', '(1,3)', '(1,4)', '(2,0)', '(2,1)',
          '(2,2)', '(2,3)', '(2,4)', '(3,0)', '(3,1)', '(3,2)', '(3,3)', '(3,4)', '(4,0)', '(4,1)', '(4,2)', '(4,3)',
          '(4,4)']

samples = {'blank': {'compare_sample_dict': {}, 'ranges': {}}}


def parse_samples(sample_file):
    sample_df = pd.read_csv(sample_file, sep="\t")
    sample_ranges_to_be_merged = {}
    compare_sample_dict = {}
    event_mapping = {
        'Homozygous Copy Loss': -1,
        'CN Loss': -1,
        'High Copy Gain': 1,
        'CN Gain': 1,
    }

    for index, row in sample_df.iterrows():
        sample = row['Sample'][8:12]
        if sample not in samples:
            samples[sample] = {}
            samples[sample]['compare_sample_dict'] = {}
            samples[sample]['ranges'] = {}
        chrom_start_end = row['Chromosome Region'].split(':')
        chromosome = chrom_start_end[0]
        if chromosome not in samples[sample]['compare_sample_dict']:
            samples[sample]['compare_sample_dict'][chromosome] = []
        start_end = chrom_start_end[1].split('-')
        mystart = int(start_end[0].replace(',', ''))
        myend = int(start_end[1].replace(',', ''))
        s_value = 0
        if row['Event'] in event_mapping:
            s_value = event_mapping[row['Event']]
        if s_value != 0:
            temp_dict = {}
            temp_dict['Start'] = mystart
            temp_dict['End'] = myend
            temp_dict['S_val'] = s_value
            samples[sample]['compare_sample_dict'][chromosome].append(temp_dict)
        if s_value != 0:
            if chromosome not in samples[sample]['ranges']:
                samples[sample]['ranges'] [chromosome] = []
            samples[sample]['ranges'] [chromosome].append(int(mystart))
            samples[sample]['ranges'] [chromosome].append(int(myend))

def calculate_arm_lengths(input_file):
    chromosome_lengths = {}
    file = open(input_file, "r")
    prev_arm = "x"
    p_found = False
    q_found = False
    new_chromosome_found = False
    previous_index = 1
    prev_chromosome = ""
    when_file_ends = 0
    for line in file:
        if line == "\n":
            continue
        values = line.split('\t')
        chromosome = values[0]
        # chromosome = chromosome.strip("chr")
        start = int(values[1])
        end = int(values[2])
        when_file_ends = end
        arm = values[3]
        centromere = values[4]
        # print(centromere)
        if chromosome not in chromosome_lengths:
            chromosome_lengths[chromosome] = {}
            new_chromosome_found = True
        if arm[0] != prev_arm[0]:
            if new_chromosome_found:
                if not p_found:
                    chromosome_lengths[chromosome]['p_start'] = start
                    p_found = True
                else:
                    chromosome_lengths[prev_chromosome]['q_end'] = previous_index
                    chromosome_lengths[chromosome]['p_start'] = start
                    q_found = False
                new_chromosome_found = False
            elif not q_found:
                chromosome_lengths[chromosome]['p_end'] = previous_index
                chromosome_lengths[chromosome]['q_start'] = start
                q_found = True
        prev_chromosome = chromosome
        prev_arm = arm
        previous_index = end
    chromosome_lengths[prev_chromosome]['q_end'] = when_file_ends
    return chromosome_lengths
# Usage





def parse_profile_data(gains_file, losses_file, i, j, sample_number):
    # if path doesn't exist return
    if not os.path.exists(gains_file):
        return
    if not os.path.exists(losses_file):
        return
    # read cytoband file to calculate chromosome arm_lengths
    chromosome_lengths = calculate_arm_lengths("cytoband.txt")
    # create data frames from gains/losses file
    gains_df = pd.read_csv(gains_file, sep="\t")
    losses_df = pd.read_csv(losses_file, sep="\t")

    # create a data frame of the sample profile data
    # sample_profile_df = pd.read_csv('sampleprofile.csv')

    #   #   #   #
    # These variable basically store every line in the file by index
    # The nth index of each array below contains the chromosome, start, end, gain, loss, p_value of the nth line
    chromosomes = []
    start = []
    end = []
    gains = []
    losses = []
    p_values = []
    #   #   #   #

    # A new section should be made every time a p_value changes
    # Thus we will consolidate the gains/losses data, to merge consecutive sections with the same p_values
    # together and store only one big range for those consecutive sections
    # The sections start/end value will be store in the profile_ranges dictionary below
    consolidated_p_vals = []
    profile_ranges = {}

    # When a sample range is a nonzero value, the range should also be considered a new section
    # Thus, we will store the nonzero sample ranges here and merge them with the profile ranges dictionary above
    sample_ranges_to_be_merged = {}

    # These dictionaries are to quickly access the p_values/s_values for a specific range
    # for later down in the code and doesn't do any data manipulation
    compare_profile_dict = {}
    compare_sample_dict = {}

    # This dictionary is the final storage container after all the data has been parsed and filtered
    # It will store chromosome as it's key and an array of arrays with the desired variables
    ranges = {}

    # This is to initially read the gains/losses file and initialize compare_profile_dict and other variables
    for (index_losses, row_losses), (index_gains, row_gains) in zip(losses_df.iterrows(), gains_df.iterrows()):
        chromosome = row_losses['Chromosome Region']
        if chromosome not in compare_profile_dict:
            compare_profile_dict[chromosome] = []

        temp_dict = {}
        temp_dict['Start'] = row_losses['Start']
        temp_dict['End'] = row_losses['End']
        temp_dict['Loss'] = row_losses['loss']
        temp_dict['Gains'] = row_gains['gain']
        temp_dict['P_val'] = row_gains['gain'] - row_losses['loss']
        compare_profile_dict[chromosome].append(temp_dict)

        chromosomes.append(row_losses['Chromosome Region'])
        start.append(row_losses['Start'])
        end.append(row_losses['End'])
        losses.append(row_losses['loss'])
        gains.append(row_gains['gain'])

    #This is to go back and calculate the p_values
    for gain, loss in zip(gains, losses):
        p_values.append(gain - loss)

    # variable to traverse the chromosomes, start, end, gains, losses array
    # it's basically a line counter
    line = 0
    #

    # for chromosome in chromosomes:
    #     for item in compare_profile_dict[chromosome]:
    #         profile_ranges[chromosome] = []
    #         profile_ranges[chromosome].append(item['Start'])
    #         profile_ranges[chromosome].append(item['End'])
    #         consolidated_p_vals.append(item['P_val'])

    for chromosome, p_val in zip(chromosomes, p_values):
        # add chromosome to dictionary if it hasn't been found yet
        if chromosome not in profile_ranges:
            profile_ranges[chromosome] = []

        # This is to check for end of file since the last range value on the last line wasn't being added
        # Manually adds the last range value
        if line == len(start) - 1:
            profile_ranges[chromosome].append(end[line])
            consolidated_p_vals.append(p_val)
        # This occurs at any other line of the file
        elif line > 0:
            # Check if p_value changed, if it did check if it's the same or different chromosome
            # If it's a different chromosome go back and append the end of the previous chromosome
            if p_values[line - 1] != p_val:
                consolidated_p_vals.append(p_val)
                if chromosomes[line - 1] == chromosome:
                    profile_ranges[chromosome].append(end[line - 1])
                if chromosomes[line - 1] != chromosome:
                    profile_ranges[chromosomes[line - 1]].append(end[line - 1])
                if end[line - 1] != start[line]:
                    profile_ranges[chromosome].append(start[line])
        # This occurs at the start of the file
        elif line == 0:
            profile_ranges[chromosome].append(start[0])
            consolidated_p_vals.append(p_val)

        line += 1
    #

    # for chromosome in chromosome_lengths:
    #     if chromosome not in profile_ranges:
    #         profile_ranges[chromosome] = []
    #     profile_ranges[chromosome].append(chromosome_lengths[chromosome]['p_start']-1)
    #     profile_ranges[chromosome].append(chromosome_lengths[chromosome]['q_start'])
    #     profile_ranges[chromosome].append(chromosome_lengths[chromosome]['q_end'])


    # Parse the sample profile csv, initializes compare_sample_dict and sample_ranges_to_be_merged if s_value isn't 0

    # for index, row in sample_profile_df.iterrows():
    #     chromosome = row['chromosome']
    #     if chromosome not in compare_sample_dict:
    #         compare_sample_dict[chromosome] = []
    #     region = row['Range']
    #     s_value = row['change']
    #     region = region.strip('(')
    #     region = region.strip(')')
    #     region = region.split(", ")
    #     mystart = region[0]
    #     myend = region[1]
    #     if s_value != 0:
    #         temp_dict = {}
    #         temp_dict['Start'] = mystart
    #         temp_dict['End'] = myend
    #         temp_dict['S_val'] = s_value
    #         compare_sample_dict[chromosome].append(temp_dict)
    #     if s_value != 0:
    #         if chromosome not in sample_ranges_to_be_merged:
    #             sample_ranges_to_be_merged[chromosome] = []
    #         sample_ranges_to_be_merged[chromosome].append(int(mystart))
    #         if end not in sample_ranges_to_be_merged[chromosome]:
    #             sample_ranges_to_be_merged[chromosome].append(int(myend))
    #
    sample_ranges_to_be_merged = samples[sample_number]['ranges']
    compare_sample_dict = samples[sample_number]['compare_sample_dict']

    # merge profile ranges and sample ranges
    for chromosome in profile_ranges.keys():
        if chromosome not in sample_ranges_to_be_merged:
            combined_range = sorted(profile_ranges[chromosome])
            ranges[chromosome] = [combined_range]
        else:
            combined_range = sorted(list(set(profile_ranges[chromosome] + sample_ranges_to_be_merged[chromosome])))
            ranges[chromosome] = [combined_range]

    # consolidate range array into pair values in the (start_range, end_range) format
    for chromosome in ranges.keys():
        new_list = []
        curr_list = ranges[chromosome][0]
        index = 0
        for start in curr_list:
            if (index + 1 < len(curr_list)):
                new_list.append((start, curr_list[index + 1]))
            index += 1

        ranges[chromosome][0] = new_list

    # This for loop goes through and does all the desired calculations for each range and store them as an array
    # The ranges dictionary is a key:value pairing, with chromosome as the key and value as an array of arrays
    # The array of arrays has the following strucuture
    # Array[0] = ranges
    # Array[1] = p_values
    # Array[2] = s_values
    # Array[3] = m_values
    # Array[4] = (s-p)^2 values (denoted as distances)
    # Array[5] = f values
    for chromosome in ranges.keys():
        p_list = []
        s_list = []
        m_list = []

        # grab the arm lengths from the chromosome_lengths dictionary
        p_arm_start = chromosome_lengths[chromosome]['p_start']
        p_arm_end = chromosome_lengths[chromosome]['p_end']
        p_arm_length = float(p_arm_end - p_arm_start + 1)
        q_arm_start = chromosome_lengths[chromosome]['q_start']
        q_arm_end = chromosome_lengths[chromosome]['q_end']
        q_arm_length = float(q_arm_end - q_arm_start + 1)
        #   #   #   #

        for range in ranges[chromosome][0]:
            eventsize = float(range[1] - range[0])
            in_between = range[0] + 1

            if p_arm_start <= in_between < p_arm_end:
                m = float(eventsize) / float(p_arm_length)
                m_list.append(m)
            elif q_arm_start <= in_between < q_arm_end:
                m = float(eventsize) / float(q_arm_length)
                m_list.append(m)
            else:
                m_list.append(-0)
            found = False
            s_found = False
            if chromosome in compare_sample_dict:
                for s_item in compare_sample_dict[chromosome]:
                    compare_start = s_item['Start']
                    compare_end = s_item['End']
                    if float(compare_start) < in_between < float(compare_end):
                        s_found = True
                        s_list.append(s_item['S_val'] * 1.0)
                        break
                if not s_found:
                    s_list.append(0.0)
            else:
                s_list.append(0.0)

            if chromosome in compare_profile_dict:
                for item in compare_profile_dict[chromosome]:
                    compare_start = item['Start']
                    compare_end = item['End']
                    if float(compare_start) < in_between < float(compare_end):
                        found = True
                        p_list.append(item['P_val'] / 100.0)
                        break
                if not found:
                    p_list.append(0)
            else:
                p_list.append(0)
        ranges[chromosome].append(p_list)
        ranges[chromosome].append(s_list)
        ranges[chromosome].append(m_list)

        distances = []
        f_values = []
        for p_val, s_val in zip(ranges[chromosome][1], ranges[chromosome][2]):
            if s_val == 0:
                distances.append(abs(p_val))
            else:
                distances.append((pow(float(s_val - p_val), 2)))

            f_value = float(abs(float(s_val)) + abs(float(p_val)) + 0.01)
            f_values.append(f_value)
        ranges[chromosome].append(distances)
        ranges[chromosome].append(f_values)

    # The rest of the code below is to output the ranges dictionary in the form of csvs by chromosome number and print
    # the final distance value

    # used to count the chromosomes added
    index = 1
    for chromosome in ranges.keys():
        curr_range = ranges[chromosome][0]
        curr_p_vals = ranges[chromosome][1]
        curr_s_vals = ranges[chromosome][2]
        curr_m = ranges[chromosome][3]
        curr_distance_vals = ranges[chromosome][4]
        curr_f_vals = ranges[chromosome][5]
        np_m_vals = np.array(curr_m)
        np_f_vals = np.array(curr_f_vals)
        curr_w_vals = np.array(curr_f_vals) * np.array(curr_m)
        ranges[chromosome].append(curr_w_vals)

        # result_df = pd.DataFrame({
        #     'chromosome': f"chr{index}",
        #     'Range': curr_range,
        #     'S_value': curr_s_vals,
        #     'P_value': curr_p_vals,
        #     'f_value': curr_f_vals,
        #     'm': curr_m,
        #     'w_value': curr_w_vals,
        #     '(s-p)^2': curr_distance_vals,
        # })
        # # result_df.reset_index(drop=True)
        # filename = f"chr{index}.csv"
        # result_df.to_csv(filename, index=True)

        index += 1

    w_sum = 0
    for chromosome in ranges.keys():
        for w in ranges[chromosome][6]:
            w_sum += w

    numerator = 0
    for chromosome in ranges.keys():
        for distance, weight in zip(ranges[chromosome][4], ranges[chromosome][6]):
            top = float(distance) * float(weight)
            numerator += top

    finalDistance = math.sqrt(float(numerator) / float(w_sum))
    return round(finalDistance, 5)


def main():
    parse_samples("all_samples.txt")
    matrix = open('matrix.csv', 'w')
    # create the csv writer
    writer = csv.writer(matrix)

    # write a row to the csv file
    writer.writerow(header)
    for sample in samples:
        row = [sample]
        for i in range(5):
            for j in range(5):
                distance = parse_profile_data(f"({i},{j})_gains.txt", f"({i},{j})_losses.txt", f"{i}", f"{j}", sample)
                row.append(str(distance))
        writer.writerow(row)

    # for sample in samples:
    #     for i in range(5):
    #         for j in range(5):
    #             parse_profile_data(f"({i},{j})_gains.txt", f"({i},{j})_losses.txt", f"{i}", f"{j}")




main()
