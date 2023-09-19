# Calculate p and q arm lengths per chromosome
# input_file should be 'cytoband.txt'
# returns a dictionary of chromosome lengths (specifically p_start, p_end, q_start, q_end)


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
        if centromere == 'acen':
            continue
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
    print(chromosome_lengths)
    return chromosome_lengths
# Usage


calculate_arm_lengths("cytoband.txt")