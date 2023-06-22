import argparse
from pathlib import Path

from Bio import SeqIO


def construct_arguments():
    desc = "Get sequence stats"
    parser = argparse.ArgumentParser(description=desc)
    help_reference = "Input sequences to calculate stats"
    parser.add_argument("--sequences" , 
                        "-s", type=str,
                        help=help_reference,
                        required=True)
    help_output = "Output file"
    parser.add_argument("--out", "-o",
                        type=str, help=help_output,
                        required=True)
    return parser


def get_options():
    parser = construct_arguments()
    return parser.parse_args()
    # sequences_fpath = Path(options.sequences)
    # output_fpath = Path(options.out)

    # return {'sequences': sequences_fpath,
    #         'output_fpath': output_fpath}

def get_sequences_info(sequences_fpath):
    records = SeqIO.parse(sequences_fpath, "fasta")
    sequences_info = {}
    for record in records:
        sequences_info[record.id] = len(record.seq)

    return sequences_info

def get_sorted_info(sequences_info):
    sorted_info =  dict(sorted(sequences_info.items(), key= lambda item: item[1]))
    return sorted_info

def count_sequences(sequences_info):
    return len(sequences_info.keys())

def get_total_length(sorted_info):
    total_length = sum(length for length in sorted_info.values())
    return total_length

def get_shortest_and_longest(sorted_info):
    sorted_info_list = list(sorted_info.items())
    shortest_seq = sorted_info_list[0]
    longest_seq = sorted_info_list[-1]
    return shortest_seq, longest_seq

def get_average_length(sorted_info):
    num_of_seqs = len(sorted_info)
    total_length = get_total_length(sorted_info)
    average_length = total_length/num_of_seqs
    return average_length

def get_n(sorted_info, total_length, n=50):
    rev_sorted_info = dict(reversed(list(sorted_info.items())))
    n_length = total_length * n/100
    length_sum = 0

    for id in rev_sorted_info:
        length_sum += rev_sorted_info[id]
        if length_sum >= n_length:
            return [rev_sorted_info[id], id]

def main():
    arguments = get_options()
    sequences_info = get_sequences_info(Path(arguments.sequences))
    sorted_info = get_sorted_info(sequences_info)

    number_of_sequences = count_sequences(sorted_info)

    total_length = get_total_length(sorted_info)

    shortest_seq, longest_seq = get_shortest_and_longest(sorted_info)

    average_length = get_average_length(sorted_info)

    with open(arguments.out, "w") as out_fhand:
        msg = "Number of sequences:\t{}"
        out_fhand.write(msg.format(number_of_sequences) + "\n")

        msg = "Total length:\t{}"
        out_fhand.write(msg.format(total_length) + "\n") 

        msg = "Shortest sequence:\t{}\t{}"
        out_fhand.write(msg.format(shortest_seq[0],shortest_seq[1]) + "\n")

        msg = "Longest sequence:\t{}\t{}"
        out_fhand.write(msg.format(longest_seq[0],longest_seq[1]) + "\n")
        
        msg = "Average length:\t{}"
        out_fhand.write(msg.format(average_length) + "\n")

        for number in [95, 90, 75, 50, 25]:
            nXX = get_n(sorted_info, total_length, n=number)
            msg = "N{}:\t{}\t{}"
            out_fhand.write(msg.format(number, nXX[0], nXX[1]) + "\n")
        
        out_fhand.flush()

if __name__ == "__main__":
    main()