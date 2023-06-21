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


def count_sequences(sequences_fpath):
    records = SeqIO.parse(sequences_fpath, "fasta")
    return len([record for record in records])


def get_sequences_info(sequences_fpath):
    records = SeqIO.parse(sequences_fpath, "fasta")
    sequences_info = {}
    for record in records:
        sequences_info[record.id] = len(record.seq)

    return sequences_info

def get_sorted_info(sequences_info):
    sorted_info =  dict(sorted(sequences_info.items(), key= lambda item: item[1]))
    return sorted_info

def get_total_length(sorted_info):
    total_length = sum(sorted_info[id] for id in sorted_info)
    return total_length

def get_shortest_sequence_id(sorted_info):
    shortest_seq_id = list(sorted_info.keys())[0]
    return shortest_seq_id

def get_longest_sequence_id(sorted_info):
    longest_seq_id = list(sorted_info.keys())[-1]
    return longest_seq_id

def get_average_length(sorted_info):
    num_of_seqs = len(sorted_info)
    total_length = get_total_length(sorted_info)
    average_length = total_length/num_of_seqs
    return average_length


def main():
    arguments = get_options()
    number_of_sequences = count_sequences(Path(arguments.sequences))
    num_msg = "Number of sequences: {}"
    print(num_msg.format(number_of_sequences))

    sequences_info = get_sequences_info(Path(arguments.sequences))
    sorted_info = get_sorted_info(sequences_info)
    total_length = get_total_length(sorted_info)
    total_msg = "Total length: {}"
    print(total_msg.format(total_length)) 

    shortest_seq_id = get_shortest_sequence_id(sorted_info)
    shortest_seq_length = sorted_info[shortest_seq_id]
    shortest_msg = "Shortest sequence: {} {}"
    print(shortest_msg.format(shortest_seq_length, shortest_seq_id))

    longest_seq_id = get_longest_sequence_id(sorted_info)
    longest_seq_length = sorted_info[longest_seq_id]
    longest_msg = "Longest sequence: {} {}"
    print(longest_msg.format(longest_seq_length, longest_seq_id))

    average_length = get_average_length(sorted_info)
    average_msg = "Average length: {}"
    print(average_msg.format(average_length))

if __name__ == "__main__":
    main()