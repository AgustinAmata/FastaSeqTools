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

def get_n(sorted_info, n=50):
    rev_sorted_info = dict(reversed(list(sorted_info.items())))
    total_length = get_total_length(sorted_info)
    n_length = total_length * n/100
    length_sum = 0

    for id in rev_sorted_info:
        length_sum += rev_sorted_info[id]
        if length_sum >= n_length:
            return [rev_sorted_info[id], id]
            break

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

    # n95 = get_n(sorted_info, n=95)
    # n95_msg = "N95: {} {}"
    # print(n95_msg.format(n95[0], n95[1]))

    # n90 = get_n(sorted_info, n=90)
    # n90_msg = "N90: {} {}"
    # print(n90_msg.format(n90[0], n90[1]))

    # n75 = get_n(sorted_info, n=75)
    # n75_msg = "N75: {} {}"
    # print(n75_msg.format(n75[0], n75[1]))

    # n50 = get_n(sorted_info)
    # n50_msg = "N50: {} {}"
    # print(n50_msg.format(n50[0], n50[1]))

    # n25 = get_n(sorted_info, n=25)
    # n25_msg = "N25: {} {}"
    # print(n25_msg.format(n25[0], n25[1]))

    for number in [95, 90, 75, 50, 25]:
        nXX = get_n(sorted_info, n=number)
        nXX_msg = "N{}: {} {}"
        print(nXX_msg.format(number, nXX[0], nXX[1]))

if __name__ == "__main__":
    main()