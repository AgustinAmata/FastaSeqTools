import argparse
import sys
from csv import DictReader
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def argument_parser():
    desc = """Extract sequences from a fasta file.
    An additional file with specific ids with the following structures can be provided:
    one-column (id), three columns (id + start + end)
    or four columns (id + start + end + new id); start and end must be 1-based"""
    parser = argparse.ArgumentParser(description=desc)

    help_input_fasta = "Input fasta file (mandatory)"
    parser.add_argument("--file", "-f", help=help_input_fasta,
                        required=True)

    help_output = "Output file name (default will create stdout)"
    parser.add_argument("--output", "-o", help=help_output,
                        default="stdout.fasta", required=False)

    help_extract_by_id = "Filename to extract sequences by id"
    parser.add_argument("--id", "-i", help=help_extract_by_id,
                        required=False)

    help_extract_by_length = """Integer to extract sequences by length.
    If --id is present, it will only work with single-column files (only the id)"""
    parser.add_argument("--length", "-l", help=help_extract_by_length,
                        type=int, required=False)

    help_extract_by_composition = """Extract sequences by
    minimum nucleotide content of the complete sequence in the fasta file. Example: GC=0.5"""
    parser.add_argument("--composition", "-c", help=help_extract_by_composition,
                        required=False)

    help_split_output = """Output is split into individual files
    (the name of each output is the id of the sequence contained in the file)"""
    parser.add_argument("--split", "-S", help=help_split_output,
                        action="store_true", required=False)

    return parser

def get_options():
    parser = argument_parser()
    return parser.parse_args()

def get_fasta_sequences(fasta_fhand):
    records = SeqIO.parse(fasta_fhand, "fasta")
    fasta_sequences = {}
    for record in records:
        fasta_sequences[record.id] = record.seq

    return fasta_sequences

def get_selected_ids(fid_list):
    selected_ids = [list(tags.values())[0] for tags in fid_list]

    return selected_ids

def get_selected_sequences(fasta_sequences, selected_ids):
    selected_sequences = {}
    for id in selected_ids:
        if id in fasta_sequences.keys():
            selected_sequences[id] = fasta_sequences[id]

    return selected_sequences

def get_trimming_positions(fid_list):
    trimming_positions = {}
    for fid in fid_list:
        fid = list(fid.values())
        id = fid[0]
        trim_pos = fid[1:3]

        trimming_positions[id] = trim_pos

    return trimming_positions

def get_new_ids(fid_list):
    new_ids = {}
    for line in fid_list:
        tags = list(line.values())
        old_id = tags[0]
        new_id = tags[3]
        new_ids[old_id] = new_id
    
    return new_ids

def get_cut_sequences(sequences_hand, length):
    cut_sequences = {}
    for id in sequences_hand.keys():
        cut_sequences[id] = sequences_hand[id][:length]

    return cut_sequences

def get_composition_sequences(sequences_hand, composition_hand):
    composition_sequences = {}

    if composition_hand.startswith("GC="):
        composition = float(composition_hand.strip("GC="))
        for id in sequences_hand.keys():
            sequence = sequences_hand[id]
            seq_GC_percentage = gc_fraction(sequence)
            
            if seq_GC_percentage >= composition:
                composition_sequences[id] = sequence

    elif composition_hand.startswith("AT="):
        composition = float(composition_hand.strip("AT="))
        for id in sequences_hand.keys():
            sequence = sequences_hand[id]
            seq_AT_percentage = 1 - gc_fraction(sequence)
            
            if seq_AT_percentage >= composition:
                composition_sequences[id] = sequence

    if composition_sequences:
        return composition_sequences
    else:
        sys.exit("No sequence passed the composition threshold, try a lower one")

def get_trimmed_sequences(sequences_hand, trimming_positions):
    trimmed_sequences = {}
    for id in sequences_hand.keys():
        start = int(trimming_positions[id][0]) - 1
        end = int(trimming_positions[id][1])
        trimmed_seq = sequences_hand[id][start:end]
        trimmed_sequences[id] = trimmed_seq

    return trimmed_sequences

def get_new_id_sequences(sequences_hand, new_ids):
    new_id_sequences = {}
    for id in sequences_hand:
        new_id = new_ids[id]
        sequence = sequences_hand[id]
        new_id_sequences[new_id] = sequence

    return new_id_sequences

def construct_individual_outputs(sequences_hand):
    for id in sequences_hand.keys():
        sequence = sequences_hand[id]

        with open(f"{id}.fasta", "w") as out_fhand:
            out_fhand.write(f">{id}\n{sequence}\n")
            out_fhand.flush()

def main():
    arguments = get_options()
    fasta_sequences = get_fasta_sequences(Path(arguments.file))

    if arguments.composition:
        fasta_sequences = get_composition_sequences(fasta_sequences,
                                                     arguments.composition)

    if arguments.id:
        fid_list = []
        for line in DictReader(open(arguments.id), delimiter="\t"):
            fid_list.append(line)

        line_length = len(fid_list[0].keys())
        selected_ids = get_selected_ids(fid_list)
        selected_sequences = get_selected_sequences(fasta_sequences,
                                                    selected_ids)
        
        if line_length == 1:
            if arguments.length:
                length = arguments.length
                cut_sequences = get_cut_sequences(selected_sequences, length)
                final_sequences = cut_sequences

            else:
                final_sequences = selected_sequences

        elif line_length == 3:
            trimming_positions = get_trimming_positions(fid_list)
            final_sequences = get_trimmed_sequences(selected_sequences,
                                                    trimming_positions)

        elif line_length == 4:
            trimming_positions = get_trimming_positions(fid_list)
            new_ids = get_new_ids(fid_list)

            trimmed_sequences = get_trimmed_sequences(selected_sequences,
                                                      trimming_positions)
            final_sequences = get_new_id_sequences(trimmed_sequences,
                                                   new_ids)

    else:
        if arguments.length:
            length = arguments.length
            cut_sequences = get_cut_sequences(fasta_sequences, length)
            final_sequences = cut_sequences

        else:
            final_sequences = fasta_sequences

    if arguments.split:
        construct_individual_outputs(final_sequences)

    else:
        with open(arguments.output, "w") as out_fhand:
            for id in final_sequences:
                out_fhand.write(f">{id}\n{final_sequences[id]}\n")
            
            out_fhand.flush()

if __name__ == "__main__":
    main()