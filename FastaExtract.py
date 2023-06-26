import argparse
import sys
from csv import DictReader
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqUtils import gc_fraction

def argument_parser():
    desc = """Extract sequences from a fasta file.
    An additional file with specific ids along with other data can be provided.
    It must follow the order: ID + start + end + orientation + new ID;
    start and end must be 1-based"""
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

    help_extract_by_length = "Integer to extract sequences by length."
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

    help_revcom = "Returns the reverse complementary of the sequences"
    parser.add_argument("--revcom", "-R", help=help_revcom,
                        action="store_true", required=False)

    help_cut= "Cuts the sequences given the coordinates in the id file"
    parser.add_argument("--cut", "-C", help=help_cut,
                        action="store_true", required=False)
    
    help_orientate= "Orientates the sequence given the strand specified in the id file"
    parser.add_argument("--orientate", "-O", help=help_orientate,
                        action="store_true", required=False)

    help_replace_ids = "Replaces ids with the ones specified in the id file"
    parser.add_argument("--replace", "-I", help=help_replace_ids,
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

def get_length_sequences(sequences_hand, length):
    length_sequences = {}
    for id in sequences_hand.keys():
        if len(sequences_hand[id]) >= length:
            length_sequences[id] = sequences_hand[id]

    if length_sequences:
        return length_sequences
    else:
        sys.exit("No sequence passed the length threshold, please select a lower one")

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
        sys.exit("No sequence passed the composition threshold, please select a lower one")

def get_cut_positions(fid_list):
    cut_positions = {}
    for fid in fid_list:
        fid = list(fid.values())
        id = fid[0]
        trim_pos = fid[1:3]

        cut_positions[id] = trim_pos

    return cut_positions

def get_cut_sequences(sequences_hand, cut_positions):
    cut_sequences = {}
    for id in sequences_hand.keys():
        start = int(cut_positions[id][0]) - 1
        end = int(cut_positions[id][1])
        cut_seq = sequences_hand[id][start:end]
        cut_sequences[id] = cut_seq

    return cut_sequences

def get_orientation(fid_list):
    orientation = {}
    for line in fid_list:
        tags = list(line.values())
        id = tags[0]
        orient = tags[3]
        orientation[id] = orient
    
    return orientation

def get_orientated_sequences(sequences_hand, orientation):
    orientated_sequences = {}
    for id in sequences_hand.keys():
        sequence = sequences_hand[id]
        if orientation[id] == "-":
            revcom_seq = reverse_complement(sequence)
            orientated_sequences[id] = revcom_seq
        
        else:
            orientated_sequences[id] = sequence

    return orientated_sequences

def get_new_ids(fid_list):
    new_ids = {}
    for line in fid_list:
        tags = list(line.values())
        old_id = tags[0]
        new_id = tags[4]
        new_ids[old_id] = new_id
    
    return new_ids

def get_new_id_sequences(sequences_hand, new_ids):
    new_id_sequences = {}
    for id in sequences_hand:
        new_id = new_ids[id]
        sequence = sequences_hand[id]
        new_id_sequences[new_id] = sequence

    return new_id_sequences

def get_revcom_sequences(sequences_hand):
    revcom_sequences = {}
    for id in sequences_hand.keys():
        revcom = reverse_complement(sequences_hand[id])
        revcom_sequences[id] = revcom

    return revcom_sequences

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

    if arguments.length:
        length = arguments.length
        fasta_sequences = get_length_sequences(fasta_sequences, length)

    if arguments.id:
        fid_list = []
        for line in DictReader(open(arguments.id), delimiter="\t"):
            fid_list.append(line)

        selected_ids = get_selected_ids(fid_list)
        selected_sequences = get_selected_sequences(fasta_sequences,
                                                    selected_ids)
        
        if arguments.cut:
                cut_positions = get_cut_positions(fid_list)
                selected_sequences = get_cut_sequences(selected_sequences,
                                                       cut_positions)

        if arguments.orientate:
            orientation = get_orientation(fid_list)
            selected_sequences = get_orientated_sequences(selected_sequences,
                                                          orientation)

        if arguments.replace:
            new_ids = get_new_ids(fid_list)
            selected_sequences = get_new_id_sequences(selected_sequences,
                                                   new_ids)
            
        final_sequences = selected_sequences
        
    else:
        final_sequences = fasta_sequences

    if arguments.revcom:
        final_sequences = get_revcom_sequences(final_sequences)

    if arguments.split:
        construct_individual_outputs(final_sequences)

    else:
        with open(arguments.output, "w") as out_fhand:
            for id in final_sequences:
                out_fhand.write(f">{id}\n{final_sequences[id]}\n")
            
            out_fhand.flush()

if __name__ == "__main__":
    main()