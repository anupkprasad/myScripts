import re
import argparse

def process_sequence(line):

    new_line = []
    inside_parentheses = False

    for char in line:
        if char == '(':
            inside_parentheses = True
            new_line.append(char)  # Keep the parentheses but remove later
        elif char == ')':
            inside_parentheses = False
            new_line.append(char)  # Keep the parentheses but remove later
        elif inside_parentheses:
            new_line.append(char.lower())  # Make lowercase inside parentheses
        else:
            new_line.append(char)  # Keep original outside parentheses

    # Join processed characters and remove any unwanted characters
    processed_sequence = ''.join(new_line)
    cleaned_sequence = re.sub(r'[\}\{\]\[\)\(\*]', '', processed_sequence)

    return cleaned_sequence

def remove_mma_patterns(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Skip lines starting with '$'
            if line.startswith('$'):
                continue
            # Directly write lines that start with '>'
            if line.startswith('>'):
                outfile.write(line)
                continue
            # Process lines that start with '{'
            if line.startswith('{'):
                processed_sequence = process_sequence(line)
                outfile.write(processed_sequence + '\n')



prefixes = [
    'EGFR',
    'EGFR_erbb',
    'tk'
]


path = "/mnt/mydrive/alpfaFold3/dynamics/Custom_MSA/efgr_egfrFamily/"

for file in prefixes:
    print(path + file)
    remove_mma_patterns(path + file + '.fasta_A.mma', path + file + '.a3m')


jhdfkjd
def main():
    parser = argparse.ArgumentParser(description='Convert an mma or cma file to an a3m file')
    parser.add_argument('input_file', type=str, help='Path to input file')
    parser.add_argument('output_file', type=str, help='Path to output file')

    args = parser.parse_args()

    remove_mma_patterns(args.input_file, args.output_file)

if __name__ == "__main__":
    main()

