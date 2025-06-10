import os
import subprocess
from absl import app, flags
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment

def run_jackhmmer(query_fasta, database_fasta, output_sto, iterations, jackhmmer_path= '/usr/bin/jackhmmer'):
    cmd = [
        jackhmmer_path,
        '-N', str(iterations),
        '-A', output_sto,
        query_fasta,
        database_fasta
    ]
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


from Bio.Seq import Seq

def remove_query_gaps(stockholm_file, fasta_output):
    aln = AlignIO.read(stockholm_file, "stockholm")
    query_record = aln[0]
    query_seq = str(query_record.seq).replace("-", "")

    aligned_records = []
    for record in aln:
        new_seq = []
        for q_char, s_char in zip(query_record.seq, record.seq):
            if q_char != "-":
                new_seq.append(s_char)
        record.letter_annotations = {}  #Clear per-letter annotations
        record.seq = Seq("".join(new_seq))  #Safely assign new sequence
        aligned_records.append(record)

    AlignIO.write(MultipleSeqAlignment(aligned_records), fasta_output, "fasta")
    print(f"Saved final MSA with ungapped query to {fasta_output}")


def convert_to_a3m(input_fasta, output_file):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    query_seq = records[0].seq
    with open(output_file, "w") as out:
        for record in records:
            a3m_seq = "".join([res if q != "-" else "" for q, res in zip(query_seq, record.seq)])
            out.write(f">{record.id}\n{a3m_seq}\n")




def main(argv):
    del argv  # unused
    tmp_dir ="/tmp/jackhmmer"
    msa_outfile = "msa.fasta"
    sto_file = os.path.join(tmp_dir, "jackhmmer_output.sto")
    FLAGS = flags.FLAGS
    run_jackhmmer(FLAGS.query_fasta, FLAGS.database_fasta, output_sto=sto_file, iterations = 3, tmp_dir=tmp_dir)
    remove_query_gaps(sto_file, msa_outfile)
    convert_to_a3m(msa_outfile, msa_outfile.replace(".fasta", ".a3m"))

if __name__ == "__main__":
    flags.mark_flag_as_required('query_fasta')
    flags.mark_flag_as_required('database_fasta')
    app.run(main)

