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

def sto_to_fasta(stockholm_file, fasta_output):
    aln = AlignIO.read(stockholm_file, "stockholm")
    with open("output.fasta", "w") as out_handle:
        AlignIO.write(aln, fasta_output, "fasta")
    print(f"Saved final MSA with ungapped query to {fasta_output}")

def fasta_to_GapRemoved_a3m(input_fasta, output_file):
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
    sto_to_fasta(sto_file, msa_outfile)
    fasta_to_GapRemoved_a3m(msa_outfile, msa_outfile.replace(".fasta", ".a3m"))

if __name__ == "__main__":
    flags.mark_flag_as_required('query_fasta')
    flags.mark_flag_as_required('database_fasta')
    app.run(main)

