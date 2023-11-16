import random
from Bio.Seq import Seq
from demux_read_index import load_sample_sheet
import click

mutate = {
    'A': list('CGT'),
    'C': list('AGT'),
    'G': list('ACT'),
    'T': list('ACG')
}
def add_mismatches(seq, n):
    for i in random.sample(range(len(seq)), k=n):
        seq = seq[:i] + random.choice(mutate[seq[i]]) + seq[i+1:]
    return seq

@click.command()
@click.argument("sample_sheet", type=click.File('rt'))
@click.argument("r1_file", type=click.File('wt'))
@click.argument("r2_file", type=click.File('wt'))
@click.option("--reads_per_sample", type=int, default=10)
@click.option("--insert_length", type=int, default=300)
@click.option("--read_length", type=int, default=150)
@click.option("--bc_mismatches", type=int, default=0)
def main(sample_sheet, r1_file, r2_file, reads_per_sample, insert_length, read_length, bc_mismatches):
    samples = load_sample_sheet(sample_sheet.name, ".")
    gen_reads = [ s for s in samples for i in range(reads_per_sample)]
    random.shuffle(gen_reads)
    insert = "".join(random.choices(list("ACGT"), k=insert_length))
    insert_rc = str(Seq(insert).reverse_complement())
    qscores = "~" * read_length

    for i, s in enumerate(gen_reads):
        r1 = (add_mismatches(s.f_idx.seq, bc_mismatches) + s.f_idx.spacer + insert)[:read_length]
        r2 = (add_mismatches(s.r_idx.seq, bc_mismatches) + s.r_idx.spacer + insert_rc)[:read_length]
        r1_file.write(f"@{i}:{s.sample_id}\n{r1}\n+{i}:{s.sample_id}\n{qscores}\n")
        r2_file.write(f"@{i}:{s.sample_id}\n{r2}\n+{i}:{s.sample_id}\n{qscores}\n")

main()

