# GeneFinder

This project aims to provide some tools to aid in genetic archeology and understanding.
Given a 5-to-3 FASTA DNA sequence, a codon translation table in NCBI format, and a mutation profile in the custom .mprof mutation file format, it will attempt to aid in determining whether the sequence has been mutating around a similar configuration or a different one.
It will do this by using a 'robustness' number from 0-1, indicating the probability that a given codon frame will express the amino acid (or start/stop codon) in question after mutation, approximated by considering all mutations down to a given relative probability, and assuming the potential for infinite mutation.
Example mutation spectrum is a slightly edited version of the data presented in Lee et al.'s "Rate and molecular spectrum of spontaneous mutations in the bacterium Escherichia coli as determined by whole-genome sequencing", [doi:10.1073/pnas.1210309109](http://www.pnas.org/content/109/41/E2774)

Anyone may use this tool for non-commercial use. This tool is provided without any warranty or guarantee.

For commercial inquiries and feedback; Contact e-mail: brian3.14159@gmail.com
