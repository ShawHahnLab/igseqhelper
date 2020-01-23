# Example Sequence Pairs

Trying to get a handle on how the read pairs are arranged.

Example #3 has short enough reads to recognize the extra material at the ends
of R1 and R2, but long enough to see how they align.

`P5_Graft` and `P5_seq` come just before R1.  R1 starts with the varying-length
randomized nucleotides, the barcode, and then the constant primer sequence.
Looking along the same 5' to 3' direction, R2 begins (at the right) just before
the C primer, followed by the reverse barcode (captured in I1), and finally the
Illumina P7 primer.

To trim:

 * R1: we should always stop *before* the C primer, so we can use that as 3'
   adapter.
 * R2: (looking 3' to 5') we should always stop *before* the P5 primer site, so
   we can use that as the 5' adapter.

For example:

    cutadapt -a TCCACCAAGGGCCCATCGGTCTTCCCCCTGGC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGA -o R1.trimmed.fastq -p R2.trimmed.fastq R1.fastq R2.fastq

With this, R2 should extend no farther than the beginning of R1, and vice versa.

<https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-11.pdf>
