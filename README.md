# TetraBin-Nucleotide-Encoding
A method for compressing FASTA files (or other Nucleotide raw data) into a more size efficient file

FASTA FILE ENCODING - ASCII:
In ASCII encoding, each character is represented by a single byte, which is a sequence of 8 bits. Since there are only 128 characters in the ASCII character set, each character can be represented by a unique combination of 7 bits (2^7 = 128).

For example, the ASCII encoding of the letter "A" is the binary sequence 01000001. This corresponds to the decimal value 65, which is the ASCII code for the letter "A". Similarly, the ASCII encoding of the letter "B" is the binary sequence 01000010, which corresponds to the decimal value 66.

Since each character is represented by a single byte in ASCII encoding, the size of an ASCII-encoded file is directly proportional to the number of characters in the file. For example, a text file containing 1000 ASCII characters will have a file size of 1000 bytes.

BINARY TETRAD ENCODING:
Tetrad	Letter	Representative Nucleotides
1111	N	
0001	A	
0010	T	
0010	U	
0100	C	
1000	G	
0011	W	A/T
1100	S	C/G
0101	M	A/C
1010	K	G/T
1001	R	A/G
0110	Y	C/T
1110	B	Not A
1011	D	Not C
0111	H	Not G
1101	V	Not T

COMPARISON OF ENCODING:
Sequence of "ATCGATCG" would normally get encoded in ASCII as "0100000101010100010000110100011101000001010101000100001101000111" (64 bits)

With the proposed encoding the same sequence would be "00010010010010000001001001001000" (32 bits)

AREAS OF EXPANSION (WIP):
.md5 checksums

FILE DESCRIPTIONS (WIP):
.fasta a 2 line file 1)Header, 2) Sequence
.fastq (Illumina), .fasta with ascii codes phred32/phred64 base calling meta-quality-data with spacer lines, a 4 line file 1) Header, 2) Quality, 3) Bases, 4) Spacer Line
.fast5 (Nanopore), .fastq with pore data There are 3 main branches of data stored in the fast5, Analysis, Raw, and UniqueGlobalKey.

.sam Sequence Alinment Map data, human readable sequence alignment file
.bam Binary Alignment Mapping file, gzip and share (loads Ref Genome .fasta with .bams)
