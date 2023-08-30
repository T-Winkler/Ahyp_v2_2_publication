#!/bin/bash

# Used to remove ambiguous bases from the reference genome to prepare as input for NextPolish

# do not remove Ns but write them into a new entry instead, creating N specific chromosomes. All entries and their order is saved
# N specific entries can in a last step be removed. To reconstruct, the saved order can be used to integrate the N chromosomes again (which record the number of Ns removed)

# Change these two parameters
INPUT=reference_genomes/Ahypochondriacus/assembly/Ahypochondriacus_459_v2.0.nospace.underscore.fa
OUTDIR=data/NextPolish/input/

mkdir -p $OUTDIR

# split fasta into new entrys based on gap character N
# print into a single line; replace stretch of Ns with new "split" entry; in the end create newline after first fasta entry
LC_ALL=C awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' $INPUT | sed 's/N*N/\n>spl\n&\n>spl\n/g' > "$OUTDIR"tmp.txt

echo "splits generated"

# start with 0
i=0

# for each line, if it is a newly created header, increase value of i by 1 and add _i to the header name
for j in $(cat "$OUTDIR"tmp.txt); do
	if [[ $j =~ .sp* ]] ; then
                i=$((i+1))
        fi
	echo $j | sed "s/>spl/>spl_"$i"/"
done > "$OUTDIR"tmp2.txt

echo "In total: "$i" splits renamed"


# save the correct order of all headers in a file to be able to restore the order later
grep ">" "$OUTDIR"tmp2.txt > "$OUTDIR"out.headers.txt

echo "headers saved"

# remove headers without sequence (Can be caused by stretch of Ns at the start of a Scaffold (see Scaffold 10))
sed -r 'N; /(>)[^\n]*\n\1/ s/[^\n]*//; P; D' "$OUTDIR"tmp2.txt | grep . | grep -i -B 1 --no-group-separator  '[ATGC]'  > "$OUTDIR"data/NextPolish/input/Ahypochondriacus_split.fasta

# remove and rename temporary files
rm "$OUTDIR"tmp.txt
mv "$OUTDIR"tmp2.txt "$OUTDIR"out.prefiltered.txt

