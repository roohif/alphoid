#!/bin/bash

blastn -task blastn -dust no -soft_masking false \
	-query X07685.fa -db /home/glenn/genome/T2T/T2T -num_threads 6 -out vs_T2T.csv \
	 -outfmt "10 qseqid qstart qend sseqid sstart send nident pident length qlen evalue"

