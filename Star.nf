#!/usr/bin/env nextflow

/*
*params input 
*/
params.fastq_1 = "$baseDir/color/*F3.fastq"
params.fastq_2 = "$baseDir/color/*F5-BC.fastq"
params.genome = "$baseDir/color/1M_hg19.fasta"
params.path_Star  ="/usr/local/bin/STAR" 
params.cpu = "4" 
params.index = null
params.help = false

//print usage
if (params.help) {
    log.info ''
    log.info 'STAR pipeline RNA SOLID'
    log.info '-----------------------'
    log.info '.'
    log.info ''
	log.info 'Usage: '
	log.info '    pass.nf --csfasta_1  .../*1.{csfasta,QV.qual} --csfasta_2 .../2.{csfasta,QV.qual} '
	log.info '            --genome GENOME_FILE '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --fastq_1                         First read cs.fasta [ex :".../*1.csfasta."'
    log.info '    --fastq_2                         Second read cs.fasta [ex :".../*2.csfasta"'
    log.info '    --index GENOME_INDEX_FILE           Index file.[Optional but must faster].'
	log.info '    --genome GENOME_FILE 				  Reference genome file(s).'
    log.info '    --path_Star                         path of tool Star (by default in the bin/).'
    log.info '    --cpu                               Numbers of core use (CPU)by default [4].'
	exit 1
}

genome_file = file(params.genome)

read_1 = Channel.fromPath(params.fastq_1).map { file -> tuple(file.baseName[0..-3], file) }
read_2 = Channel.fromPath(params.fastq_2).map { file -> tuple(file.baseName[0..-6], file) }


fastq = read_1.combine(read_2,by:0)



path_Star= params.path_Star
cpu = params.cpu

/*
*Code the F3 color to 0 --> A , 1--> T , 2 --> C , 3 --> G , . --> N
*Code the F5-BC color to 0 --> T , 1--> A , 2 --> G , 3 --> C , . --> N
*/

process codage_fastq{
    tag{id}


    input: 
    set id ,file(read_1), file(read_2) from fastq


    output: 
    set id , file ("*_1.fastq") ,file ("*_2.fastq")into codage_fastq



    """
    sed -e '/^[ATCG]/  s/0/A/g' -e '/^[ATCG]/  s/1/T/g' -e '/^[ATCG]/  s/2/C/g' -e '/^[ATCG]/  s/3/G/g' -e '/^[ATCG]/  s/\\./N/g' ${read_1} > ${id}1.fastq
    sed -e '/^[ATCG]/  s/0/T/g' -e '/^[ATCG]/  s/1/A/g' -e '/^[ATCG]/  s/2/G/g' -e '/^[ATCG]/  s/3/C/g' -e '/^[ATCG]/  s/\\./N/g' ${read_2} > ${id}2.fastq 
    """
    //sed -e '/^[ATCG]/  s/0/A/g' -e '/^[ATCG]/  s/1/T/g' -e '/^[ATCG]/  s/2/C/g' -e '/^[ATCG]/  s/3/G/g' -e '/^[ATCG]/  s/\\./N/g' ${read_2} > ${id}2.fastq
    // on code le read 2 en foward pour l'alignement comme pour ilumina ( RF) solid ( RR ou FF)
}





if(params.index == null){
    
    /*
    *Convert the sequence fasta to color format 
    */
    process conversion_to_csfasta{

        input: 
            file genome from genome_file

        output:
            file "genome_color.csfasta" into genome_ref_color 

        """
        python $baseDir/tool/fasta_to_csfasta_v3.py ${genome} genome_color.csfasta
        """

    }

    /*
    *Code the sequence color to 0 --> A , 1--> T , 2 --> C , 3 --> G , . --> N
    */
    process codage_genome{


        input: 
            file genome from genome_ref_color 

        output: 
            file "genome_ref_codage.fasta" into codage_genome

        """
        sed -e '/^>/! s/0/A/g' -e '/>/! s/1/T/g' -e '/>/! s/2/C/g' -e '/>/! s/3/G/g' -e '/>/! s/\\./N/g' $genome > genome_ref_codage.fasta
        """
    }

    /*
    *Build index with the genome reference convert in color 
    *--genomeDir    name of directory with the index genome
    *--runMode      Mode choose  
    *--runThreadN   number of cores
    *--genomeFastaFiles     Name of input file in fasta 
    *--genomeChrBinNbits  Option if big numbers of reference in file fasta (>...)
    */
    process buildIndex{

        input:
        file genome from codage_genome
    
        output: 
        file "STARgenome" into genomeIndex


    """
        mkdir STARgenome
        ${path_Star} --runThreadN ${cpu} \
             --runMode genomeGenerate \
             --genomeDir STARgenome \
             --genomeFastaFiles $genome \
             --limitGenomeGenerateRAM 12000000000 \

    """
        
    //    --genomeChrBinNbits 12 \ if big number of reference 

    }



    }
else{
    genomeIndex = Channel.fromPath(params.index)
}

//genomeIndex.subscribe { println "$it"}

/*
*Mapping the paired sequence with the genome reference convert 
*OPTION :
*--genomeDir             name of directory with the index genome
*--readFilesIn           name of read file 
*--runThreadN            number of cores
*--outFileNamePrefix     Name of output file 
*--outSAMunmapped        add the unmapped read in BAM/SAM file 
*--outFilterMultimapNmax read alignments will be output only if the read maps fewer than this value,
*                        otherwise no alignments will be output
*--limitBAMsortRAM       maximum available RAM for sorting BAM.
*--clip3pNbases          number(s) of bases to clip from 3p of each mate.  If one value is given, it
*                        will be assumed the same for both mates.
*--chimSegmentMin        minimum length of chimeric segment length
*/


process mapping {
    tag{id}
    publishDir "result/Star"
    input:
    file genome from genomeIndex.first()
    set id , file (read_1) , file (read_2) from codage_fastq


    output:
    file "${id}STAR_mapping" into mappedReads 

    """
    ${path_Star} --runThreadN ${cpu} \
         --genomeDir ${genome} \
         --readFilesIn ${read_1} ${read_2}\
         --outFileNamePrefix $id \
         --outSAMtype BAM SortedByCoordinate\
         --outSAMunmapped Within \
         --clip3pNbases 1 \
         --limitBAMsortRAM 12000000000 \
         --chimSegmentMin 1 \
         --outFilterMultimapNmax 200 \

    mkdir ${id}STAR_mapping
    mv ${id}Aligned* ${id}STAR_mapping/.
    mv ${id}Log* ${id}STAR_mapping/.
    mv ${id}Chimeric* ${id}STAR_mapping/.
    """

}


