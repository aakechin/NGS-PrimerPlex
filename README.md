# NGS-PrimerPlex is a high-throughput tool for mupltiplex primer design
It includes three Python-scripts:
* **getGeneRegions.py** - script that converts list of genes and their parts into genome coordinates for which the next script designs primers.
* **NGS_primerplex.py** - script that designs primers for the user defined genome coordinates.
* **addSeqToPrimers.py** - script that adds adapter sequences for all of the designed primers.
## Docker use (RECOMMENDED!)
NGS-PrimerPlex can be run as a [Docker](https://www.docker.com) image. In this way you only need to [install Docker](https://docs.docker.com/install/) (for windows 7 users [this install steps](https://docs.docker.com/toolbox/toolbox_install_windows/) should be performed). If you have "VD-x, VD-t error", you need to turn on virtualization in BIOS CPU section.

And then download image of NGS-PrimerPlex that contains all necessary modules and reference files:

`docker pull aakechin/ngs-primerplex`

After that, run downloaded image:

`docker run -it --entrypoint 'bash' aakechin/ngs-primerplex`

You will be in the main image directory, where you can find folders with all necessary scripts and reference genome files. Now, you need to unzip reference genome FASTA-file and you will be able to start example primer design or upload to the container (that you've obtained from the image) your list of genes or genome regions:

```
cd NGS-PrimerPlex/ 
gunzip hg19/ucsc.hg19*
python3 getGeneRegions.py -email info@gmail.com -glf example_gene_list_file.txt -ref hg19/ -rf example_gene_list.regions.csv 
python3 NGS_primerplex.py -regions example_gene_list_file.regions.csv -wgref hg19/ucsc.hg19.fasta -minampllen 140 -optampllen 150 -maxampllen 150 -minprimerlen 16 -optprimerlen 25 -maxprimerlen 32 -minprimermelt 60 -optprimermelt 64 -maxprimermelt 68 -minprimergc 23 -maxprimergc 75 -minprimerendgc 0 -maxprimerendgc 4 -maxprimerpolyn 7 -subst 1 -th 8 -run example 
```
This will give you primers that could be designed with the defined parameters. Then, you can use generated file with draft primers as -draft argument and defining less strict parameters for primer design.

For uploading files to a container, you can use the following commands in another terminal (**open new docker terminal**):

`docker ps -a`

to know, what is your container's name. It will be written in the last column.

`docker cp <file or directory that you want to copy to the container> <container name>:/<folder where you want to put your files>`

And then you can run analysis with command like above.

## Tips
1. If you have problems with use of _**vim on Windows docker**_, type the command `:set term=cygwin` in the vim, and it will work fine.
1. If you obtain error _**segmentation fault**_ in the docker container, increase memory and processor performance provided to your virtual machine in the Virtual Box settings.
### Primers could not be designed with the defined parameters
Thanks to the function of draft primers you can subsequently design primers with less and less stringent parameters. Below, the most frequently parameters that you need to change, are listed.
1. Look at sequences that NGS-PrimerPlex outputed for regions for which primers could not be designed. In the most cases, you need to change GC-content of your primers. Leave optimal primer GC as you want, but change minimal and maximal GC-content.
1. Change minimal and maximal GC-content of the primers' ends (5 last nucleotides). For the most complex regions they can be set to 0 and 5, respectively.
1. Change maximal length of primer poly-N (-maxprimerpolyn). It can be enlarged to 7-8 nucleotides.
1. If primers could not be joined into amplified blocks, increase -maxoverlap to allow overlapping of the neiborhing amplicons by their studied regions (between left and right primers). Also, if you have more time, you can try to increase -primernum1 parameter that will lead to designing more primers for each studied position.
1. IMPORTANT! If you have draft-file with primers that were previously filtered by specificity and covering SNPs and you want to join them into multiplexes, you NEED to use -blast option in order to check primers for forming non-specific amplicons between primers from different pairs!

## Requirements
All of the Python-scripts listed above work under Python3+ and require the following additional Python-modules:
* biopython
* argparse
* primer3-py (as program for choosing primer pairs for one region, NGS-PrimerPlex uses primer3-py Python package)
* logging
* pysam
* xlrd
* xlsxwriter
* networkx (version==1.11, newer versions has different sintaxis)

They can be installed with pip:

`sudo pip3 install biopython argparse primer3-py logging pysam xlrd xlsxwriter "networkx==1.11"`

Also, for searching non-target primer hybridization, it uses BWA, so you need to install it with e.g.:

`sudo apt-get install bwa`

Then, download reference genome sequence (e.g. [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) or [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) human genome version), convert it to FASTA-file with [twoBitToFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa) (if it is not in this format) and index it with BWA:

```
twoBitToFa hg19.2bit ucsc.hg19.fa
bwa index ucsc.hg19.fasta
```

It will take some time. If you want to automatically extract genome regions for genes needed, you will have to also download GenBank-files for each of chromosome for genome version that you are going to use, e.g. from [NCBI Genome database](https://www.ncbi.nlm.nih.gov/genome/). 

If you want to check primers for crossing SNPs, download dbSNP VCF-file for the correspondent version of human genome. For hg19:
```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi
```
For hg38:
```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi
```

And, finally, you can run your primer design:

```
cd NGS-PrimerPlex/ 
python3 getGeneRegions.py -email info@gmail.com -glf example_gene_list_file.txt -ref hg19/ -rf example_gene_list.regions.csv 
python3 NGS_primerplex.py -regions example_gene_list_file.regions.csv -wgref hg19/ucsc.hg19.fasta -minampllen 140 -optampllen 150 -maxampllen 150 -minprimerlen 16 -optprimerlen 25 -maxprimerlen 32 -minprimermelt 60 -optprimermelt 64 -maxprimermelt 68 -minprimergc 23 -maxprimergc 75 -minprimerendgc 0 -maxprimerendgc 4 -maxprimerpolyn 7 -subst 1 -th 8 -run example 
```

Unfortunately, native use of NGS-PrimerPlex is possible only undex Linux or Mac OS. For Windows users, we suggest to use docker version that we recommend also for Linux and Mac OS users, because it already includes all pre-installed packages and files required.

Below three scripts of NGS-PrimerPlex are described in details.
## getGeneRegions.py
This script takes names of genes and numbers of their exons or positions in CDS and makes region file for NGS_primerplex.py. It has the following arguments:
```
  -h, --help            show this help message and exit
  --geneListFile GENELISTFILE, -glf GENELISTFILE
                        file with list of genes. Format is: GENE EXONS CODONS
  --refDir REFDIR, -ref REFDIR
                        directory with reference files
  --organism ORGANISM, -org ORGANISM
                        common name of organism (human, sheep, etc). Default:
                        human
  --resultFile RESULTFILE, -rf RESULTFILE
                        file for results
  --intron-nucleotides INTRONSIZE, -intron INTRONSIZE
                        number of nucleotides from intron to take. Default: 2
  --include-noncoding, -noncoding
                        use this parameter, if you want to include 5'- and
                        3'-non-coding regions of mRNA
```
An example gene list file is included into the repository (example_gene_list_file.txt). User has the following opportunities to define regions to be studied:
* Whole gene - then only gene name should be written in the first column.
* One or list of exons of some gene - then user should write name of gene into the 1st column and list of exons into the 2nd column seperated by comma (without spaces). E.g.: 1,2,3.
* Range of exons of some gene - then user should write name of gene into the 1st column and range from one exon to another separated by "-", e.g.: 1-5.
* Mix of two previous variants, e.g.: 1,4,5-10,15.
* One or list of codons of some gene - then user should write name of gene into the 1st column and list of codons into the 3rd column seperated by comma (without spaces). E.g.: 100,202,303.
* Range of codons of some gene - then user should write name of gene into the 1st column and range from one codon to another separated by "-", e.g.: 1-5.
* Mix of two previous variants, e.g.: 1,4,5-10,15.
* Mix of defined exons and codons.
```
Directory with reference files means that in some directory genbank-files (GB-files) for all of the chromosomes of the reference genome should be located. The program reads these GB-files and determines coordinates of genes, their exons, introns, and codons. 
Number of nucleotides from intron to take means that by default NGS-primerplex extracts only exon coordinates and two nucleotides from neighbouring introns.
Argument -noncoding is necessary for including also non-coding exons when user defines only name of gene to study.
```
## NGS_primerplex.py
This is the main script of this tool. It takes list of genome regions for which user needs to design primers. It has the following format:

| chromosome | region start | region end | amplicon name | desired multiplex numbers (optional) | type of primers (left, right or both, optional) | use this region as one amplicon (optional) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1000000 | 1000100 | RANDOM_REGION | 1,2,3 | B | W |

In this file user can manually define the following features of primer design:
* _**desired multiplex numbers**_ - list of primer pools to which user wants to put amplicons for this region. For exampple, if all amplicons should be splited onto 4 multiplex reactions, write "1,2,3,4" to this column for each region. If some regions should be put into two multiplex reactions and other into two other ones, write "1,2" for first regions and "3,4" for other ones. If you don't need splitting onto multiplex reactions, leave this column empty.
* _**type of primers (left, right or both)**_. This can be useful, if someone creates amplicon-based library for detection of gene fusions and primers only from one side are needed. If you need both primers, you can leave it empty. In other cases, use "L", "R", "B".
* sometimes, especially, when long deletions are detected, you need to amplify _**whole region as one amplicon**_, e.g. EGFR deletion in the exon 19. For such regions, you can specify in this column "W", and NGS-PrimerPlex will try to design primers for whole such region. In other cases leave it empty, and NGS-PrimerPlex will design primers for each position of this region distinctly and choose the best variant. Then it can cover whole region by one as well as several primer pairs.
### Other features of NGS-PrimerPlex can be used by defining parameters:
1. _**Multi-step primer design**_ for NGS-panels with regions that have structure difficult for primer design. For example, you have regions with low GC-content. If you design primers in one step, you need to allow low GC-content for all primers designed. But if you do it in several steps, you can initially start designing primers with strict parameters. The program will construct primers for regions with normal GC-content and save them into "draft primers" file. For regions with too low or too high GC-content (or other parameters like Tm, poly-N length, number of non-target hybridizations etc.) it will raise an error. After that you can use parameter _-draft_ with choosing file that was generated by NGS-PrimerPlex (it ends with "all_draft_primers.xls") and choosing less strict parameters. Then the program will not design primers for regions for which primers has been constructed earlier. This will significantly speed up choosing acceptable primer design parameters.
    - Similar approach can be used when some region has repeats in genome or contain high-frequent SNPs. In this case, you can use "draft primers" that have already passed filtering by genome mapping and/or checking for covering SNPs. This files end with "all_draft_primers_after_specificity.xls" and "all_draft_primers_after_SNPs.xls". The last one will have primers that passed also through mapping genome, if this was chosen. Remember, that for distributing primer pairs among multiplexes, you still need -blast option to check for forming non-specific amplicons between primers from different primer pairs!
2. _**Checking primers for non-target hybridization**_. The program will check all of designed primers for number of non-target hybridizations in genome with defined number of substitutions/insertions/deletions. Primer pairs can be filtered out by number of such non-target sites for one of the primers (e.g. to exclude primers that are complement to genome repeats) and by forming non-target amplicons. Forming of non-target amplicons is checked both for each primer pairs designed and for primers from different pairs while combining them into multiplex reactions. For performing such analysis, use parameters -blast and -subst. The last one controls number of substitutions, insertions and deletions allowed while searching for primer sequences in genome.
3. _**Checking primers for covering high-frequent SNPs**_. The program will check each primer for covering variable sites of genome by 3'-end nucleotide. This is necessary for decreasing influence of sample SNPs onto amplification efficiency. For performing such analysis, use parameters -snps and -gv (--genome-version).
4. _**Designing primers for embedded PCR**_. The program can design both internal and external primers for one run or design only external primers for previously designed internal primers. Use of external primers can increase sensitivity of amplicon-based NGS-panel. For doing it, use parameters -embedded, -minprimershift, -optextampllen, and -maxextampllen.

Other parameters of NGS-primerplex.py are listed below (the most of parameters have default values):
```
  -h, --help            show this help message and exit
  --regions-file REGIONSFILE, -regions REGIONSFILE
                        file with regions for amplification in the following f
                        ormat:Chromosome{Tab}Start_Position{Tab}End_Position{T
                        ab}Amplicon_Name{Tab} Desired_Multiplex_Numbers(option
                        al){Tab}Type_Of_Primers(only left/only
                        right/both)(optional){Tab}Use_Whole_Region(optional)
  --primers-file PRIMERSFILE, -primers PRIMERSFILE
                        file with previously designed internal primers. Use
                        this parameter, if you want only to design external
                        primers
  --draft-primers DRAFTFILE, -draft DRAFTFILE
                        file with internal primers previously designed for
                        part of input regions. The program will design primers
                        for the left regions
  --whole-genome-ref WHOLEGENOMEREF, -wgref WHOLEGENOMEREF
                        file with INDEXED whole-genome reference sequence
  --min-amplicon-length MINAMPLLEN, -minampllen MINAMPLLEN
                        minimal length of amplicons. Default: 75
  --max-amplicon-length MAXAMPLLEN, -maxampllen MAXAMPLLEN
                        maximal length of amplicons. Default: 100
  --optimal-amplicon-length OPTAMPLLEN, -optampllen OPTAMPLLEN
                        optimal length of amplicons. Default: 90
  --min-primer-length MINPRIMERLEN, -minprimerlen MINPRIMERLEN
                        minimal length of primers. Default: 18
  --max-primer-length MAXPRIMERLEN, -maxprimerlen MAXPRIMERLEN
                        maximal length of primers. Default: 25
  --optimal-primer-length OPTPRIMERLEN, -optprimerlen OPTPRIMERLEN
                        optimal length of primers. Default: 23
  --min-primer-melting-temp MINPRIMERMELT, -minprimermelt MINPRIMERMELT
                        minimal melting temperature of primers, degrees
                        Celsius. Default: 62
  --max-primer-melting-temp MAXPRIMERMELT, -maxprimermelt MAXPRIMERMELT
                        maximal melting temperature of primers, degrees
                        Celsius. Default: 66
  --optimal-primer-melting-temp OPTPRIMERMELT, -optprimermelt OPTPRIMERMELT
                        optimal melting temperature of primers, degrees
                        Celsius. Default: 64
  --min-primer-gc MINPRIMERGC, -minprimergc MINPRIMERGC
                        minimal acceptable GC-content for primers. Default: 25
  --max-primer-gc MAXPRIMERGC, -maxprimergc MAXPRIMERGC
                        maximal acceptable GC-content for primers. Default: 60
  --optimal-primer-gc OPTPRIMERGC, -optprimergc OPTPRIMERGC
                        optimal acceptable GC-content for primers. Default: 40
  --min-primer-end-gc MINPRIMERENDGC, -minprimerendgc MINPRIMERENDGC
                        minimal acceptable number of G or C nucleotides within
                        last 5 nucleotides of 3'-end of primers. Default: 1
  --max-primer-end-gc MAXPRIMERENDGC, -maxprimerendgc MAXPRIMERENDGC
                        maximal acceptable number of G or C nucleotides within
                        last 5 nucleotides of 3'-end of primers. Default: 3
  --opt-primer-end-gc OPTPRIMERENDGC, -optprimerendgc OPTPRIMERENDGC
                        optimal number of G or C nucleotides within last 5
                        nucleotides of 3'-end of primers. Default: 2
  --max-primer-poly-n MAXPRIMERPOLYN, -maxprimerpolyn MAXPRIMERPOLYN
                        maximal acceptable length of some poly-N in primers.
                        Default: 3
  --max-primer-compl-end-th MAXPRIMERCOMPLENDTH, -maxprimercomplendth MAXPRIMERCOMPLENDTH
                        maximal Tm for complementarity of 3'-ends of primers.
                        Default: 15
  --max-primer-compl-any-th MAXPRIMERCOMPLANYTH, -maxprimercomplanyth MAXPRIMERCOMPLANYTH
                        maximal Tm for any complementarity of primers.
                        Default: 30
  --max-primer-hairpin-th MAXPRIMERHAIRPINTH, -maxprimerhairpinth MAXPRIMERHAIRPINTH
                        maximal melting temperature of primer hairpin
                        structure. Default: 40
  --max-primer-nonspecific MAXPRIMERNONSPEC, -maxprimernonspec MAXPRIMERNONSPEC
                        maximal number of nonspecific regions to which primer
                        can hybridizes. Default: 1000
  --max-amplicons-overlap MAXOVERLAP, -maxoverlap MAXOVERLAP
                        maximal length of overlap between two amplified blocks
                        (it does not include primers). Default: 5
  --primers-number1 PRIMERNUM1, -primernum1 PRIMERNUM1
                        number of primer that user wants to get on the 1st
                        stage. The more this value, the more precise the
                        choice of primers, but the longer the design time.
                        Default: 5
  --auto-adjust-parameters, -autoadjust
                        use this parameter if you want NGS-PrimerPlex to
                        automatically use less stringent parameters if no
                        primer were constructed for some region
  --tries-to-get-best-combination TRIESTOGETCOMBINATION, -tries TRIESTOGETCOMBINA
                        number of of tries to get the best primer combination.
                        More the value, better combination will be, but this
                        will take more time. Default: 1000
  --return-variants-number RETURNVARIANTSNUM, -returnvariantsnum RETURNVARIANTSNUM
                        number of multiplexes variants that user wants to get
                        after all analyses and filters. Default: 1
  --embedded-amplification, -embedded
                        use this parameter if you want to create NGS-panel
                        with embedded amplification
  --min-internal-primer-shift MINPRIMERSHIFT, -minprimershift MINPRIMERSHIFT
                        minimal shift of external primer from the 3'-end of
                        internal primer. Default: 5
  --opt-external-amplicon-length OPTEXTAMPLLEN, -optextampllen OPTEXTAMPLLEN
                        optimal length of the external amplicons. Default: 110
  --max-external-amplicon-length MAXEXTAMPLLEN, -maxextampllen MAXEXTAMPLLEN
                        maximal length of the external amplicons. Default: 130
  --do-blast, -blast    use this parameter if you want to perform Blast-
                        analysis of constructed primers
  --substititutions-num SUBSTNUM, -subst SUBSTNUM
                        accepted number of substitutions for searching primers
                        in genome. Default: 2
  --max-nonspecific-amplicon-length MAXNONSPECLEN, -maxnonspeclen MAXNONSPECLEN
                        maximal length of nonspecific amplicons that the
                        program should consider. For example, if you design
                        primers for DNA from serum, you can set it as 150.
                        Default: 200
  --snps, -snps         use this parameter if you want to check that 3'-ends
                        of your primers do not cover any SNPs with high
                        frequency
  --dbsnp-vcf DBSNPVCFFILE, -dbsnp DBSNPVCFFILE
                        VCF-file (may be gzipped) with dbSNP variations
  --snp-freq SNPFREQ, -freq SNPFREQ
                        minimal frequency of SNP in whole population to
                        consider it high-frequent SNP. Default: 0.05
  --nucletide-number-to-check NUCNUMTOCHECK, -nucs NUCNUMTOCHECK
                        Number of nucleotides from 3`-end to check for
                        covering SNPs. Default: None and the program will
                        check all nucleotides
  --min-multiplex-dimer-dg1 MINMULTDIMERDG1, -minmultdimerdg1 MINMULTDIMERDG1
                        minimal acceptable value of free energy of primer
                        dimer formation with hybridized 3'-end in one
                        multiplex in kcal/mol. Default: -6
  --min-multiplex-dimer-dg2 MINMULTDIMERDG2, -minmultdimerdg2 MINMULTDIMERDG2
                        minimal acceptable value of free energy of primer
                        dimer formation in one multiplex in kcal/mol. Default:
                        -10
  --threads THREADS, -th THREADS
                        number of threads. Default: 2
  --run-name RUNNAME, -run RUNNAME
                        name of program run. It will be used in the output
                        file names
  --skip-uncovered, -skip
                        use this parameter if you want to skip some targets
                        for which primers can not be designed with defined
                        parameters
  --monovalent-concentration MVCONC, -mv MVCONC
                        Concentration of monovalent cations, commonly K+ or
                        NH4+, in mM. Default: 50
  --divalent-concentration DVCONC, -dv DVCONC
                        Concentration of divalent cations, commonly Mg2+, in
                        mM. Default: 3
  --dntp-concentration DNTPCONC, -dntp DNTPCONC
                        Total concentration of dNTPs. If you have each dNTP
                        with concantration 0.2 mM, then total is 0.8 mM.
                        Default: 0.8
  --primer-concentration PRIMERCONC, -primerconc PRIMERCONC
                        Concentration of each primer, in nM. Default: 250
```
## addSeqToPrimers.py
This script adds adapter sequences to all designed primers. As an input it uses NGS-primerplex.py output file and file with adapter sequences. An example file with adapter sequences is included into the repository. This script outputs sequences into new XLS-file listing all designed primers with names and adapter sequences added.
All parameters are listed below:
```
  -h, --help            show this help message and exit
  --input INPUT, -in INPUT
                        input XLS-file with designed primers
  --tags-file TAGSFILE, -tags TAGSFILE
                        text file with tags that we want to add to each
                        primer. Default: "/NGS-
                        PrimerPlex/kplex_for_primers.txt"
```
## Citation
Manuscript is preparing. Now you can cite to this repository.
