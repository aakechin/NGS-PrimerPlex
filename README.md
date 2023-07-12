# NGS-PrimerPlex is a high-throughput tool for mupltiplex primer design
It includes four Python-scripts:
* **getGeneRegions.py** - script that converts list of genes and their parts into genome coordinates for which the next script designs primers.
* **NGS_primerplex.py** - script that designs primers for the user defined genome coordinates.
* **addSeqToPrimers.py** - script that adds adapter sequences for all of the designed primers.
* **convertToDraftFile.py** - script that converts the main NGS-PrimerPlex output file into the draft-primers file.
## Docker use (RECOMMENDED!)
NGS-PrimerPlex can be run as a [Docker](https://www.docker.com) image. In this way you only need to [install Docker](https://docs.docker.com/install/) (for windows 7 users [this install steps](https://docs.docker.com/toolbox/toolbox_install_windows/) should be performed). If you have "VD-x, VD-t error", you need to turn on virtualization in BIOS CPU section.

Now users have two options of NGS-PrimerPlex use in docker: (1) with already uploaded human reference genome hg19 version, (2) without any reference genomes. The 1st variant is idead for use with hg19 genome, but you will have to download about 8 Gb of data. The 2nd is ideal for use with other reference genomes, including other organisms. In this case you will have to download about only 0.5 Gb, but you also have to download reference genome files manually (see section "Reference genome for other organisms than human").

To use 1st variant, download docker image of NGS-PrimerPlex with the following command:

`docker pull aakechin/ngs-primerplex:full_1.3.4`

To use the 2nd variant, download docker image of NGS-PrimerPlex in the following way:

`docker pull aakechin/ngs-primerplex:1.3.4`

__Windows users__ will also have to change some default settings of the Virtual Machines. For Windows 7 it can be done in the Oracle VM VirtualBox, for Windows 10 users in the Docker Settings:
1. Go to settings of the default virtual machine (VM) -> Shared folders -> Add your hardware drives that will be used in docker, giving them names like '/C/' for disk 'C:\' and turning on options 'Auto-mount' and 'Make Permanent' (if the last one is available).
1. Go to 'System' settings -> change 'Base Memory' value to at least 3 Gb (more the better).
1. Switch to the 'Processor' tab -> Change 'Processors' value to at least 3 (more the better) -> Change 'Execution Cap' value to 100%.
Save new settings (click OK), turn off your VM (if it was turned on previously) and restart docker with runnning 'Docker Quickstart Terminal'

At this step, users also have two options of NGS-PrimerPlex use: in the command-line and with GUI.

### Command-line version
If you downloaded version with previously uploaded hg19 reference (aakechin/ngs-primerplex:latest), you will have to gunzip reference genome FASTA-file:

`docker run -it --entrypoint 'bash' --name ngs_primerplex_ref -v '<directory where you are going to design new primers>:<name of this directory in the container>' aakechin/ngs-primerplex:latest`, where -v option lets you to mount some of your local directory to the virtual machine (container). This command will put you into the virtual machine command line. Note, that Windows users can only mount folders from drives that were shared and they should be written as '/C/...'

`gunzip NGS-PrimerPlex/hg19/ucsc.hg19.fasta*.gz`

The last command will take some time. After that, you can run testing of NGS-PrimerPlex (for version without uploaded reference genome you initially need to prepare your reference genome, see **"Reference genome"**):

`python3 /NGS-PrimerPlex/test.py`

All of the tests should be completed successfully. If you met any errors, report about it in the Issues at the GitHub here, please.

Now, you will be able to start example primer design or your own list of genes from folder that was mounted to the container (with -v version, and also in shared folders for Windows users): 

```
cd /NGS-PrimerPlex
python3 getGeneRegions.py -glf example_gene_list_file.txt -ref hg19/ -rf example_gene_list_file.regions.csv 
python3 NGS_primerplex.py -regions example_gene_list_file.regions.csv -ref hg19/ucsc.hg19.fasta -blast -snps -dbsnp hg19/common_all_20180423_hg19.vcf.gz
```
This will give you primers that could be designed with the default parameters. The default parameters are defined in such a way that a user can surely obtain designed primers for the example. For a subsequent use of the program, we recommend to use more stringent parameters. Then, you can use generated file with draft primers as -draft argument and defining less strict parameters for primer design.

### GUI-version
To use GUI-version of NGS-PrimerPlex you need to download also NGS-PrimerPlex from GitHub and install Python and some additional Python modules. To install all of it automatically, NGS-PrimerPlex main package (from GitHub) contains two scrips:
* **install_GUI_linux.sh** - for Linux users
* **install_GUI_windows.bat** - for Windows users. But Windows users will have to install Python3+ manually, downloading it from python.org. During installation, remember to switch "Add Python to PATH".

Run script that is dedicated for your case. After that you can run GUI-version from the command line from NGS-PrimerPlex folder:
* `python main.py` - for Windows users
* `python3 main.py` - for Linux users

If you downloaded version with reference hg19 genome, press 'Prepare hg19 reference' and wait until this button become disabled. In the GUI-version you can choose files and run all steps maximally intuitive.

## Requirements
Non-docker version is available only for Linux and iOS users. To install automatically all of the requirements, run the following commands:
```
chmod +x install_for_linux.sh
./install_for_linux.sh
```
Also additional Python-modules can be installed manually:
* biopython
* argparse
* primer3-py (as program for choosing primer pairs for one region, NGS-PrimerPlex uses primer3-py Python package)
* pysam
* xlrd
* xlsxwriter
* networkx (version==1.11, newer versions have different sintaxis)
* numpy

They can be installed with pip:

`sudo pip3 install biopython argparse primer3-py pysam xlrd xlsxwriter "networkx==1.11" numpy`

Also, for searching non-target primer hybridization, it uses BWA, so you will also need to install it manually with e.g.:

`sudo apt-get install bwa`

### Reference genome
For genome other than human, go to the next Chapter **"Reference genome for other organisms than human"**.

If you use non-docker version of NGS-PrimerPlex or docker-version without uploaded hg19 reference genome, download it (e.g. [hg19](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit) or [hg38](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit) human genome version), convert it to FASTA-file with [twoBitToFa](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa) (if it is not in this format) and index it with BWA:

```
twoBitToFa hg19.2bit ucsc.hg19.fa
bwa index ucsc.hg19.fasta
```

It will take some time. If you want to automatically extract genome regions for genes needed, you will have to also download GenBank-files for each of chromosome for genome version that you are going to use, e.g. from [NCBI Genome database](https://www.ncbi.nlm.nih.gov/genome/). Each GenBank-file should be named as this chromosome is called in the reference genome FASTA-file or as it is ordered in the reference FASTA-file. For example, for the above hg19 version chromosome 1 GenBank-file can be named as chr1.gb or 2.gb (because in the reference genome chrM is written as the 1st chromosome and chr1 as the 2nd).

If you want to check primers for crossing SNPs, download dbSNP VCF-file for the correspondent version of human genome. For hg19 (by default, it is already downloaded to the docker image):
```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz.tbi
```
For hg38:
```
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz.tbi
```
#### Reference genome for other organisms than human
To prepare your own reference genome, you need to prepare one FASTA-file with whole reference genome and one directory (it can be the same as for FASTA-file) with GenBank-files for each of an organism chromosome. For example, to prepare reference genome for Arabidopsis thaliana, download reference genome FASTA-file from [Genome database of NCBI](https://www.ncbi.nlm.nih.gov/genome/?term=arabidopsis+thaliana): click "genome" in the line "Download sequences in FASTA format for __genome__". Extract downloaded archive.

Then, download each chromosome of A. thaliana in GenBank format. To do it, go to the bottom of the Genome database page for A. thaliana and click on each chromosome in format like "NC_003070.9". On the opened page (https://www.ncbi.nlm.nih.gov/nuccore/NC_003070.9 should be opened) at the right menu "Customize view" choose "Customize" and then "Show sequence". Press "Update view". At the top of the page click "Send to:" -> File -> "Create File". Each GenBank-file should be named as this chromosome is called in the reference genome FASTA-file or as it is ordered in the reference FASTA-file. For example, for the A. thaliana chromosome 1 GenBank-file can be named as NC_003070.9.gb or as 1.gb; chloroplast genome as NC_000932.1.gb or 7.gb.

And, finally, you can run your primer design as it is written for your variant of NGS-PrimerPlex (see above).

## Tips
1. If during redistribution of primer pairs among multiplex reactions, only 1-6 pairs could not be redistributed, try to run designing again with -draft option and increasing -returnvatiantsnum option. Because, choosing primer pair combinations is quite stochastic, sometimes several restarts of design process can give good results.
1. If you have problems with use of _**vim on Windows docker**_, type the command `:set term=cygwin` in the vim, and it will work fine.
1. If you obtain error _**segmentation fault**_ in the docker container, increase memory and processor performance provided to your virtual machine in the Virtual Box settings.
1. If some of the designed primers don't work or work with low efficiency, you can remove them from the output NGS-PrimerPlex file (\*info.xls) and begin designing again with -draft option.
### Primers could not be designed with the defined parameters
Thanks to the function of draft primers you can subsequently design primers with less and less stringent parameters. Below, the most frequently parameters that you need to change, are listed.
1. Look at sequences that NGS-PrimerPlex outputed for regions for which primers could not be designed. In the most cases, you need to change GC-content of your primers. Leave optimal primer GC as you want, but change minimal and maximal GC-content.
1. Change minimal and maximal GC-content of the primers' ends (5 last nucleotides). For the most complex regions they can be set to 0 and 5, respectively.
1. Change maximal length of primer poly-N (-maxprimerpolyn). It can be enlarged to 7-8 nucleotides.
1. If primers could not be joined into amplified blocks, increase -maxoverlap to allow overlapping of the neiborhing amplicons by their studied regions (between left and right primers). Also, if you have more time, you can try to increase -primernum1 parameter that will lead to designing more primers for each studied position.
1. IMPORTANT! If you have draft-file with primers that were previously filtered by specificity and covering SNPs and you want to join them into multiplexes, you NEED to use -blast option in order to check primers for forming non-specific amplicons between primers from different pairs!
1. If you gonna to design primers for embedded PCR, use higher number of returned variants (>=10), because not for all internal primers external primers can be designed.
1. If you want to change distribution of your already designed primers into multiplex reactions, you can make draft XLS-file from the result (*_info.xls*) file and run NGS-PrimerPlex again. Or, if some of primer pairs couldn't be sorted to any of the multiplexes, you can remove such primer pairs from the "*_info.xls*" file, convert it to draft-file and run NGS-PrimerPlex to design new primers for this regions and to resort all primers into new multiplex reactions.

## Description of all scripts
Below three scripts of NGS-PrimerPlex are described in details.
### getGeneRegions.py
This script takes names of genes and numbers of their exons or positions in CDS and makes regions-file for NGS_primerplex.py. It has the following arguments:
```
  -h, --help            show this help message and exit
  --geneListFile GENELISTFILE, -glf GENELISTFILE
                        file with list of genes. Format is: GENE EXONS CODONS
  --refDir REFDIR, -ref REFDIR
                        directory with reference files
  --reference-genome WHOLEGENOMEREF, -wgref WHOLEGENOMEREF
                        file with INDEXED whole-genome reference sequence
  --resultFile RESULTFILE, -rf RESULTFILE
                        file for results
  --intron-nucleotides INTRONSIZE, -intron INTRONSIZE
                        number of nucleotides from intron to take. Default: 2
  --include-noncoding, -noncoding
                        use this parameter, if you want to include 5'- and
                        3'-non-coding regions of mRNA
```
Two example gene list files are included into the repository (example_gene_list_file.txt and example_gene_list_file2.txt). For each of them the repository also contains output files for this script: for hg19 and hg38 versions of human genome. Note, that for EGFR default value "NotW" is manually replaced with "W", because this is extended deletion of 15 nucleotides in the exon 19. User has the following opportunities to define regions to be studied:
* Whole gene - then only gene name should be written in the first column.
* One or list of exons of some gene - then user should write name of gene into the 1st column and list of exons into the 2nd column seperated by comma (without spaces). E.g.: 1,2,3.
* Range of exons of some gene - then user should write name of gene into the 1st column and range from one exon to another separated by "-", e.g.: 1-5.
* Mix of two previous variants, e.g.: 1,4,5-10,15.
* One or list of codons of some gene - then user should write name of gene into the 1st column and list of codons into the 3rd column seperated by comma (without spaces). E.g.: 100,202,303.
* Range of codons of some gene - then user should write name of gene into the 1st column and range from one codon to another separated by "-", e.g.: 1-5.
* Mix of two previous variants, e.g.: 1,4,5-10,15.
* Mix of defined exons and codons.

Directory with reference files means that in some directory genbank-files (GB-files) for all of the chromosomes of the reference genome should be located. The program reads these GB-files and determines coordinates of genes, their exons, introns, and codons. 
Number of nucleotides from intron to take means that by default NGS-primerplex extracts only exon coordinates and two nucleotides from neighbouring introns.
Argument -noncoding is necessary for including also non-coding exons when user defines only name of gene to study.

### NGS_primerplex.py
This is the main script of this tool. It takes list of genome regions for which user needs to design primers. It has the following format:

| chromosome | region start | region end | amplicon name | desired multiplex numbers (optional) | type of primers (left - L, right -R or both -B, optional) | use this region as one amplicon (optional) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 1000000 | 1000100 | RANDOM_REGION | 1,2,3 | B | W |

In this file user can manually define the following features of primer design:
* _**desired multiplex numbers**_ - list of primer pools to which user wants to put amplicons for this region. For exampple, if all amplicons should be splited onto 4 multiplex reactions, write "1,2,3,4" to this column for each region. If some regions should be put into two multiplex reactions and other into two other ones, write "1,2" for first regions and "3,4" for other ones. If you don't need splitting onto multiplex reactions, leave this column empty. By default, getGeneRegions writes here "1".
* _**type of primers (left, right or both)**_. This can be useful, if someone creates amplicon-based library for detection of gene fusions and primers only from one side are needed. If you need both primers, you can leave it empty. In other cases, use "L", "R", or "B".
* sometimes, especially, when long deletions are detected, you need to amplify _**whole region as one amplicon**_, e.g. EGFR deletion in the exon 19. For such regions, you can specify in this column "W", and NGS-PrimerPlex will try to design primers for whole such region. In other cases leave it empty or write something different (like NotW as it is written by default), and NGS-PrimerPlex will design primers for each position of this region distinctly and choose the best variant. Then it can cover whole region by one as well as several primer pairs.
### Other features of NGS-PrimerPlex can be used by defining parameters:
1. _**Multi-step primer design**_ for NGS-panels with regions that have structure difficult for primer design. For example, you have regions with low GC-content. If you design primers in one step, you need to allow low GC-content for all primers designed. But if you do it in several steps, you can initially start designing primers with strict parameters. The program will construct primers for regions with normal GC-content and save them into "draft primers" file. For regions with too low or too high GC-content (or other parameters like Tm, poly-N length, number of non-target hybridizations etc.) it will raise an error. After that you can use parameter _-draft_ with choosing file that was generated by NGS-PrimerPlex (it ends with "all_draft_primers.xls") and choosing less strict parameters. Then the program will not design primers for regions for which primers have been constructed earlier. This will significantly speed up choosing acceptable primer design parameters.
    - Similar approach can be used when some region has repeats in genome or contain high-frequent SNPs. In this case, you can use "draft primers" that have already passed filtering by genome mapping and/or checking for covering SNPs. This files end with "all_draft_primers_after_specificity.xls" and "all_draft_primers_after_SNPs.xls", respectively. The last one will have primers that passed also through mapping genome, if this was chosen. Remember, that for distributing primer pairs among multiplexes, you still need -blast option to check for forming non-specific amplicons between primers from different primer pairs!
    - Another example of using draft-file is a redistribution of primer pairs among multiplex reactions if some primers couldn't be sorted. In this case, you can remove such primer pairs from the *_info.xls* file, convert it to draft file and run NGS-PrimerPlex again, with -blast option, too.
2. _**Checking primers for non-target hybridization**_. The program will check all of designed primers for number of non-target hybridizations in genome with defined number of substitutions/insertions/deletions. Primer pairs can be filtered out by number of such non-target sites for one of the primers (e.g. to exclude primers that are complement to genome repeats) and by forming non-target amplicons. Forming of non-target amplicons is checked both for each primer pairs designed and for primers from different pairs while combining them into multiplex reactions. For performing such analysis, use parameters -blast and -subst. The last one controls number of substitutions, insertions and deletions allowed while searching for primer sequences in genome.
3. _**Checking primers for covering high-frequent SNPs**_. The program will check each primer for covering variable sites of genome by whole sequence or only by 3'-end. Number of nucleotides from 3'-end is controlled by the -nucs parameter. If you want to check whole primer, simply do not use this argument. This check is necessary for decreasing influence of sample SNPs onto amplification efficiency. For performing such analysis, use parameters -snps and -dbsnp with defining path to the VCF-file of dbSNP database (see detailed information above).
4. _**Designing primers for embedded PCR**_. The program can design both internal and external primers for one run or design only external primers for previously designed internal primers. Use of external primers can increase sensitivity of amplicon-based NGS-panel. For doing it, use parameters -embedded, -minprimershift, -optextampllen, and -maxextampllen.
5. _**Designing primers for anchored PCR**_. Anchored PCR allows detecting gene fusions in mRNA samples. After mRNA fragmentation, reversion of it, synthesis of the second strand for cDNA, and adapter ligation, one can amplify exon juntions with two primers: one gene-specific primer and one primer for adapter sequence. In this case, user can detect different gene fusions, including previously unknown. You can try it with the example file _example_gene_list_file2.txt_.

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
  --reference-genome WHOLEGENOMEREF, -ref WHOLEGENOMEREF
                        file with INDEXED whole-genome reference sequence
  --adapter-for-left LEFTADAPTER, -ad1 LEFTADAPTER
                        adapter for left primers. Use it, if you want to
                        preserve formation of second structures with adapter
                        sequences (optional)
  --adapter-for-right RIGHTADAPTER, -ad2 RIGHTADAPTER
                        adapter for right primers. Use it, if you want to
                        preserve formation of second structures with adapter
                        sequences (optional)
  --min-amplicon-length MINAMPLLEN, -minampllen MINAMPLLEN
                        minimal length of amplicons. Default: 100
  --max-amplicon-length MAXAMPLLEN, -maxampllen MAXAMPLLEN
                        maximal length of amplicons. Default: 110
  --optimal-amplicon-length OPTAMPLLEN, -optampllen OPTAMPLLEN
                        optimal length of amplicons. Default: 110
  --min-primer-length MINPRIMERLEN, -minprimerlen MINPRIMERLEN
                        minimal length of primers. Default: 16
  --max-primer-length MAXPRIMERLEN, -maxprimerlen MAXPRIMERLEN
                        maximal length of primers. Default: 28
  --optimal-primer-length OPTPRIMERLEN, -optprimerlen OPTPRIMERLEN
                        optimal length of primers. Default: 23
  --min-primer-melting-temp MINPRIMERMELT, -minprimermelt MINPRIMERMELT
                        minimal melting temperature of primers, degrees
                        Celsius. Default: 60
  --max-primer-melting-temp MAXPRIMERMELT, -maxprimermelt MAXPRIMERMELT
                        maximal melting temperature of primers, degrees
                        Celsius. Default: 68
  --optimal-primer-melting-temp OPTPRIMERMELT, -optprimermelt OPTPRIMERMELT
                        optimal melting temperature of primers, degrees
                        Celsius. Default: 64
  --min-primer-gc MINPRIMERGC, -minprimergc MINPRIMERGC
                        minimal acceptable GC-content for primers. Default: 20
  --max-primer-gc MAXPRIMERGC, -maxprimergc MAXPRIMERGC
                        maximal acceptable GC-content for primers. Default: 80
  --optimal-primer-gc OPTPRIMERGC, -optprimergc OPTPRIMERGC
                        optimal acceptable GC-content for primers. Default: 40
  --min-primer-end-gc MINPRIMERENDGC, -minprimerendgc MINPRIMERENDGC
                        minimal acceptable number of G or C nucleotides within
                        last 5 nucleotides of 3'-end of primers. Default: 0
  --max-primer-end-gc MAXPRIMERENDGC, -maxprimerendgc MAXPRIMERENDGC
                        maximal acceptable number of G or C nucleotides within
                        last 5 nucleotides of 3'-end of primers. Default: 5
  --opt-primer-end-gc OPTPRIMERENDGC, -optprimerendgc OPTPRIMERENDGC
                        optimal number of G or C nucleotides within last 5
                        nucleotides of 3'-end of primers. Default: 2
  --max-primer-poly-n MAXPRIMERPOLYN, -maxprimerpolyn MAXPRIMERPOLYN
                        maximal acceptable length of some poly-N in primers.
                        Default: 8
  --max-primer-compl-end-th MAXPRIMERCOMPLENDTH, -maxprimercomplendth MAXPRIMERCOMPLENDTH
                        maximal Tm for complementarity of 3'-ends of primers.
                        Default: 25
  --max-primer-compl-any-th MAXPRIMERCOMPLANYTH, -maxprimercomplanyth MAXPRIMERCOMPLANYTH
                        maximal Tm for any complementarity of primers.
                        Default: 35
  --max-primer-hairpin-th MAXPRIMERHAIRPINTH, -maxprimerhairpinth MAXPRIMERHAIRPINTH
                        maximal melting temperature of primer hairpin
                        structure. Default: 40
  --max-primer-nonspecific MAXPRIMERNONSPEC, -maxprimernonspec MAXPRIMERNONSPEC
                        maximal number of nonspecific regions to which primer
                        can hybridizes. Default: 10000
  --max-amplicons-overlap MAXOVERLAP, -maxoverlap MAXOVERLAP
                        maximal length of overlap between two amplified blocks
                        (it does not include primers). Default: 50
  --primers-number1 PRIMERNUM1, -primernum1 PRIMERNUM1
                        number of primer that user wants to get on the 1st
                        stage. The more this value, the more precise the
                        choice of primers, but the longer the design time.
                        Default: 50
  --auto-adjust-parameters, -autoadjust
                        use this parameter if you want NGS-PrimerPlex to
                        automatically use less stringent parameters if no
                        primer were constructed for some region
  --tries-to-get-best-combination TRIESTOGETCOMBINATION, -tries TRIESTOGETCOMBINATION
                        number of of tries to get the best primer combination.
                        More the value, better combination will be, but this
                        will take more time. Default: 10000
  --return-variants-number RETURNVARIANTSNUM, -returnvariantsnum RETURNVARIANTSNUM
                        number of multiplexes variants that user wants to get
                        after all analyses and filters. Default: 10
  --embedded-amplification, -embedded
                        use this parameter if you want to create NGS-panel
                        with embedded amplification
  --min-internal-primer-shift MINPRIMERSHIFT, -minprimershift MINPRIMERSHIFT
                        minimal shift of external primer from the 3'-end of
                        internal primer. Default: 5
  --opt-external-amplicon-length OPTEXTAMPLLEN, -optextampllen OPTEXTAMPLLEN
                        optimal length of the external amplicons. Default: 150
  --max-external-amplicon-length MAXEXTAMPLLEN, -maxextampllen MAXEXTAMPLLEN
                        maximal length of the external amplicons. Default: 150
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
  --gui, -gui           this parameter is only automatically used by GUI of
                        the application
```
### addSeqToPrimers.py
This script adds adapter sequences to all designed primers. As an input it uses NGS-primerplex.py output file and file with adapter sequences. Example files with adapter sequences are included into the repository. This script outputs sequences into new XLS-file listing all designed primers with names and adapter sequences added.
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
### convertToDraftFile.py
This script converts the main NGS-PrimerPlex output file (*_info.xls* file) into draft file for subsequent use them as draft-file for redesigning some primers or redistributing them into another multiplex sets. All parameters are listed below:
```
  -h, --help            show this help message and exit
  --input INFILE, -in INFILE
                        input XLS-file with primers designed by NGS-PrimerPlex
  --output OUTFILE, -out OUTFILE
                        file for outputing draft file for NGS-PrimerPlex
```
## Citation
Kechin A, Borobova V, Boyarskikh U, Khrapov E, Subbotin S, Filipenko M (2020) NGS-PrimerPlex: High-throughput primer design for multiplex polymerase chain reactions. PLoS Comput Biol 16(12): e1008468. https://doi.org/10.1371/journal.pcbi.1008468
