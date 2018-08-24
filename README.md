# pyRNAseq
Public RNAseq data analysis pipline based on python2 and R. The basic workflow is get the sra file based on GSM ID. Those sra file can be transfered to fastq with fastq-dump from [sratoolkit](https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/). The software [salmon](https://combine-lab.github.io/salmon/getting_started/) can quantify based on refseq ID. Differential expression analysis is implemented by DESeq2.

### Requirements
- python >= 2.6 and <=2.8
- R with DESeq2 installed
- [salmon](https://combine-lab.github.io/salmon/getting_started/)

### Installation and setup of pyRNAseq

	git clone https://github.com/WChangson/pyRNAseq.git
	cd pyRNAseq/source/

You can see a _config.json_ here. Open the file and set up the path of fastq-dump in the first line, the path of salmon in the second line, and the path of R  script you want to use in the third line. If you install those software globally, your can keep them as _fastq-dump_, _salmon_, _Rscript_. The fourth and fifth line is for the salmon index of hg38 and mm10.

	cd ..
	python setup.py install

### Usage of pyRNAseq
	pyRNAseq -h
This will list all the arguments of pyRNAseq
#### Options
##### -d/--dmat
The design matrix which should include two columns and delaminated by tab. The first column is GSM IDs of RNAseq samples. The second column is the design of experiment. Usually control group shoud be upper and treatment group should be at the bottom. For example:

GSM10001 | Control    |
GSM10002 | Control    |
GSM10003 | Treatment  |
GSM10004 | Treatment  |

##### -o/--output
Output directory of pyRNAseq results, including fastq files, read count files, and differential expression file. If this file do not exists, pyRNAseq will make it atomatically.

##### -s/--species [hg38, mm10]
We only support hg38 or mm10 for salmon's pseudo-alignment.

##### -sal/--salmon
Whether to run salmon.

##### -de/--difexpr
Whether to do differential expression analysis with DESeq2.
