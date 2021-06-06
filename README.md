# Pangolin

### Installation
* Install PyTorch: https://pytorch.org/
* Install other dependencies: 
  ```
  pip install pyvcf gffutils
  conda install -c conda-forge -c bioconda pandas pybedtools pyfaidx
* Install Pangolin: 
  ```
  git clone https://github.com/tkzeng/Pangolin.git
  cd Pangolin
  pip install .
  ```

### Usage

1. Create an annotation database from a GTF file, such as one from https://www.gencodegenes.org/, using `scripts/create_db.py`. Example usage:
   ```
   python scripts/create_db.py gencode.v34.annotation.gtf
   # output: gencode.v34.annotation.db  
   ```

2. Run Pangolin on a VCF or CSV file containing a list of variants. Under default settings, the maximum splice gain and loss scores within 50 bases of the variant, along with their positions, will be reported. Only substitutions and simple insertions/deletions are currently supported.
   Example usage:
   ```
   pangolin examples/brca.vcf GRCh37.primary_assembly.genome.fa gencode.v34.annotation.db brca_pangolin
   ```
   See full options below:
    
   ```
   positional arguments:
     variant_file          VCF or CSV file with a header (see COLUMN_IDS option).
     reference_file        FASTA file containing a reference genome sequence.
     annotation_file       gffutils database file. Can be generated using create_db.py.
     output_file           Prefix for output file. Will be a VCF/CSV if variant_file is VCF/CSV.
   
   optional arguments:
     -h, --help            show this help message and exit
     -c COLUMN_IDS, --column_ids COLUMN_IDS
                           (If variant_file is a CSV) Column IDs for: chromosome, variant position, reference bases, and alternative bases. Separate IDs by commas. (Default: CHROM,POS,REF,ALT)
     -m {False,True}, --mask {False,True}
                           If True, splice gains at annotated splice sites and splice losses at unannotated splice sites will be set to 0.
     -s SCORE_CUTOFF, --score_cutoff SCORE_CUTOFF
                           Output all sites with absolute predicted change in score >=cutoff, instead of only the maximum loss/gain sites.
     -d DISTANCE, --distance DISTANCE
                           Number of bases on either side of the variant for which splice scores should be calculated. (Default: 50)
   ```

### Recommendations