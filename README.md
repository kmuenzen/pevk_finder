# Downloading required Python modules
PEVK_Finder uses several modules that are not automatically installed with most versions of Python (Biopython and Statistics). The easiest way to install these packages is with the **pip** tool. If you do not already have pip installed on your machine, run the following commands:

1. Download the get-pip.py file. The file will be downloaded to the current directory, and you may delete this file after step #2 (below) is complete:
```bash
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
```
2. Run get-pip.py. For most users, it is recommended to install packages for the current user with the --user flag. If you have admin priveleges on your machine, you may also use the sudo command and enter the admin password:
```bash
python get-pip.py --user
```
or
```bash
sudo python get-pip.py
```

Once pip is installed on your machine, install the required Python modules using the following commands. The same sudo/user rules apply:

1. Install Biopython
```bash
pip install biopython --user
```
2. Install Statistics
```bash
pip install statistics --user
```


# Running the PEVK_Finder Program
## Inputs and Outputs
### Required Input
The PEVK_Finder program requires the filepath to a complete TTN nucleotide sequence (in fasta format), or the path to a directory containing one or more TTN sequences.

### Flags
**-i (--input)**<br/>
Path to the fasta file(s) that contain the full TTN nucleotide sequence.<br/>
<br/>
**-w (--window_length)**<br/>
The integer value of the sliding window length that PEVK_Finder will use to identify likely PEVK exons. **Default: 10**<br/>
<br/>
**-r (--pevk_ratio)**<br/>
The float value of the minimum PEVK ratio that PEVK_Finder will use to identify likely PEVK exons. **Default: 0.54**<br/>
<br/>
**-l (--exon_length)**<br/>
The integer value of the minimum PEVK exon length that PEVK_Finder will use to identify likely PEVK exons. **Default: 12**<br/>
<br/>
**-o (--outpath)**<br/>
The path to the directory where the output files will be stored. **Default: ./data/**<br/>
<br/>
**-v (--verbose)**<br/>
When verbose is True, will emit messages about script progress. **Default: True**<br/>
<br/>
**-p (--optional_outputs)**<br/>
When optional_outputs is True, will create the two optional output files containing exon length, PEVK ratio and coordinate information. **Default: False**<br/>
<br/>

### Output
Four output files per input sequence, containing exon sets with the following characteristics:

1. (species_name_)PEVK_exons_unbounded_AA(_w_r_l).fasta (Fasta file): ALL predicted PEVK exons as
    amino acid sequences, sorted by increasing TTN coordinates. **For all output files, coordinates in sequence descriptions
    are relative to full TTN DNA sequence.**
2. (species_name_)PEVK_exons_unbounded_NT(_w_r_l).fasta (Fasta file): ALL predicted PEVK exons as
    nucleotide sequences, sorted by increasing TTN coordinates.
3. (species_name_)PEVK_exons_bounded_AA(_w_r_l).fasta (Fasta file): IQR +- (1.5 x IQR) predicted PEVK exons as
    amino acid sequences, sorted by increasing TTN coordinates.
4. (species_name_)PEVK_exons_bounded_NT(_w_r_l).fasta (Fasta file): IQR +- (1.5 x IQR) predicted PEVK exons as
    nucleotide sequences, sorted by increasing TTN coordinates.

Two optional output files per input sequence:

1. (species_name_)exon_coordinates(_w_r_l).csv (CSV file): The start and end coordinates of each predicted PEVK exon.
1. (species_name_)exon_lengths_and_ratios(_w_r_l).csv (CSV file): The exon length and PEVK ratio for each predicted PEVK exon.

## Running PEVK_Finder in the Command Line

To run PEVK_Finder on **one TTN nucleotide sequence** using default parameters, run the following command from the pevk_finder_public directory:
```bash
python -W ignore pevk_finder_v_1.py -i [species_name]_ttn.fasta
```

To run PEVK_Finder on **multiple TTN sequences** using default parameters, run the following command from the pevk_finder_public directory:
```bash
python -W ignore pevk_finder_v_1.py -i [directory_with_tnn_seqs]
```
