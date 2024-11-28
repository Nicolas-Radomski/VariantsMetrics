# Usage
The repository VariantsMetrics provides two R scripts (VariantsMetricsReference:1.0.R or VariantsMetricsWorkflows:1.0.R) to compute metrics (i.e. TP, TN, FP, FN, recall, precision, specificity, accuracy and Fscore) comparing genotypes ('undetected' stands for 'undetected genotype' and 'missing' stands for 'missing data') from expected (i.e. reference file or first workflow file) and predicted (i.e. workflow file or second workflow file) variants (i.e. SNP and InDels).
# Dependencies
The R scripts VariantsMetricsReference:1.0.R and VariantsMetricsWorkflows:1.0.R were prepared and tested with the R version 4.4.0 and Ubuntu 20.04 LTS Focal Fossa.
- require(remotes) # version 2.5.0
- library(argparse) # version 2.2.3
- the argparse library requires a sufficient Python binary (i.e. python3.12)
# Explaination
## Glossary
```
MD stands for missing data
MDr stands for missing data in the reference
MDw stands for missing data in the workflow
MDw1 stands for missing data in the first workflow
MDw2 stands for missing data in the second workflow
TP stands for true positive
TN stands for true negative
FP stands for false positive
FN stands for false negative
```
## Behavior
```
The workflows flag missing data as "missing" and the undetected positions as "undetected".
The users can also flag "missing" and "undefined" by themself.
```
## Controlled input variants
```
Reference      Workflow         Results    Positions (S1439.014)
1st Workflow   2nd Workflow
genotype       genotype         TP         254
genotype       other genotype   FP         653
undetected     genotype         FP         748
missing        genotype         MD         25014
genotype       undetected       FN         35124
genotype       missing          MD         7589
undetected     undetected       TN         47000
undetected     missing          MD         5012
missing        missing          MD         43000
missing        undetected       MD         45000
```
# Expected Input
## input_Reference_Variants.tsv
```
POS	GENOTYPE
254	AAA
653	T
7589	C
25014	
35124	-
43000	
45000	
47000	undetected
```
## input_Workflow_Variants.tsv
```
POS	S1239.014	S1339.020	S1339.010	S1442.017
254     AAA             AAA            AAA            AAA
653     G               T              G              T
748     G               undetected     G              undetected
5012                    undetected     missing        undetected
7589                    C              missing        C
25014   G                              G              missing
35124   undetected      -              undetected     -
43000                                  missing        missing
45000   undetected     undetected      undetected     undetected
47000   undetected     undetected      undetected     undetected
```
## input_First_Workflow_Variants.tsv
```
POS	S1239.014	S1442.017	S1339.020	S1339.010
254	AAA	        AAA	        AAA	        AAA
653	T	        T	        T	        T
7589	C	        C	        C	        C
25014				
35124	-	        -	        -	        -
43000				
45000				
47000	undetected	undetected	undetected	undetected
```
## input_Second_Workflow_Variants.tsv
```
POS     S1239.014	S1339.020	S1339.010	S1442.017
254     AAA            AAA             AAA            AAA
653     G              T               G              T
748     G              undetected      G              undetected
5012                   undetected      missing        undetected
7589                   C               missing        C
25014   G                              G              missing
35124   undetected     -               undetected     -
43000                                  missing        missing
45000   undetected     undetected      undetected     undetected
47000   undetected     undetected      undetected     undetected
```
# Usage
### VariantsMetricsReference:1.0.R
```
usage: /context-VariantsMetricsReference/VariantsMetricsReference:1.0.R
       [-h] -r CHARACTERS -w CHARACTERS -g INTEGER [-o CHARACTERS]
       [-b LOGICAL]

This script computes performance metrics (i.e. TP, TN, FP, FN, recall,
precision, specificity, accuracy and Fscore) comparing genotypes ('undetected'
stands for 'undetected genotype' and 'missing' stands for 'missing data') from
expected (i.e. reference file) and predicted (i.e. workflow file) variants
(i.e. SNP and InDels).

options:
  -h, --help            show this help message and exit
  -r CHARACTERS, --reference CHARACTERS
                        Reference input file with an absolute or relative path
                        (tab-separated values). First column: positions of
                        variants (header: whatever). Second column: profiles
                        of genotypes (header: whatever). [MANDATORY]
  -w CHARACTERS, --workflow CHARACTERS
                        Workflow input file with an absolute or relative path
                        (tab-separated values). First column: positions of
                        variants (header: whatever). Other columns: profiles
                        of genotypes (header: sample identifiers). [MANDATORY]
  -g INTEGER, --genome INTEGER
                        Size of the reference genome (bases). [MANDATORY]
  -o CHARACTERS, --prefix CHARACTERS
                        Absolute or relative output path with or without
                        output file prefix. [OPTIONAL, default output_]
  -b LOGICAL, --backup LOGICAL
                        Save an external representation of R objects (i.e.
                        saved_data.RData) and a short-cut of the current
                        workspace (i.e. saved_images.RData). [OPTIONAL,
                        default False]
```
### VariantsMetricsWorkflows:1.0.R
```
usage: /context-VariantsMetricsWorkflows/VariantsMetricsWorkflows:1.0.R
       [-h] -w1 CHARACTERS -w2 CHARACTERS -g INTEGER [-o CHARACTERS]
       [-b LOGICAL]

This script computes performance metrics (i.e. TP, TN, FP, FN, recall,
precision, specificity, accuracy and Fscore) comparing genotypes ('undetected'
stands for 'undetected genotype' and 'missing' stands for 'missing data') from
expected (i.e. first workflow file) and predicted (i.e. second workflow file)
variants (i.e. SNP and InDels).

options:
  -h, --help            show this help message and exit
  -w1 CHARACTERS, --workflow1 CHARACTERS
                        First workflow input file with an absolute or relative
                        path (tab-separated values). First column: positions
                        of variants (header: whatever). Other columns:
                        profiles of genotypes (header: sample identifiers
                        without duplicates and identical to those from the
                        '-w2/--workflow2' argument). [MANDATORY]
  -w2 CHARACTERS, --workflow2 CHARACTERS
                        Second workflow input file with an absolute or
                        relative path (tab-separated values). First column:
                        positions of variants (header: whatever). Other
                        columns: profiles of genotypes (header: sample
                        identifiers without duplicates and identical to those
                        from the '-w1/--workflow1' argument). [MANDATORY]
  -g INTEGER, --genome INTEGER
                        Size of the reference genome (bases). [MANDATORY]
  -o CHARACTERS, --prefix CHARACTERS
                        Absolute or relative output path with or without
                        output file prefix. [OPTIONAL, default output_]
  -b LOGICAL, --backup LOGICAL
                        Save an external representation of R objects (i.e.
                        saved_data.RData) and a short-cut of the current
                        workspace (i.e. saved_images.RData). [OPTIONAL,
                        default False]
```
# Usage
## Import the GitHub repository
```
git clone https://github.com/Nicolas-Radomski/VariantsMetrics.git
cd VariantsMetrics
```
## Import the Docker image
```
docker pull nicolasradomski/variantsmetricsreference:1.0
```
## Launch with R script
### VariantsMetricsReference:1.0.R
```
Rscript VariantsMetricsReference:1.0.R -r input_Reference_Variants.tsv -w input_Workflow_Variants.tsv -g 50000 -o output_Reference_
```
### VariantsMetricsWorkflows:1.0.R
```
Rscript VariantsMetricsWorkflows:1.0.R -w1 input_First_Workflow_Variants.tsv -w2 input_Second_Workflow_Variants.tsv -g 50000 -o output_Workflows_
```
## Launch with Docker
### VariantsMetricsReference:1.0.R
```
docker run --rm --name nicolas -u $(id -u):$(id -g) -v $(pwd):/wd nicolasradomski/variantsmetricsreference:1.0 -r input_Reference_Variants.tsv -w input_Workflow_Variants.tsv -g 50000 -o output_Reference_DockerHub_
```
### VariantsMetricsWorkflows:1.0.R
```
docker run --rm --name nicolas -u $(id -u):$(id -g) -v $(pwd):/wd nicolasradomski/variantsmetricsworkflows:1.0 -w1 input_First_Workflow_Variants.tsv -w2 input_Second_Workflow_Variants.tsv -g 50000 -o output_Workflows_DockerHub_
```
# Expected output
### VariantsMetricsReference:1.0.R
```
  samples MDr MDw MD TP    TN FP FN  recall precision  specificity accuracy  Fscore
S1239.014   3   3  5  1 49991  2  1 0.50000   0.33333      0.99996  0.99994 0.40000
S1339.010   3   3  5  1 49991  2  1 0.50000   0.33333      0.99996  0.99994 0.40000
S1339.020   3   2  3  4 49993  0  0 1.00000   1.00000      1.00000  1.00000 1.00000
S1442.017   3   2  3  4 49993  0  0 1.00000   1.00000      1.00000  1.00000 1.00000
```
### VariantsMetricsWorkflows:1.0.R
```
  samples MDw1 MDw2 MD TP    TN FP FN  recall precision  specificity accuracy  Fscore
S1239.014    3    3  5  1 49991  2  1 0.50000   0.33333     0.99996  0.99994 0.40000
S1339.010    3    3  5  1 49991  2  1 0.50000   0.33333     0.99996  0.99994 0.40000
S1339.020    3    2  3  4 49993  0  0 1.00000   1.00000     1.00000  1.00000 1.00000
S1442.017    3    2  3  4 49993  0  0 1.00000   1.00000     1.00000  1.00000 1.00000
```
# Illustration
![workflow figure](https://github.com/Nicolas-Radomski/VariantsMetrics/blob/main/illustration.png)
# References
- Docker Hub: https://hub.docker.com/r/nicolasradomski/variantsmetricsreference
- Docker Hub: https://hub.docker.com/r/nicolasradomski/variantsmetricsworkflows
# Acknowledgment
The GENPAT-IZSAM Staff
# Author
Nicolas Radomski
