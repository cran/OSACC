\name{CCdata}
\docType{data}
\alias{CCdata}
\title{A case-control with marker data set}
\description{A tab-delimited file without a header line. It has a standard pre-MAKEPED pedigree file format, even though it includes
only data from unrelated individuals. Only a single OSACC covariate may be included in a single CCdata file, and all markers have to b
e biallelic. The file contains the following fields in each row:
Fam_ID: Family ID, which could be the same for all individuals or could take on a different value for each individual. It is not used b
y the program, which simply assumes that each row contains independent observations.

Person_ID: individual ID.
 
Father_ID: This is typically zero for case-control data and is not used by the program.
 
Mother_ID:  This is typically zero for case-control data and is not used by the program.

Sex: Male=1, Female=2
 
Phenotype: Case=2, Control=1. Individuals that cannot be classified as either a case or a control should not be included in the CCdata 
file. 
 
Covariate: Value of the covariate used in the OSACC sorting procedure. Any numeric value is allowed, and -9999 codes for missing data. 
 
mk1_1 mk1_2 mk2_1 mk2_2 .......: Allele 1 at Marker 1, Allele 2 at Marker 1, Allele 1 at Marker 2, Allele 2 at Marker 2 etc. Only biall
elic markers are allowed  (coded as 1, 2, or 0 for missing data). The order of the markers in the CCdata file has to match that in the 
Mrklist file. 

}
\usage{CCdata}
\format{A tab-delimited file without a header line}
\keyword{CCdata}
