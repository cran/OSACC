\name{OSACC2}
\alias{OSACC2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordered Subset Analysis for Case-Control Studies (two)
}
\description{
The Ordered Subset Association Analysis for case-control studies incorporates covariate values of cases only (OSACC1-hom, OSACC1-risk) or of cases and controls (OSACC2-hom, OSACC2-risk). The default test statistic is the chi-squared test for the 2x2 contingency table of affection status vs. allele. Options for genotype-based Cochran-Armitage trend test and the Wald chi-squared coefficient test in a logistic regression model with additively coded genotypes are available. Other tests may be performed by editing the source code.
}
\usage{OSACC2(Mrklist, CCdata, outputfile, OSACC2.risk=TRUE,
 sort.cov="ascending", checked.pval=0.05, perm.reps=10, statest="chisq")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mrklist}{A tab-delimited file that lists all SNP markers included in the CCdata file without a header line. Each row in this file represents one marker and includes three fields:

marker name

chromosome number

map position (physical or genetic).
}
  \item{CCdata}{A tab-delimited file without a header line. It has a standard pre-MAKEPED pedigree file format, even though it includes only data from unrelated individuals. Only a single OSACC covariate may be included in a single CCdata file, and all markers have to be biallelic. The file contains the following fields in each row:
 
Fam_ID: Family ID, which could be the same for all individuals or could take on a different value for each individual. It is not used by the program, which simply assumes that each row contains independent observations.

Person_ID: individual ID.
 
Father_ID: This is typically zero for case-control data and is not used by the program.
 
Mother_ID:  This is typically zero for case-control data and is not used by the program.

Sex: Male=1, Female=2
 
Phenotype: Case=2, Control=1. Individuals that cannot be classified as either a case or a control should not be included in the CCdata file. 
 
Covariate: Value of the covariate used in the OSACC sorting procedure. Any numeric value is allowed, and -9999 codes for missing data. 
 
mk1_1 mk1_2 mk2_1 mk2_2 .......: Allele 1 at Marker 1, Allele 2 at Marker 1, Allele 1 at Marker 2, Allele 2 at Marker 2 etc. Only biallelic markers are allowed  (coded as 1, 2, or 0 for missing data). The order of the markers in the CCdata file has to match that in the Mrklist file. 

}
  \item{outputfile}{A file name for storing program output.

}
  \item{OSACC2.risk}{A logical value indicating OSACC2-hom test or OSACC2-risk test, the default is OSACC2.risk=TRUE.
}
  \item{sort.cov}{The default is "ascending", please specify sort.cov="descending" otherwise.
}
  \item{checked.pval}{p-value threshold for performing permutations, the default is 0.05.
}
  \item{perm.reps}{number of permutations, the default is 1000.
}
  \item{statest}{The default is statest="chisq", it is chi-square test for 2x2 contingency table of affection status and number of two allele. The option 'statest="trend"' implements the genotype-based Cochran-Armitage trend test. The option 'statest="logist"' implements the Wald chi-square test for the genotype model coefficient in a logistic regression, the code for  genotyoe is in additive way (0,1,2),  it is a continous variable . If there is  more covariate need in logistic regression, you need modified R code.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{marker}{marker name}
\item{pval.raw}{ smallest observed p-value (from specified test statistic) obtained for any data subset, based on specified covariate order}
\item{pval.permp}{ empirical p-value from OSACC permutation test. This p-value should be used to evaluate whether the identified subset shows significantly increased association evidence, compared to the evidence obtained from an analysis of the entire dataset.}
\item{sd.pval}{standard deviation of permutation test p-value, which depends on the number of permutations that were performed. If pval.raw is larger than the specified threshold for performing a permutation test (ie, no permutation test was performed for this marker), this value will be 1.} 
\item{freq.case}{allele frequency in selected subset of cases.}
\item{freq.control}{allele frequency in selected subset of controls.}
\item{cuttrait}{covariate cutoff value for selected subset from both case and control.}
\item{OR}{allelic odds ratio for selected subset of individuals.}
\item{lowerCI,upperCI}{95 percent confidence interval for OR.}
\item{total_permutation}{the total number of  permutations performed for the marker.}
\item{numaff}{number of cases with non-missing SNP and covariate data in the entire dataset.}
\item{numaff.subset}{number of cases with non-missing SNP and covariate data in selected subset.}
\item{numunaff}{number of controls with non-missing SNP and covariate data in the entire dataset.}
\item{numunaff.subset}{number of controls with non-missing SNP and covariate data in selected subset.}
}
\author{Xuejun Qin, Elizabeth R. Hauser, Silke Schmidt, (Center for Human Genetics, Duke University)
}
\references{Xuejun Qin, Elizabeth R. Hauser, Silk Schmidt (2010),Ordered Subset Analysis for Case-Control Studies, Genetic Epidemiology 34(5), 407-417. 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%%OSACC2{\link{help}}, ~~~
}
\examples{
%%data(Mrklist,package="OSACC")
%%data(Ccdata,package="OSACC")
\donttest{OSACC2("Mrklist.txt", "CCdata.txt", "temp.txt", OSACC2.risk=TRUE, 
		sort.cov="ascending", checked.pval=0.05,  perm.reps=5, statest="chisq")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
