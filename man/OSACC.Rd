\name{OSACC}
\alias{OSACC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ordered Subset Analysis for Case-Control Studies
}
\description{
This OSACC function is written in C++ code and implemented in an R function. This function runs much faster than the other two functions OSACC1, OSACC2, but only the allelic chi-squared test statistic for the 2x2 contingency table of affection status vs. allele is available. The OSACC function can perform the OSACC1-risk, OSACC1-hom, OSACC2-risk and OSACC2-hom tests.
}
\usage{
OSACC(Mrklist, CCdata, outputfile, OSACC.test = "OSACC1", 
      OSACC.risk = TRUE, sort.cov = "ascending", 
      checked.pval = 0.05, perm.reps = 1000, mark.out = 0)
}


\arguments{
  \item{Mrklist}{A tab-delimited file that lists all SNP markers included in the CCdata file without a header line. Each row in this file r
epresents one marker and includes three fields: marker name, chromosome number, map position (physical or genetic).
}
  \item{CCdata}{A tab-delimited file without a header line. It has a standard pre-MAKEPED pedigree file format, even though it includes onl
y data from unrelated individuals. Only a single OSACC covariate may be included in a single CCdata file, and all markers have to be bialle
lic. The file contains the following fields in each row:
 
Fam_ID: Family ID, which could be the same for all individuals or could take on a different value for each individual. It is not used by th
e program, which simply assumes that each row contains independent observations.

Person_ID: individual ID.
 
Father_ID: This is typically zero for case-control data and is not used by the program.
 
Mother_ID:  This is typically zero for case-control data and is not used by the program.

Sex: Male=1, Female=2
 
Phenotype: Case=2, Control=1. Individuals that cannot be classified as either a case or a control should not be included in the CCdata file
. 
 
Covariate: Value of the covariate used in the OSACC sorting procedure. Any numeric value is allowed, and -9999 codes for missing data. 
 
mk1_1 mk1_2 mk2_1 mk2_2 .......: Allele 1 at Marker 1, Allele 2 at Marker 1, Allele 1 at Marker 2, Allele 2 at Marker 2 etc. Only biallelic
 markers are allowed  (coded as 1, 2, or 0 for missing data). The order of the markers in the CCdata file has to match that in the Mrklist 
file. 

}
  \item{outputfile}{A file name for storing program output.

}

  \item{OSACC.test}{There are two options, one is "OSACC1" which uses covariates of cases only, the other is "OSACC2" which uses covariates of cases and controls. The default is OSACC.test="OSACC1"
}
  \item{OSACC.risk}{A logical value indicating OSACC1-risk test or OSACC1-hom test,the default is OSACC1.risk=TRUE. OSACC1.risk=FALSE will perform the OSACC1-hom test.
}
  \item{sort.cov}{Sort covariate values in "ascending" order, or in "descending" order (case only), the default sort.cov="ascending".
}
  \item{checked.pval}{p-value threshold for performing permutations, the default is 0.05.
}
  \item{perm.reps}{number of permutations, the default is 1000.
}
  \item{mark.out}{The options are 0,1,2. This specifies the criteria for generating marker-specific output that lists the identified subset of cases (OSACC1) or cases and controls (OSACC2). The first (0) does not generate any subset listings for any marker, regardless of OSACC permutation values; the second (1) generates subset listings only for markers with permutation p-value<0.05; the third (2) uses a more stringent criterion by generating subset listings only for markers with Bonferroni-corrected permutation p-value<0.05. The default mark.out=o.
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{marker}{marker name}
\item{pval.raw}{ smallest observed p-value (from specified test statistic) obtained for any data subset, based on specified covariate order}
\item{pval.permp}{ empirical p-value from OSACC permutation test. This p-value should be used to evaluate whether the identified subset shows significantly increased association evidence, compared to the evidence obtained from an analysis of the entire dataset.}
\item{sd.pval}{standard deviation of permutation test p-value, which depends on the number of permutations that were performed. If pval.raw is larger than the specified threshold for performing a permutation test (ie, no permutation test was performed for this marker),this value will be 1.} 
\item{freq.case}{allele frequency in selected subset of cases.}
\item{freq.control}{allele frequency in selected subset of controls. For OSACC1, this is the allele frequency for all controls.}
\item{cuttrait}{covariate cutoff value for selected cases subset.}
\item{OR}{allelic odds ratio for selected subset of individuals.}
\item{lowerCI,upperCI}{95 percent confidence interval for OR.}
\item{total_permutation}{the total number of  permutations performed for the marker.}
\item{numaff}{number of cases with non-missing SNP and covariate data in the entire dataset.}
\item{numaff.subset}{number of cases with non-missing SNP and covariate data in selected subset.}
\item{numunaff}{number of controls with non-missing SNP and covariate data in the entire dataset.}
\item{numunaff.subset}{number of controls with non-missing SNP and covariate data in selected subset (same as numunaff for OSACC1).}
}
\references{
Xuejun Qin, Elizabeth R. Hauser, Silke Schmidt (2010), Ordered Subset Analysis for Case-Control Studies, Genetic Epidemiology 34(5), 407-417. 
}
\author{Xuejun Qin, Elizabeth R. Hauser, Silke Schmidt, (Center for Human Genetics, Duke University)
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%%data(Mrklist,package="OSACC")
%%data(CCdata,package="OSACC")
\donttest{OSACC("Mrklist.txt","CCdata.txt","outputfile",OSACC.test="OSACC1",
		OSACC.risk=TRUE,sort.cov="ascending",check.pval=0.05,
		perm.reps=1000,mark.out=0)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 
