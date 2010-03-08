OSACC <-
function(Mrklist,CCdata,outputfile,OSACC.test="OSACC1",OSACC.risk=TRUE,sort.cov="ascending",checked.pval=0.05,perm.reps=1000,mark.out=0)
{

if (OSACC.test=="OSACC1") n1<-1 else 
if(OSACC.test=="OSACC2") n1<-0
else 
{
print ("You need select right OSACC test 'OSACC1' or 'OSACC2'")
break
}
if (OSACC.risk) n2<-1 else n2<-0
if(sort.cov=="ascending") n3<-0 else
if(sort.cov=="descending") n3<-1
else
{
print ("You need select a way to sort the covariate either by 'ascending' or 'descending'")
break
}
n4<-checked.pval
n5<-perm.reps
n6<-mark.out
filenames<-c(Mrklist,CCdata,outputfile)
param<-c(n1,n2,n3,n4,n5,n6)
 .C("OSACC",as.character(filenames),as.double(param),PACKAGE="OSACC")
 res<-read.table(outputfile, sep="\t",header=TRUE)
allinfo<-list(CCdata, Mrklist,res)
names(allinfo)<-c("CCdata","Mrklist","output")
return(allinfo)
}

