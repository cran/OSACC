OSACC1 <-
function(Mrklist,CCdata,outputfile,OSACC1.risk=TRUE,sort.cov="ascending",checked.pval=0.05, perm.reps=10,statest="chisq")
{
print("OSACC version 1.0",quote=FALSE)

typetest<-c("logistic","trend","chisq")
statest<-match.arg(statest,typetest)

typesort<-c("ascending","descending")

sort.cov<-match.arg(sort.cov,typesort)

if(sort.cov=="ascending")
{ covsort<-FALSE
}else
if(sort.cov=="descending")
{
 covsort<-TRUE
}
if(OSACC1.risk) permway<-"OSACC1-risk" else permway<-"OSACC1-hom"

print(paste("You are doing",permway,"with",sort.cov,"covariate","by",statest,"test",sep=" "),quote=FALSE)

totaldata<-read.table(CCdata,header=FALSE,sep="\t")

map<-read.table(Mrklist,header=FALSE,sep="\t")

nummark<-dim(map)[1]
markname<-paste(c("allele1-","allele2-"),sort(c(1:nummark,1:nummark)),sep="")
names(totaldata)<-c("fam","ind","fid","mid","sex","status","trait",markname)


cleantrait<-totaldata[totaldata$trait!=-9999&!is.na(totaldata$trait),]
totaldata$trait<-NA
cleandata<-rbind(totaldata[totaldata$status==1,],cleantrait[cleantrait$status==2,])

selinfo<-NULL
for(j in 1:nummark)
{     
print(paste(map[j,1],"...",sep=" "),quote=FALSE) 
mkdata<-cleandata[,c("fam","ind","fid","mid","sex","status","trait",markname[2*j-1],markname[2*j])]
names(mkdata)<-c("fam","ind","fid","mid","sex","status","trait","allele1","allele2") 


#unaffected: remove unaffected with missing genotype or covariate, skip to nex mark if there are less than 20 unaffected.
dataunaff <-mkdata[mkdata$status ==1&(mkdata$allele1 +mkdata$allele2)!= 0, ]
if(dim(dataunaff)[1]<20) next

#affected: remove affected with missing genotype or covariate, skip to nex mark if there are less than 20 affected.
dataaff <-mkdata[mkdata$status == 2&(mkdata$allele1 +mkdata$allele2)!= 0,] 
if(dim(dataaff)[1]<20) next


#combine affected and unaffected together
datatest <- rbind(dataunaff,dataaff)
       	
numaff<-length(dataaff$fam)	#number of cases
numunaff<-length(dataunaff$fam) #number of controls
numind<-length(datatest$fam)   #number of all individuals.
ranks<-rank(datatest$trait,na.last="keep")
rkmkdata<-cbind(ranks,datatest)
names(rkmkdata)<-c("ranks","fam","ind","fid","mid","sex","status","trait","allele1","allele2")

pvalraw <-0.9999 # initial pvalraw


# A list of the unique ranks for loops
if(covsort)
{
uniqrk <- unique(sort(na.omit(ranks), decreasing =covsort))# sorted in descending
uniqrk <- ifelse(uniqrk <= numaff-20,uniqrk,NA) # remove subsets with less than 20 indivs
uniqrk <- na.omit(uniqrk) # remove missing values created by last line
} else
{
uniqrk <- unique(sort(na.omit(ranks), decreasing =covsort))# sorted in ascending
uniqrk <- ifelse(uniqrk >=20,uniqrk,NA) # remove subsets with less than 20 indivs
uniqrk <- na.omit(uniqrk) # remove missing values created by last line
}
for (ithrk in uniqrk)
{
if(covsort) updatedsubset<-rkmkdata[rkmkdata$ranks >= ithrk|rkmkdata$status==1,] else 
updatedsubset<-rkmkdata[rkmkdata$ranks <= ithrk|rkmkdata$status==1,]
alletable<-table(cbind(updatedsubset$status,updatedsubset$status),cbind(updatedsubset$allele1,updatedsubset$allele2))
if(dim(alletable)[1]*dim(alletable)[2]<4) next  # skip to next maker as it can not form 2X2 alletableency table
sc<-colSums(alletable)
sr<-rowSums(alletable)
E <- outer(sr, sc, "*")/sum(alletable)    # calculate the expect frquency for each cell to decide which test to be used
if(any(E<5)) pval<-fisher.test(alletable)$p.value else
 if(statest=="chisq")
{
pval<-chisq.test(alletable,correct=FALSE)$p.value
}else
if(statest=="logistic")
{
updatedsubset$add<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
updatedsubset$affection<- updatedsubset$status-1
logregfull<-glm(affection ~ add,family = binomial,data=updatedsubset)
restest<-summary(logregfull)
pval<-restest$coefficients[2,4]
}else
if(statest=="trend")
{
updatedsubset$snpgeno<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
genotrend<-table(updatedsubset$snpgeno,updatedsubset$status)
trendgenores<-tabletrend(genotrend)
pval<-trendgenores$p.value.bi
}
if(pval<pvalraw) 
{
   pvalraw<-pval
   cutrank<-ithrk
}
}
# to calucalte odds ratio etc. for the dataset which has smallest pvalue
if(covsort)
{ 
selsubset<-rkmkdata[rkmkdata$ranks >= cutrank|rkmkdata$status==1,]
cuttrait<-min(selsubset$trait,na.rm=TRUE) 
} else
{
selsubset <-rkmkdata[rkmkdata$ranks <= cutrank|rkmkdata$status==1,]
cuttrait<-max(selsubset$trait,na.rm=TRUE)
}
alletable<-table(cbind(selsubset$status,selsubset$status),cbind(selsubset$allele1,selsubset$allele2))
freq.control<- alletable[1,1]/sum(alletable[1,]) # freq in controls
freq.case <- alletable[2,1]/sum(alletable[2,]) # freq in cases
OR <-(freq.case-freq.control*freq.case)/(freq.control-freq.control*freq.case)
subnumaff <-sum(selsubset$status-1)
subnumunaff<-sum(2-selsubset$status)
sqlogOR <- (sum(1/(freq.control*subnumunaff*2)+1/(freq.case*subnumaff*2)+1/((1-freq.control)*subnumunaff*2)+1/((1-freq.case)*subnumaff*2)))^0.5 
upperCI <- OR*exp(sqlogOR*1.96)
lowerCI <- OR/exp(sqlogOR*1.96)
numaff.subset<-subnumaff
numunaff.subset<-subnumunaff
pval.raw<-pvalraw


##permutation test
countperm <-0
permresults <-0


if(pval.raw<checked.pval)
{
       
	   if(OSACC1.risk)
	   {

		for (i in 1:perm.reps) 
           	{ 
			samplenums <-sample(c(1:numind))  
          	 	samplemkdata<-rkmkdata[samplenums,]
 	  
          		permnewdata<-cbind(rkmkdata[,c("ranks","trait","status")],samplemkdata[,c("allele1","allele2")])
			names(permnewdata)<-c("ranks","trait","status","allele1","allele2")
           		pvalsample <-0.9999 # reset for each permutation

           		for (ithrk in uniqrk) 
           		{ 
             		  if(covsort)  updatedsubset<-permnewdata[permnewdata$ranks >= ithrk|permnewdata$status==1,] else 
			  updatedsubset<-permnewdata[permnewdata$ranks <= ithrk|permnewdata$status==1,]
              		  alletable<-table(cbind(updatedsubset$status,updatedsubset$status),cbind(updatedsubset$allele1,updatedsubset$allele2))
 			  if(dim(alletable)[1]*dim(alletable)[2]<4) next  # skip to next maker as it can not form 2X2 alletableency table
			  sc<-colSums(alletable)
			  sr<-rowSums(alletable)
			  E <- outer(sr, sc, "*")/sum(alletable)    # calculate the expect frquency for each cell to decide which test to be used
			  if(any(E<5)) pval<-fisher.test(alletable)$p.value else
 			  if(statest=="chisq")
			  {
			  pval<-chisq.test(alletable,correct=FALSE)$p.value
			  }else
			  if(statest=="logistic")
			  {
			  updatedsubset$add<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
			  updatedsubset$affection<- updatedsubset$status-1
			  logregfull<-glm(affection ~ add,family = binomial,data=updatedsubset)
			  restest<-summary(logregfull)
			  pval<-restest$coefficients[2,4]
			  }else
			  if(statest=="trend")
			  {
			  updatedsubset$snpgeno<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
			  genotrend<-table(updatedsubset$snpgeno,updatedsubset$status)
			  trendgenores<-tabletrend(genotrend)
			  pval<-trendgenores$p.value.bi
			  }
               		  pvalsample<-ifelse(pval<pvalsample,pval,pvalsample)
			}
			  countperm<-countperm +1 
             		  permresults[countperm] <- pvalsample 
             
		} 

         
		index<-ifelse(permresults<pval.raw,1,0)
         	pval.perm<- sum(index,na.rm =TRUE)/(countperm)  # countperm is real number of  permutation to get the pvalue
         	sd.pval <- sqrt(pval.perm*(1-pval.perm)/countperm)
	}

	if(!(OSACC1.risk))
	{
		for (i in 1:perm.reps) 
        	{ 
			affwithrank<-rkmkdata[rkmkdata$status==2,]
			samplenums <-sample(c(1:numaff)) 
			sampleaff<-affwithrank[samplenums,] 
          		permaff<-cbind(sampleaff[,c("ranks","trait")],affwithrank[,c("status","allele1","allele2")])
			names(permaff)<-c("ranks","trait","status","allele1","allele2")
			unaffwithrank<-rkmkdata[rkmkdata$status==1,c("ranks","trait","status","allele1","allele2")]
          		permnewdata<-rbind(permaff,unaffwithrank)
	  		names(permnewdata)<-c("ranks","trait","status","allele1","allele2")

           		pvalsample <-0.9999 # reset for each permutation

           		for (ithrk in uniqrk) 
           		{ 
             		  if(covsort)  updatedsubset<-permnewdata[permnewdata$ranks >= ithrk|permnewdata$status==1,] else 
			  updatedsubset<-permnewdata[permnewdata$ranks <= ithrk|permnewdata$status==1,]
              		  alletable<-table(cbind(updatedsubset$status,updatedsubset$status),cbind(updatedsubset$allele1,updatedsubset$allele2))
 			  if(dim(alletable)[1]*dim(alletable)[2]<4) next  #if the permuation data cannot form a 2X2 alletableency table, it skip to next 
			   sc<-colSums(alletable)
			   sr<-rowSums(alletable)
			   E <- outer(sr, sc, "*")/sum(alletable)    # calculate the expect frquency for each cell to decide which test to be used
			   if(any(E<5)) pval<-fisher.test(alletable)$p.value else
			   if(statest=="chisq")
			   {
			   pval<-chisq.test(alletable,correct=FALSE)$p.value
			   }else
			   if(statest=="logistic")
			   {
			   updatedsubset$add<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
			   updatedsubset$affection<- updatedsubset$status-1
			   logregfull<-glm(affection ~ add,family = binomial,data=updatedsubset)
			   restest<-summary(logregfull)
			   pval<-restest$coefficients[2,4]
			   }else
			   if(statest=="trend")
			   {
			   updatedsubset$snpgeno<-(2- updatedsubset$allele1)+(2- updatedsubset$allele2)
			   genotrend<-table(updatedsubset$snpgeno,updatedsubset$status)
			   trendgenores<-tabletrend(genotrend)
			   pval<-trendgenores$p.value.bi
			   }
			  pvalsample<-ifelse(pval<pvalsample,pval,pvalsample)
           		}
 			  countperm<-countperm +1  # increment counter
             		  permresults[countperm] <- pvalsample # store best p
		}

         	index<-ifelse(permresults<pval.raw,1,0)
         	pval.perm<- sum(index,na.rm =TRUE)/countperm  # countperm is real number of  permutation to get the pvalue
         	sd.pval <- sqrt(pval.perm*(1-pval.perm)/countperm)
	}
} else 
{
pval.perm<-pval.raw 
sd.pval<-1
}

seldatainfo<-data.frame(map[j,1],pval.raw,pval.perm,sd.pval,freq.case,freq.control,cuttrait,OR,lowerCI,upperCI,perm.reps,numaff,numaff.subset,numunaff,numunaff.subset)
selinfo<-rbind(selinfo,seldatainfo)

}
names(selinfo)<-c("marker","pval.raw","pval.perm","sd.pval","freq.case","freq.control","cuttrait","OR","lowerCI","upperCI","total_permutation","numaff","numaff.subset","numunaff","numunaff.subset")
write.table(selinfo,outputfile,row.names =FALSE,quote=FALSE,sep="\t")

}

