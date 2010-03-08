#include "OSACC.h"


extern "C" void OSACC(char* namefile[], double parinput[])
//main()
{
//The main program is in c++ , the header change in compatible to c
	 
	//input parameters
	// float checked_pval=0.20;
     //int sort_order=0,totperm=200;//these parameter will be in makefile.
	double checked_pval;
	int sort_order,totperm;

	 
	FILE *res1;
	FILE *res2;
	FILE *res3=fopen(namefile[2],"w");
	FILE *LogFile =fopen("OSACC.log", "w");
	float igp[BUFFER],ind[BUFFER],isex[BUFFER];
	float istatus[BUFFER],status[BUFFER],pstatus[BUFFER];
	float itrait[BUFFER],trait[BUFFER],afftrait[BUFFER];
	float iallele1[BUFFER],allele1[BUFFER],pallele1[BUFFER],affallele1[BUFFER];
	float iallele2[BUFFER],allele2[BUFFER],pallele2[BUFFER],affallele2[BUFFER];
	float ranktrait[BUFFER];
	float uniqtrait[BUFFER];
	float chisqval,pval,permpval,df,cff,cmv,itemp,tempnum,cuttrait,empval,sdempval,bonfpval;
	double minpval,minpermpval;
	int irownum=0,num_row,num_uniqtrait,numaff,numunaff,numrow_cell=2,numcol_cell=2,num_mrk,iout;
	int i,j,k,cutrank,icell,jcell,imrk,iperm,num_smallpval,num_sucessperm,index,whole_cov,permu_way,pull_marker;
	static int cellfrq[3][3];
	int **frq;

	int *nrowp=&numrow_cell,*ncolp=&numcol_cell;
	double expected = -1.0;
	double percnt = 100.0;
	double emin = 0.0;
	double prt = 0.0;
	double *expectedp = &expected;
	double *percntp = &percnt;
	double *eminp = &emin;
	double *prtp = &prt;
	double pvalue;
	double *pvaluep = &pvalue;
	int workspace = 1000000;
	int *workspacep = &workspace;
	double table[4],minnumcell=100;
	int sumcell,sumrowcell[2],sumcolcell[2],currentnumaff,currentnumunaff;
	static int *midpoint[3]={&cellfrq[0][0],&cellfrq[1][0],&cellfrq[2][0]};
	double freq_cntrol,freq_case,OR,selogOR,upperCI,lowerCI,pval_raw;
	char *marker_names[BUFFER];


	    if(parinput[3]>1 ||parinput[3]<0)
	      { 
		cout<<"\n\n You input wrong parameter, it should be in (0 1]\n";
		 exit(0);
	      }

	    else checked_pval=parinput[3];
	    
      

	    if(parinput[4]!=(int)parinput[4])
	      { 
		cout<<"\n\n Warning: You input a non-integer parameter ";
	      }
	    totperm=(int)parinput[4];
	      

	    if( parinput[0]==0 && parinput[1]==0 && parinput[2]==0)
	      {
		cout<<"\n\nYou are doing OSACC2-hom test by ascending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==0 && parinput[1]==0 && parinput[2]==1)
	      {
		cout<<"\n\nYou are doing OSACC2-hom test by descending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==0 && parinput[1]==1 && parinput[2]==0)
	      {
		cout<<"\n\nYou are doing OSACC2-risk test by ascending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==0 && parinput[1]==1 && parinput[2]==1)
	      {
		cout<<"\n\nYou are doing OSACC2-risk test by descending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==1 && parinput[1]==0 && parinput[2]==0)
	      {
		cout<<"\n\nYou are doing OSACC1-hom test by ascending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==1 && parinput[1]==0 && parinput[2]==1)
	      {
		cout<<"\n\nYou are doing OSACC1-hom test by descending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==1 && parinput[1]==1 && parinput[2]==0)
	      {
		cout<<"\n\nYou are doing OSACC1-risk test by ascending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else if( parinput[0]==1 && parinput[1]==1 && parinput[2]==1)
	      {
		cout<<"\n\nYou are doing OSACC1-risk test by descending covariate.";
		cout<<"\nThe total permutation is "<< totperm <<" for those markers with p-value less than "<< checked_pval;
		cout<<"\n\n\n";
		
	      }
	    else {
	      cout<<"\n\n You input wrong parameter\n\n";
	    exit(0);
	    }
			   whole_cov=(int)parinput[0];
			   permu_way=(int)parinput[1];
			   sort_order=(int)parinput[2];
			   pull_marker=(int)parinput[5];
			   

 
	  

	num_mrk=read_dat_data(namefile[0],marker_names);	
cout<<"\nTotal marker: "<<num_mrk<<"\n";
 if (num_mrk<=0)
   {
     cout<<"\n\n There is no marker in your data file\n";
		 exit(0);
   }

if(whole_cov==1)
{       

fprintf(res3,"mark\tpval.raw\tpval.perm\tsd.pval\tfreq.case\tfreq.control\tcuttrait\tOR\tlowerCI\tupperCI\ttotal_permutation\tnumaff\tnumaff.subset\tnumunaff\tnumunaff.subset\n");


for(imrk=0;imrk<num_mrk;imrk++)
{

	index=2*imrk;
	make_ped_file(namefile[1],index,&irownum,num_mrk);

	res1=fopen("./temp_ped_file.ped","r");
	
	cout<<"Marker "<<marker_names[imrk]<<" ...\n";
	//	cout<<"read temp ped file\n";


	num_uniqtrait=0, numaff=0,numunaff=0;
	minpval=0.99;
	cellfrq[0][0]=0;//no use
	cellfrq[0][1]=0;//no use
	cellfrq[0][2]=0;//no use
	cellfrq[1][0]=0;//no use
	cellfrq[1][1]=0;//use
	cellfrq[1][2]=0;//use
	num_row=0;
	for (i=0;i<irownum;i++)
	{  
	//fread(&osaccdata[i],sizeof (struct osacc_type),1,res1);
	fscanf(res1,"%f\t",&igp[i]);
	fscanf(res1,"%f\t",&ind[i]);
	fscanf(res1,"%f\t",&isex[i]);
	fscanf(res1,"%f\t",&istatus[i]);// get data set for one SNP data informat staus, trait, and two alleles
	fscanf(res1,"%f\t",&itrait[i]);// unaffected  coded  1, affected coded 2.
	fscanf(res1,"%f\t",&iallele1[i]);// missing value for trait is -9999. missing value for allele is 0.
	fscanf(res1,"%f\t",&iallele2[i]);
	fscanf(res1,"\n");
	if(iallele1[i]!=0&&iallele2[i]!=0)
	{
		if(itrait[i]!=-9999||istatus[i]==1)
		{
			status[num_row]=istatus[i];
			trait[num_row]=itrait[i];
			allele1[num_row]=iallele1[i];
			allele2[num_row]=iallele2[i];
			num_row++;
		}
	}
	}

	fclose(res1);
	r.init();
	myperm.init(irownum,10);

 for (i=0;i<num_row;i++)
   { 
     if(status[i]==1&&allele1[i]==1&&allele2[i]==1)
      
       {
		cellfrq[1][1]+=2;
		numunaff++;
       }
     if(status[i]==1&&allele1[i]==1&&allele2[i]==2)
       {
		cellfrq[1][1]++;
		cellfrq[1][2]++;
		numunaff++;
       }
     if(status[i]==1&&allele1[i]==2&&allele2[i]==1)
       {
		cellfrq[1][1]++;
		cellfrq[1][2]++;
		numunaff++;
       }
     if(status[i]==1&&allele1[i]==2&&allele2[i]==2)
       {
		cellfrq[1][2]+=2;
		numunaff++;
       }
	if(status[i]==2&&trait[i]!=-9999&&allele1[i]!=0)
		{ 
		afftrait[numaff]=trait[i];
		affallele1[numaff]=allele1[i];
		affallele2[numaff]=allele2[i];
		numaff++;
		}
	}
 // cout<<"before sort cov\n";

 if (numaff<20)
   { cout<<"There is less than 20 affected after removing missing data, skip to next marker\n";
     continue;
   }

if (numunaff<20)
   { cout<<"There is less than 20 unaffected after removing missing data, skip to next marker\n";
     continue;
   }
 cout<<"number of affected  for  this marker "<<numaff<<"\n";
 cout<<"\n";

	sort3(numaff,afftrait-1, affallele1-1, affallele2-1);
	//cout<<"after sort cov\n";
	if(sort_order==0)
	{
		crank(numaff,afftrait-1,ranktrait-1,uniqtrait-1,&num_uniqtrait);

		for (j=0;j<num_uniqtrait;j++)
		{

			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use


			for (k=0;afftrait[k]<=uniqtrait[j]&&k<numaff;k++)
			{
				if( affallele1[k]==1&&affallele2[k]==1)
				cellfrq[2][1]+=2;
				if(affallele1[k]==1&&affallele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if( affallele1[k]==2&&affallele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(affallele1[k]==2&&affallele2[k]==2)
				cellfrq[2][2]+=2;
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			minnumcell=100;
			for (icell=0;icell<2;icell++)
			{
				for(jcell=0;jcell<2;jcell++)
				{
					tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
					if(tempnum<minnumcell) minnumcell=tempnum;
				}
			}
			if (minnumcell<5)
			{
				table[0]=cellfrq[1][1];
				table[1]=cellfrq[2][1];
				table[2]=cellfrq[1][2];
				table[3]=cellfrq[2][2];
				fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
				pval=pvalue;
			}
			else {
			frq=midpoint;
			cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
			}
			if(pval<minpval)
			{
				minpval=pval;
				cutrank=j;
			}
		}
			cuttrait=uniqtrait[cutrank];
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use


			for (k=0;afftrait[k]<=uniqtrait[cutrank]&&k<numaff;k++)
			{
				if( affallele1[k]==1&&affallele2[k]==1)
				cellfrq[2][1]+=2;
				if(affallele1[k]==1&&affallele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if( affallele1[k]==2&&affallele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(affallele1[k]==2&&affallele2[k]==2)
				cellfrq[2][2]+=2;
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			freq_cntrol=(double)cellfrq[1][1]/sumrowcell[0] ;// freq in controls
			freq_case=(double)cellfrq[2][1]/sumrowcell[1] ;//req in cases
			OR =(freq_case-freq_cntrol*freq_case)/(freq_cntrol-freq_cntrol*freq_case);  // Odd ratio
 
			currentnumaff=k;
			selogOR = pow(1/(freq_cntrol*numunaff*2)+1/(freq_case*currentnumaff*2)+1/((1-freq_cntrol)*numunaff*2)+1/((1-freq_case)*currentnumaff*2),0.5); 
			upperCI = OR*exp(selogOR*1.96); // get CI
			lowerCI = OR/exp(selogOR*1.96);
			pval_raw=minpval;
 

	}


	if(sort_order==1)
	{
		for (i=0;i<(int)numaff/2;i++)
		{
			SWAP(afftrait[i],afftrait[numaff-i-1]);
			SWAP(affallele1[i],affallele1[numaff-i-1]);
			SWAP(affallele2[i],affallele2[numaff-i-1]);
		}

		crank(numaff,afftrait-1,ranktrait-1,uniqtrait-1,&num_uniqtrait);

		for (j=0;j<num_uniqtrait;j++)
		{

			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use


			for (k=0;afftrait[k]>=uniqtrait[j]&&k<numaff;k++)
			{
				if( affallele1[k]==1&&affallele2[k]==1)
				cellfrq[2][1]+=2;
				if(affallele1[k]==1&&affallele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if( affallele1[k]==2&&affallele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(affallele1[k]==2&&affallele2[k]==2)
					cellfrq[2][2]+=2;
			}
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			minnumcell=100;
			for (icell=0;icell<2;icell++)
			{
				for(jcell=0;jcell<2;jcell++)
				{
					tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
					if(tempnum<minnumcell) minnumcell=tempnum;
				}
			}
			if (minnumcell<5)
			{	table[0]=cellfrq[1][1];
				table[1]=cellfrq[2][1];
				table[2]=cellfrq[1][2];
				table[3]=cellfrq[2][2];
				fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
				pval=pvalue;
			}
			else {
			frq=midpoint;
			cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
			}
			if(pval<minpval)
			{
				minpval=pval;
				cutrank=j;
			}
		}
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use

			cuttrait=uniqtrait[cutrank];
			for (k=0;afftrait[k]>=uniqtrait[cutrank]&&k<numaff;k++)
			{
				if( affallele1[k]==1&&affallele2[k]==1)
				cellfrq[2][1]+=2;
				if(affallele1[k]==1&&affallele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if( affallele1[k]==2&&affallele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(affallele1[k]==2&&affallele2[k]==2)
				cellfrq[2][2]+=2;
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			freq_cntrol=(double)cellfrq[1][1]/sumrowcell[0] ;// freq in controls
			freq_case=(double)cellfrq[2][1]/sumrowcell[1] ;//req in cases
			OR =(freq_case-freq_cntrol*freq_case)/(freq_cntrol-freq_cntrol*freq_case);  // Odd ratio
 
			currentnumaff=k;
			selogOR = pow(1/(freq_cntrol*numunaff*2)+1/(freq_case*currentnumaff*2)+1/((1-freq_cntrol)*numunaff*2)+1/((1-freq_case)*currentnumaff*2),0.5); 
			upperCI = OR*exp(selogOR*1.96); // get CI
			lowerCI = OR/exp(selogOR*1.96);
			pval_raw=minpval;
 

	}

	//	cout<<"before permutation\n";
if(pval_raw<checked_pval)
{

	num_smallpval=0;//count the number of reps permuation which has  smaller  pvalue than raw data
	num_sucessperm=0;// count the number  of  reps permuation which has success to permuation.

	for (iperm=0;iperm<totperm;iperm++)
	{
	
	 numaff=0;
	 // if(permu_way==1)
	 //  permutation1(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 if(permu_way==0)
	   permutation2(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 // if(permu_way==3)
	 //  permutation3(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 if(permu_way==1)
	   permutation4(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
			
	  //	permutation(allele1,allele2,pallele1,pallele2,num_row);
	      
		cellfrq[0][0]=0;//no use
		cellfrq[0][1]=0;//no use
		cellfrq[0][2]=0;//no use
		cellfrq[1][0]=0;//no use
		cellfrq[1][1]=0;//use
		cellfrq[1][2]=0;//use
	
		minpermpval=0.98;
		
		for (i=0;i<num_row;i++)
		{ 
			if(pstatus[i]==1&&pallele1[i]==1&&pallele2[i]==1)
      
			{
				cellfrq[1][1]+=2;
			}
			if(pstatus[i]==1&&pallele1[i]==1&&pallele2[i]==2)
			{
				cellfrq[1][1]++;
				cellfrq[1][2]++;
			}
			if(pstatus[i]==1&&pallele1[i]==2&&pallele2[i]==1)
			{
				cellfrq[1][1]++;
				cellfrq[1][2]++;
			}
			if(pstatus[i]==1&&pallele1[i]==2&&pallele2[i]==2)
			{
				cellfrq[1][2]+=2;
			}
			if(pstatus[i]==2)
			{ 
				affallele1[numaff]=pallele1[i];
				affallele2[numaff]=pallele2[i];
				numaff++;
			}
		}
 
		if(sort_order==0)
		{

			for (j=0;j<num_uniqtrait;j++)
			{

				cellfrq[2][0]=0;//no use
				cellfrq[2][1]=0;//use
				cellfrq[2][2]=0;//use


				for (k=0;afftrait[k]<=uniqtrait[j]&&k<numaff;k++)
				{
					if( affallele1[k]==1&&affallele2[k]==1)
						cellfrq[2][1]+=2;
					if(affallele1[k]==1&&affallele2[k]==2)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if( affallele1[k]==2&&affallele2[k]==1)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(affallele1[k]==2&&affallele2[k]==2)
						cellfrq[2][2]+=2;
				}

		
					sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
					sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
					sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
					sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
					sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
					minnumcell=100;
					for (icell=0;icell<2;icell++)
					{
						for(jcell=0;jcell<2;jcell++)
						{
							tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
							if(tempnum<minnumcell) minnumcell=tempnum;
						}
					}
					if (minnumcell<5)
					{
						table[0]=cellfrq[1][1];
						table[1]=cellfrq[2][1];
						table[2]=cellfrq[1][2];
						table[3]=cellfrq[2][2];
						fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
						permpval=pvalue;
					}
					else {
						frq=midpoint;
						cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
						permpval=pval;
					}
					if(permpval<minpermpval)
					{
						minpermpval=permpval;
					}
		
			}
		}


		if(sort_order==1)
		{
	

			for (j=0;j<num_uniqtrait;j++)
			{

				cellfrq[2][0]=0;//no use
				cellfrq[2][1]=0;//use
				cellfrq[2][2]=0;//use


				for (k=0;afftrait[k]>=uniqtrait[j]&&k<numaff;k++)
				{
					if( affallele1[k]==1&&affallele2[k]==1)
						cellfrq[2][1]+=2;
					if(affallele1[k]==1&&affallele2[k]==2)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if( affallele1[k]==2&&affallele2[k]==1)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(affallele1[k]==2&&affallele2[k]==2)
						cellfrq[2][2]+=2;
				}
				sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
				sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
				sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
				sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
				sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
				minnumcell=100;
				for (icell=0;icell<2;icell++)
				{
					for(jcell=0;jcell<2;jcell++)
					{
						tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
						if(tempnum<minnumcell) minnumcell=tempnum;
					}
				}
				if (minnumcell<5)
				{	table[0]=cellfrq[1][1];
					table[1]=cellfrq[2][1];
					table[2]=cellfrq[1][2];
					table[3]=cellfrq[2][2];
					fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
					permpval=pvalue;
				}
				else {
					frq=midpoint;
					cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
					permpval=pval;
				}
			

				if(permpval<minpermpval)
				{
					minpermpval=permpval;
				}
			}
		}

	num_sucessperm++;
	if(minpermpval<=pval_raw) num_smallpval++;

	}	


	empval=(float)num_smallpval/(float)num_sucessperm;
	sdempval=(float)sqrt(empval*(1-empval)/(float)num_sucessperm);
}
else 
{
	empval=pval_raw;
	sdempval=1;
}

 char	filenames[40]="OSACC1_";
 strcat(filenames,marker_names[imrk]);
 if(pull_marker==1 && empval<0.05)
   {
     res2=fopen(filenames,"w"); 
     for (iout=0;iout<irownum;iout++)
       {
	 if(iallele1[iout]!=0&&iallele2[iout]!=0)
	   {
	     if(sort_order==0)
	       {
		 if(itrait[iout]<=cuttrait||istatus[iout]==1)
		   {
		     fprintf(res2,"%d\t",(int)igp[iout]);
		     fprintf(res2,"%d\t",(int)ind[iout]);
		     fprintf(res2,"%d\t",(int)isex[iout]);
		     fprintf(res2,"%d\t",(int)istatus[iout]);
		     fprintf(res2,"%f\t",itrait[iout]);
		     fprintf(res2,"%d\t",(int)iallele1[iout]);
		     fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		   }
	       }
	     if(sort_order==1)
	       {
		 if(itrait[iout]>=cuttrait||istatus[iout]==1)
		   {
		     fprintf(res2,"%d\t",(int)igp[iout]);
		     fprintf(res2,"%d\t",(int)ind[iout]);
		     fprintf(res2,"%d\t",(int)isex[iout]);
		     fprintf(res2,"%d\t",(int)istatus[iout]);
		     fprintf(res2,"%f\t",itrait[iout]);
		     fprintf(res2,"%d\t",(int)iallele1[iout]);
		     fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		   }
	       }
	   }
       }fclose(res2);
   }
else if(pull_marker==2)
  {
    bonfpval=0.05/num_mrk;
    if (empval<bonfpval)   
      { 
	res2=fopen(filenames,"w"); 
      for (iout=0;iout<irownum;iout++)
	{
	  if(iallele1[iout]!=0&&iallele2[iout]!=0)
	    {
	      if(sort_order==0)
		{
		  if(itrait[iout]<=cuttrait||istatus[iout]==1)
		    {
		      fprintf(res2,"%d\t",(int)igp[iout]);
		      fprintf(res2,"%d\t",(int)ind[iout]);
		      fprintf(res2,"%d\t",(int)isex[iout]);
		      fprintf(res2,"%d\t",(int)istatus[iout]);
		      fprintf(res2,"%f\t",itrait[iout]);
		      fprintf(res2,"%d\t",(int)iallele1[iout]);
		      fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		    }
		}
	      
	      if(sort_order==1)
		{
		  if(itrait[iout]>=cuttrait||istatus[iout]==1)
		    {
		      fprintf(res2,"%d\t",(int)igp[iout]);
		      fprintf(res2,"%d\t",(int)ind[iout]);
		      fprintf(res2,"%d\t",(int)isex[iout]);
		      fprintf(res2,"%d\t",(int)istatus[iout]);
		      fprintf(res2,"%f\t",itrait[iout]);
		      fprintf(res2,"%d\t",(int)iallele1[iout]);
		      fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		    }
		}
	    }
	}fclose(res2);
      }
  }
 
                fprintf(res3,"%s\t",marker_names[imrk]);
		fprintf(res3,"%10e\t",pval_raw);
		fprintf(res3,"%10e\t",empval);
		fprintf(res3,"%f\t",sdempval);
		fprintf(res3,"%f\t",freq_case);
		fprintf(res3,"%f\t",freq_cntrol);
		fprintf(res3,"%f\t",cuttrait);
		fprintf(res3,"%f\t",OR);
		fprintf(res3,"%f\t",lowerCI);
		fprintf(res3,"%f\t",upperCI);
                fprintf(res3,"%d\t",totperm);
		fprintf(res3,"%d\t",numaff);
		fprintf(res3,"%d\t",currentnumaff);
		fprintf(res3,"%d\t",numunaff);
                fprintf(res3,"%d\n",numunaff);

}

}
else if(whole_cov==0)
{ 


fprintf(res3,"mark\tpval.raw\tpval.perm\tsd.pval\tfreq.case\tfreq.control\tcuttrait\tOR\tlowerCI\tupperCI\ttotal_permutation\tnumaff\tnumaff.subset\tnumunaff\tnumunaff.subset\n");

for(imrk=0;imrk<num_mrk;imrk++)
{

	index=2*imrk;
	make_ped_file(namefile[1],index,&irownum,num_mrk);

	//cout<<"Marker "<<imrk+1<<" ...\n";
	res1=fopen("./temp_ped_file.ped","r");
	cout<<"Marker "<<marker_names[imrk]<<" ...\n";	

	num_uniqtrait=0, numaff=0,numunaff=0;
	minpval=0.99;
	
	num_row=0;
	for (i=0;i<irownum;i++)
	{  
	//fread(&osaccdata[i],sizeof (struct osacc_type),1,res1);	
	  fscanf(res1,"%f\t",&igp[i]);
	  fscanf(res1,"%f\t",&ind[i]);
	  fscanf(res1,"%f\t",&isex[i]);
	  fscanf(res1,"%f\t",&istatus[i]);// get data set for one SNP data informat staus, trait, and two alleles
	  fscanf(res1,"%f\t",&itrait[i]);// unaffected  coded  1, affected coded 2.
	  fscanf(res1,"%f\t",&iallele1[i]);// missing value for trait is -9999. missing value for allele is 0.
	  fscanf(res1,"%f\t",&iallele2[i]);
	  fscanf(res1,"\n");
	if(iallele1[i]!=0&&iallele2[i]!=0)
	{
		if(itrait[i]!=-9999)
		{
			status[num_row]=istatus[i];
			trait[num_row]=itrait[i];
			allele1[num_row]=iallele1[i];
			allele2[num_row]=iallele2[i];
			num_row++;
			if(istatus[i]==2) numaff++;
			if(istatus[i]==1) numunaff++;
		}
	}
	}

	fclose(res1);
	r.init();
 	myperm.init(irownum,10);
 cout<<"Number of sampled individuals for this marker: "<<num_row<<"\n";
 cout<<"\n";
  if (numaff<20)
   { cout<<"There is less than 20 affected after removing missing data, skip to next marker\n";
     continue;
   }

if (numunaff<20)
   { cout<<"There is less than 20 unaffected after removing missing data, skip to next marker\n";
     continue;
   }

	sort3wholedata(num_row,trait-1,status-1,allele1-1,allele2-1);

	if(sort_order==0)

	{
	
	crankwholedata(num_row,trait-1,ranktrait-1,status-1,uniqtrait-1,&num_uniqtrait);
	
		for (j=0;j<num_uniqtrait;j++)
		{
			cellfrq[0][0]=0;//no use
			cellfrq[0][1]=0;//no use
			cellfrq[0][2]=0;//no use
			cellfrq[1][0]=0;//no use
			cellfrq[1][1]=0;//use
			cellfrq[1][2]=0;//use
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use


			for (k=0;trait[k]<=uniqtrait[j]&&k<num_row;k++)
			{
				if(status[k]==2&&allele1[k]==1&&allele2[k]==1)
				cellfrq[2][1]+=2;
				if(status[k]==2&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==2)
					cellfrq[2][2]+=2;

				if(status[k]==1&&allele1[k]==1&&allele2[k]==1)
      
				{
					cellfrq[1][1]+=2;
		
				}
				if(status[k]==1&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[1][2]+=2;
				}
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			minnumcell=100;
			for (icell=0;icell<2;icell++)
			{
				for(jcell=0;jcell<2;jcell++)
				{
					tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
					if(tempnum<minnumcell) minnumcell=tempnum;
				}
			}
			if (minnumcell<5)
			{
				table[0]=cellfrq[1][1];
				table[1]=cellfrq[2][1];
				table[2]=cellfrq[1][2];
				table[3]=cellfrq[2][2];
				fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
				pval=pvalue;
			}
			else {
			frq=midpoint;
			cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
			}
			if(pval<minpval)
			{
				minpval=pval;
				cutrank=j;
			}
		}
			cuttrait=uniqtrait[cutrank];
			cellfrq[1][0]=0;//no use
			cellfrq[1][1]=0;//use
			cellfrq[1][2]=0;//use
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use
currentnumunaff=0;
currentnumaff=0;

			for (k=0;trait[k]<=uniqtrait[cutrank]&&k<num_row;k++)
			{
				
				
				
				if(status[k]==2&&allele1[k]==1&&allele2[k]==1)
				{
					cellfrq[2][1]+=2;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[2][2]+=2;
					currentnumaff++;
				}

				if(status[k]==1&&allele1[k]==1&&allele2[k]==1)
      
				{
					cellfrq[1][1]+=2;

					currentnumunaff++;
				}
				if(status[k]==1&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
					currentnumunaff++;
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
					currentnumunaff++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[1][2]+=2;
					currentnumunaff++;
				}
				
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			freq_cntrol=(double)cellfrq[1][1]/sumrowcell[0] ;// freq in controls
			freq_case=(double)cellfrq[2][1]/sumrowcell[1] ;//req in cases
			OR =(freq_case-freq_cntrol*freq_case)/(freq_cntrol-freq_cntrol*freq_case);  // Odd ratio

			selogOR = pow(1/(freq_cntrol*currentnumunaff*2)+1/(freq_case*currentnumaff*2)+1/((1-freq_cntrol)*currentnumunaff*2)+1/((1-freq_case)*currentnumaff*2),0.5); 
			upperCI = OR*exp(selogOR*1.96); // get CI
			lowerCI = OR/exp(selogOR*1.96);
			pval_raw=minpval;
 

	}


	if(sort_order==1)
	{
		for (i=0;i<(int)num_row/2;i++)
		{
			SWAP(trait[i],trait[num_row-i-1]);
			SWAP(status[i],status[num_row-i-1]);
			SWAP(allele1[i],allele1[num_row-i-1]);
			SWAP(allele2[i],allele2[num_row-i-1]);
		}

		crankwholedata(num_row,trait-1,ranktrait-1,status-1,uniqtrait-1,&num_uniqtrait);

		for (j=0;j<num_uniqtrait;j++)
		{
			cellfrq[0][0]=0;//no use
			cellfrq[0][1]=0;//no use
			cellfrq[0][2]=0;//no use
			cellfrq[1][0]=0;//no use
			cellfrq[1][1]=0;//use
			cellfrq[1][2]=0;//use
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use


			for (k=0;trait[k]>=uniqtrait[j]&&k<num_row;k++)
			{
				if(status[k]==2&&allele1[k]==1&&allele2[k]==1)
				cellfrq[2][1]+=2;
				if(status[k]==2&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==2)
					cellfrq[2][2]+=2;

				if(status[k]==1&&allele1[k]==1&&allele2[k]==1)
      
				{
					cellfrq[1][1]+=2;
		
				}
				if(status[k]==1&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[1][2]+=2;
				}
			}
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			minnumcell=100;
			for (icell=0;icell<2;icell++)
			{
				for(jcell=0;jcell<2;jcell++)
				{
					tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
					if(tempnum<minnumcell) minnumcell=tempnum;
				}
			}
			if (minnumcell<5)
			{	table[0]=cellfrq[1][1];
				table[1]=cellfrq[2][1];
				table[2]=cellfrq[1][2];
				table[3]=cellfrq[2][2];
				fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
				pval=pvalue;
			}
			else {
			frq=midpoint;
			cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
			}
			if(pval<minpval)
			{
				minpval=pval;
				cutrank=j;
			}
		}
			cellfrq[1][0]=0;//no use
			cellfrq[1][1]=0;//use
			cellfrq[1][2]=0;//use
			cellfrq[2][0]=0;//no use
			cellfrq[2][1]=0;//use
			cellfrq[2][2]=0;//use
			cuttrait=uniqtrait[cutrank];
			currentnumunaff=0;
			currentnumaff=0;
			for (k=0;trait[k]>=uniqtrait[cutrank]&&k<num_row;k++)
			{
					if(status[k]==2&&allele1[k]==1&&allele2[k]==1)
				{
					cellfrq[2][1]+=2;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[2][1]++;
					cellfrq[2][2]++;
					currentnumaff++;
				}
				if(status[k]==2&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[2][2]+=2;
					currentnumaff++;
				}

				if(status[k]==1&&allele1[k]==1&&allele2[k]==1)
      
				{
					cellfrq[1][1]+=2;

					currentnumunaff++;
				}
				if(status[k]==1&&allele1[k]==1&&allele2[k]==2)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
					currentnumunaff++;
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==1)
				{
					cellfrq[1][1]++;
					cellfrq[1][2]++;
					currentnumunaff++;
				
				}
				if(status[k]==1&&allele1[k]==2&&allele2[k]==2)
				{
					cellfrq[1][2]+=2;
					currentnumunaff++;
				}
			}

		
			sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
			sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
			sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
			sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
			sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
			freq_cntrol=(double)cellfrq[1][1]/sumrowcell[0] ;// freq in controls
			freq_case=(double)cellfrq[2][1]/sumrowcell[1] ;//req in cases
			OR =(freq_case-freq_cntrol*freq_case)/(freq_cntrol-freq_cntrol*freq_case);  // Odd ratio
 
		
			selogOR = pow(1/(freq_cntrol*currentnumunaff*2)+1/(freq_case*currentnumaff*2)+1/((1-freq_cntrol)*currentnumunaff*2)+1/((1-freq_case)*currentnumaff*2),0.5); 
			upperCI = OR*exp(selogOR*1.96); // get CI
			lowerCI = OR/exp(selogOR*1.96);
			pval_raw=minpval;
 

	}


if(pval_raw<checked_pval)
{

	num_smallpval=0;//count the number of reps permuation which has  smaller  pvalue than raw data
	num_sucessperm=0;// count the number  of  reps permuation which has success to permuation.

	for (iperm=0;iperm<totperm;iperm++)
	{
	
	  //	 if(permu_way==1)
	  //	 permutation1(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 if(permu_way==0)
	   permutation2(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 // if(permu_way==3)
	 //permutation3(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 if(permu_way==1)
	   permutation4(allele1,allele2,pallele1,pallele2,status,pstatus,num_row);
	 
	    minpermpval=0.98;
		
	
 
		if(sort_order==0)
		{

			for (j=0;j<num_uniqtrait;j++)
			{
				cellfrq[0][0]=0;//no use
				cellfrq[0][1]=0;//no use
				cellfrq[0][2]=0;//no use
				cellfrq[1][0]=0;//no use
				cellfrq[1][1]=0;//use
				cellfrq[1][2]=0;//use
				cellfrq[2][0]=0;//no use
				cellfrq[2][1]=0;//use
				cellfrq[2][2]=0;//use


				for (k=0;trait[k]<=uniqtrait[j]&&k<num_row;k++)
				{
					if(pstatus[k]==2&&pallele1[k]==1&&pallele2[k]==1)
						cellfrq[2][1]+=2;
					if(pstatus[k]==2&&pallele1[k]==1&&pallele2[k]==2)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(pstatus[k]==2&&pallele1[k]==2&&pallele2[k]==1)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(pstatus[k]==2&&pallele1[k]==2&&pallele2[k]==2)
						cellfrq[2][2]+=2;

					if(pstatus[k]==1&&pallele1[k]==1&&pallele2[k]==1)
      
					{
						cellfrq[1][1]+=2;
		
					}
					if(pstatus[k]==1&&pallele1[k]==1&&pallele2[k]==2)
					{
						cellfrq[1][1]++;
						cellfrq[1][2]++;
				
					}
					if(pstatus[k]==1&&pallele1[k]==2&&pallele2[k]==1)
					{
						cellfrq[1][1]++;
						cellfrq[1][2]++;
				
					}
					if(pstatus[k]==1&&pallele1[k]==2&&pallele2[k]==2)
					{
						cellfrq[1][2]+=2;
					}
				}

		
					sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
					sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
					sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
					sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
					sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
					minnumcell=100;
					for (icell=0;icell<2;icell++)
					{
						for(jcell=0;jcell<2;jcell++)
						{
							tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
							if(tempnum<minnumcell) minnumcell=tempnum;
						}
					}
 
					if (minnumcell<5)
					{
						table[0]=cellfrq[1][1];
						table[1]=cellfrq[2][1];
						table[2]=cellfrq[1][2];
						table[3]=cellfrq[2][2];
						fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
						permpval=pvalue;
					}
					else {
				
						frq=midpoint;
						cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
						permpval=pval;
					      
					}
					if(permpval<minpermpval)
					{
						minpermpval=permpval;
					}
		

			}
		}


		if(sort_order==1)
		{
	

			for (j=0;j<num_uniqtrait;j++)
			{
				cellfrq[0][0]=0;//no use
				cellfrq[0][1]=0;//no use
				cellfrq[0][2]=0;//no use
				cellfrq[1][0]=0;//no use
				cellfrq[1][1]=0;//use
				cellfrq[1][2]=0;//use
				cellfrq[2][0]=0;//no use
				cellfrq[2][1]=0;//use
				cellfrq[2][2]=0;//use


				for (k=0;trait[k]>=uniqtrait[j]&&k<num_row;k++)
				{
					if(pstatus[k]==2&&pallele1[k]==1&&pallele2[k]==1)
						cellfrq[2][1]+=2;
					if(pstatus[k]==2&&pallele1[k]==1&&pallele2[k]==2)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(pstatus[k]==2&&pallele1[k]==2&&pallele2[k]==1)
					{
						cellfrq[2][1]++;
						cellfrq[2][2]++;
					}
					if(pstatus[k]==2&&pallele1[k]==2&&pallele2[k]==2)
						cellfrq[2][2]+=2;

					if(pstatus[k]==1&&pallele1[k]==1&&pallele2[k]==1)
      
					{
						cellfrq[1][1]+=2;
		
					}
					if(pstatus[k]==1&&pallele1[k]==1&&pallele2[k]==2)
					{
						cellfrq[1][1]++;
						cellfrq[1][2]++;
				
					}
					if(pstatus[k]==1&&pallele1[k]==2&&pallele2[k]==1)
					{
						cellfrq[1][1]++;
						cellfrq[1][2]++;
				
					}
					if(pstatus[k]==1&&pallele1[k]==2&&pallele2[k]==2)
					{
						cellfrq[1][2]+=2;
					}
				}
				sumcell=cellfrq[1][1]+cellfrq[1][2]+cellfrq[2][1]+cellfrq[2][2];
				sumrowcell[0]=cellfrq[1][1]+cellfrq[1][2];
				sumrowcell[1]=cellfrq[2][1]+cellfrq[2][2];
				sumcolcell[0]=cellfrq[1][1]+cellfrq[2][1];
				sumcolcell[1]=cellfrq[1][2]+cellfrq[2][2];
				minnumcell=100;
				for (icell=0;icell<2;icell++)
				{
					for(jcell=0;jcell<2;jcell++)
					{
						tempnum=(double)sumrowcell[icell]*sumcolcell[jcell]/sumcell;
						if(tempnum<minnumcell) minnumcell=tempnum;
					}
				}
				if (minnumcell<5)
				{	table[0]=cellfrq[1][1];
					table[1]=cellfrq[2][1];
					table[2]=cellfrq[1][2];
					table[3]=cellfrq[2][2];
					fexact(nrowp,ncolp,table,nrowp,expectedp,percntp,eminp,prtp,pvaluep,workspacep);
					permpval=pvalue;
				}
				else {
					frq=midpoint;
					cntab1(frq,numrow_cell,numcol_cell,&chisqval,&df,&pval,&cmv,&cff);
					permpval=pval;
				}
			

				if(permpval<minpermpval)
				{
					minpermpval=permpval;
				}
			}
		}

	num_sucessperm++;
	if(minpermpval<=pval_raw) num_smallpval++;

	}	


	empval=(float)num_smallpval/(float)num_sucessperm;
	sdempval=(float)sqrt(empval*(1-empval)/(float)num_sucessperm);
}
else 
{
	empval=pval_raw;
	sdempval=1;
}
 char	filenames[40]="OSACC2_";
 strcat(filenames,marker_names[imrk]);
 if(pull_marker==1 && empval<0.05)
   {
     res2=fopen(filenames,"w"); 
     for (iout=0;iout<irownum;iout++)
       {
	 if(iallele1[iout]!=0&&iallele2[iout]!=0)
	   {
	     if(sort_order==0)
	       {
		 if(itrait[iout]<=cuttrait)
		   {
		     fprintf(res2,"%d\t",(int)igp[iout]);
		     fprintf(res2,"%d\t",(int)ind[iout]);
		     fprintf(res2,"%d\t",(int)isex[iout]);
		     fprintf(res2,"%d\t",(int)istatus[iout]);
		     fprintf(res2,"%f\t",itrait[iout]);
		     fprintf(res2,"%d\t",(int)iallele1[iout]);
		     fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		   }
	       }
	     if(sort_order==1)
	       {
		 if(itrait[iout]>=cuttrait)
		   {
		     fprintf(res2,"%d\t",(int)igp[iout]);
		     fprintf(res2,"%d\t",(int)ind[iout]);
		     fprintf(res2,"%d\t",(int)isex[iout]);
		     fprintf(res2,"%d\t",(int)istatus[iout]);
		     fprintf(res2,"%f\t",itrait[iout]);
		     fprintf(res2,"%d\t",(int)iallele1[iout]);
		     fprintf(res2,"%d\n",(int)iallele2[iout]);
		   }
	       }
	   }
       }fclose(res2);
   }
else if(pull_marker==2)
  {
    bonfpval=0.05/num_mrk;
    
    if (empval<bonfpval) 
    {
      res2=fopen(filenames,"w"); 
      for (iout=0;iout<irownum;iout++)
	{
	  if(iallele1[iout]!=0&&iallele2[iout]!=0)
	    {
	      if(sort_order==0)
		{
		  if(itrait[iout]<=cuttrait)
		    {
		      fprintf(res2,"%d\t",(int)igp[iout]);
		      fprintf(res2,"%d\t",(int)ind[iout]);
		      fprintf(res2,"%d\t",(int)isex[iout]);
		      fprintf(res2,"%d\t",(int)istatus[iout]);
		      fprintf(res2,"%f\t",itrait[iout]);
		      fprintf(res2,"%d\t",(int)iallele1[iout]);
		      fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		    }
		}
	      
	      if(sort_order==1)
		{
		  if(itrait[iout]>=cuttrait)
		    {
		      fprintf(res2,"%d\t",(int)igp[iout]);
		      fprintf(res2,"%d\t",(int)ind[iout]);
		      fprintf(res2,"%d\t",(int)isex[iout]);
		      fprintf(res2,"%d\t",(int)istatus[iout]);
		      fprintf(res2,"%f\t",itrait[iout]);
		      fprintf(res2,"%d\t",(int)iallele1[iout]);
		      fprintf(res2,"%d\n",(int)iallele2[iout]);
		     
		    }
		}
	    }
	}
    fclose(res2);
    }
  }


                fprintf(res3,"%s\t",marker_names[imrk]);
		fprintf(res3,"%10e\t",pval_raw);
		fprintf(res3,"%10e\t",empval);
		fprintf(res3,"%f\t",sdempval);
		fprintf(res3,"%f\t",freq_case);
		fprintf(res3,"%f\t",freq_cntrol);
		fprintf(res3,"%f\t",cuttrait);
		fprintf(res3,"%f\t",OR);
		fprintf(res3,"%f\t",lowerCI);
		fprintf(res3,"%f\t",upperCI);
                fprintf(res3,"%d\t",totperm);
                fprintf(res3,"%d\t",numaff);
                fprintf(res3,"%d\t",currentnumaff);
                fprintf(res3,"%d\t",numunaff);
		fprintf(res3,"%d\n",currentnumunaff);



	}
}
fclose(res3);

fclose(LogFile);
cout<<"\n Done!\n";
return;
 
}
void permutation1(float allele1[], float allele2[], float pallele1[],float pallele2[], float status[], float pstatus[],int numind)
{
	unsigned  int *result;
	int i,tempint;

	result = myperm.newperm(r.ran2());				// feed some new random number as argument to get a different/random permutation
	//myperm.display(result);							// just so you see what you got...for debugging only
	for (i=0;i<numind;i++)
	{
	  while(result[i]>=numind)
	    {
	      tempint=result[i];
	      result[i]=result[tempint];
	    }
		pallele1[i]=allele1[result[i]];
		pallele2[i]=allele2[result[i]];
		pstatus[i]=status[i];
	}
	
	//delete [] result;
	free(result);								// don't forget to free the memory once it is no longer need. Otherwise you'll have a memory leak
}													// which, depending on run time, can lead to nasty crashes/segfault.

void permutation2(float allele1[], float allele2[], float pallele1[],float pallele2[], float status[], float pstatus[],int numind)
{
	unsigned  int *result;
	int i,tempint;

	result = myperm.newperm(r.ran2());	       	// feed some new random number as argument to get a different/random permutation
      						      
	for (i=0;i<numind;i++)
	{
	  while(result[i]>=numind)
	    {
	      tempint=result[i];
	      result[i]=result[tempint];
	    }
		pallele1[i]=allele1[result[i]];
		pallele2[i]=allele2[result[i]];
		pstatus[i]=status[result[i]];

	}
	
	//delete [] result;
	free(result);		       

}					       // which, depending on run time, can lead to nasty crashes/segfault.

void permutation3(float allele1[], float allele2[], float pallele1[],float pallele2[], float status[], float pstatus[],int numind)
{
	unsigned  int *result;
	int i,tempint;
	       
	result = myperm.newperm(r.ran2());	       	// feed some new random number as argument to get a different/random permutation
	for (i=0;i<numind;i++)
	{
	  while(result[i]>=numind)
	    {
	      tempint=result[i];
	      result[i]=result[tempint];
	    }
		pstatus[i]=status[result[i]];
		pallele1[i]=allele1[i];
		pallele2[i]=allele2[i];
	}
	
	//delete [] result;
	free(result);		       	// don't forget to free the memory once it is no longer need. Otherwise you'll have a memory leak
}						  // which, depending on run time, can lead to nasty crashes/segfault.


void permutation4(float allele1[], float allele2[], float pallele1[],float pallele2[], float status[], float pstatus[],int numind)
{
	unsigned  int *result;
	int i,tempint;

	result = myperm.newperm(r.ran2());	       	// feed some new random number as argument to get a different/random permutation
      						      
	for (i=0;i<numind;i++)
	{
	  while(result[i]>=numind)
	    {
	      tempint=result[i];
	      result[i]=result[tempint];
	    }
		pallele1[i]=allele1[result[i]];
		pallele2[i]=allele2[result[i]];
	}
	
	//delete [] result;
	free(result);		       
	result = myperm.newperm(r.ran2());	       	// feed some new random number as argument to get a different/random permutation
	for (i=0;i<numind;i++)
	{
	  while(result[i]>=numind)
	    {
	      tempint=result[i];
	      result[i]=result[tempint];
	    }
		pstatus[i]=status[result[i]];
	}
	
	//delete [] result;
	free(result);		       	// don't forget to free the memory once it is no longer need. Otherwise you'll have a memory leak
}						  // which, depending on run time, can lead to nasty crashes/segfault.









int read_dat_data(char datafile[],char *marker_names[])
{


ifstream dat(datafile);
  char line[BUFFER];
  char *token;
  int index = 0;
  int total_markers;
char *marker_name; 
total_markers=0;

 if(dat.fail()) {
     cerr << "Error opening file " << datafile << "\nexiting..." << endl;
     ofstream LogFile("OSACC.log", ios::app);
     LogFile<< "Error opening file " << datafile<< "\nexiting..." << endl;
     LogFile.close();
     exit(0);
  }

 //dat.getline(line,BUFFER);

  while(dat.getline(line,BUFFER))
  { 
	 marker_name = new char[100];

	token = strtok(line,"\t");
	if(token!=NULL||token[0]=='R'||token[0]=='r') 	
	  {	strcpy(marker_name,token);
		marker_names[index]=marker_name;
		index++;
		total_markers++;
	  }
		else
		strcpy(marker_name," ");
	//	cout<<marker_name<<"\n";
	
	int num_token = 0;
  	while(token!=NULL) {
		num_token ++;
		token = strtok(NULL,"\t");

	
		}
	if(num_token!=3)
	{
	cout<<"The line "<<total_markers<<" maybe missing.The format should  have  three columns 'marker_name','chrom','location'\nexiting..." << endl;
	ofstream LogFile("OSACC.log", ios::app);
     LogFile<<"The format should  have  three columns 'chrom','marker_name','location\nexiting..." << endl;
     LogFile.close();
     exit(0);
		
	}

  }
  // cout<<total_markers<<"\n";
  
  dat.close();
  return total_markers;
 
}


void make_ped_file(char *pedfile,int index,int *numind, int num_mrk)

{
  // cout<<index;
  // cout<<"in";
  int i;	
  if(index==0)
	 {
	    ifstream ped_check(pedfile);
	    char *token;
	    int num_token;
	    char line[BUFFER];
	    if(ped_check.fail()) {
	      cerr << "Error opening file " << pedfile << "\nexiting..." << endl;
	      ofstream LogFile("OSACC.log", ios::app);
	      LogFile<< "Error opening file " << pedfile << "\nexiting..." << endl;
	      LogFile.close();
	      exit(0);
	    }

	    int num_row = 0;
	    while (ped_check.getline(line,BUFFER))
	      {
		token = strtok(line,"\t");
		 num_token = 0;

		while(token!=NULL) 
		  {
		    num_token++;
		    token = strtok(NULL,"\t");
		  }
		if(num_token!=7+2*(num_mrk)) 
		  {
		  cout<<"Error: Number of markers mismatch the number of columns in "<<pedfile<<" row "<<num_row+1<<".\n";
		  ofstream LogFile("OSACC.log", ios::app);
		  LogFile<< " Number of markers mismatch the number of columns in"<<pedfile<<".\n";
		  LogFile.close();
		  exit(0);
		  }
		num_row++;
	      }
	    ped_check.close();
	  }
	//cout<<num_row<<"\n";


 int num_row = 0;

char *read_token = new char[4];

	ifstream ped(pedfile);
	ofstream pedFile("./temp_ped_file.ped");


 while(ped>>read_token) {
	pedFile<<read_token<<"\t";
	ped >>read_token;
	pedFile<<read_token<<"\t";
    for(i=0;i<2;i++) {
      ped >>read_token;
    }
     	ped>>read_token;
    	pedFile<<read_token<<"\t";
	ped>>read_token;
	pedFile<<read_token<<"\t";
	ped>>read_token;
	pedFile<<read_token<<"\t";
   for(i=0;i<index;i++) {
		ped>>read_token;
   }
		ped>>read_token;
		pedFile<<read_token<<"\t";
		ped>>read_token;
		pedFile<<read_token<<"\n";
    ped.ignore(1000000,'\n');
num_row++;
  }
*numind=num_row;
//cout<<"out";
  ped.close();
  pedFile.close();
  free(read_token);
}

