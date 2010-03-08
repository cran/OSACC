void crankwholedata(unsigned long n,float v[],float w[],float s[],float u[],int *k)
{//sorted v[] and return rank w[] , unique value u[]
  unsigned long i, j=1,ji,jt,cutnum,mint;
	float rank;
	int m,numunaff=0,numaff=0;
	cutnum=n;
	for (i=1;i<=n;i++)
	{
		w[i]=v[i];
		if(s[i]==1) numunaff++;
		if(s[i]==2) numaff++;
		if(numaff>=20&&numunaff>=20)
		  {
		    mint=i;
		    if(cutnum>mint) cutnum=mint;
		  }

	}
	
	while (j < n) {
		
		 if(j<=cutnum) m=0;
		if (w[j+1] != w[j]) {
			m++;		  		 
			u[m]=w[j];
 
			w[j]=j;
			++j;
		
		
		} else {
		  m++;
		  u[m]=w[j];
		  

			for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
			rank=0.5*(j+jt-1);
			for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
			j=jt;
		}
	}
	if (j == n)
	  {m++;
 
	    u[m]=w[j];
	    w[n]=n;
	  }
	*k=m;
}
/* (C) Copr. 1986-92 Numerical Recipes Software "25Bi..<!. */
