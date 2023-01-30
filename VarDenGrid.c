#include<stdio.h> // For I/O stream.
#include<stdlib.h> // for memory allocation.
#include<string.h> // to convert string to number.
#include<ctype.h> // for isdigit() fn.
#include<math.h> // For pow() Fn.
int grid=6;

float Max(float a[],int n) // Fn. returns maximum of array.
{
    int i; float max=a[1];
    for(i=1; i<=n; i++)
        { if(a[i]>max) { max=a[i]; } }
    return max;
}

float Min(float a[],int n) // Fn. returns minimum of array.
{
    int i; float min=a[1];
    for(i=1; i<=n; i++)
        { if(a[i]<min) { min=a[i]; } }
    return min;
}

float Scale(float num,float max,float min) // Scale the num. along column.
{
    float scale=(float)((num-min)/(max-min));
    return grid*scale;
}

void main(int argc,char *argv[])
{
    int i=1,index,row=1,col=1,j=1,k,l,p,q,p1,q1,m=0,n=0,t=0,include=1,cubeindex=0,size=0,psize,maxsize=0,mostdence,clusterIndex=0,SCALE; float tmp1,tmp2,max[20],min[20],a[800][20],b[800],d,min_inter,max_intra,DI,sum=0; char c,*s,*s1,temp[10]; FILE *f;
    s=(char*)malloc(sizeof(char)*(strlen(argv[1])+20));
    s1=(char*)malloc(sizeof(char)*(strlen(argv[1])+20));
    strcpy(s,argv[1]); strcat(s,".data");
    strcpy(s1,"dist_matrix_"); strcat(s1,s);
    f=fopen(s,"r"); // open file to read.
    fseek(f,0,2); // move to EOF.
    p=ftell(f); // 'p' gives number of character in file(f).
    fseek(f,0,0); // move to starting of file.
    for(q=0; q<p; q++) // read file and collect numerical data in 2-D array.
        { c=getc(f); // 'c' gives current character in file(f).
        if(isdigit(c) || c=='.') { temp[m++]=c; t=1; }
        else if(!isdigit(c) && c!='.' && t==1) { temp[m]='\0'; a[i][j]=atof(temp); m=0; j++; t=0; }
        if(c=='?') { a[i][j]=-1; j++; }
        if(c=='\n') { i++; l=j; j=1; }
    }

    for(q=1; q<l; q++) // Find max. and min. of all cols.
        {for(p=1; p<i; p++)
            { b[p]=a[p][q]; }
        max[q]=Max(b,i-1); 
        min[q]=Min(b,i-1); }

    for(q=1; q<l; q++) // Fill missing value.
        {for(p=1; p<i; p++)
            { if(a[p][q]<0) { a[p][q]=max[q]+1; } } }

    for(q=1; q<l; q++) // Find max. and  min. of filled data along cols.
        {for(p=1; p<i; p++)
            { b[p]=a[p][q]; }
        max[q]=Max(b,i-1);
        min[q]=Min(b,i-1); }

    fclose(f);

    float scale[i+2][l+2],dist[i+1][i+1],sqsum=0; 
    for(q=1; q<i; q++) /* Take scalled data bw 0-10 in a particular array with last col. ele(s)=-1 and zeroth col ele. will keep info about hupercube inedx in which point belongs. */
        { for(p=1; p<l; p++)
            { scale[q][p]=Scale(a[q][p],max[p],min[p]); }
          scale[q][p]=-1; scale[q][0]=0; }
   
    for(p=1; p<i; p++) // Fill distance matrix of scalled dataset.
        { for(q=1; q<i; q++)
            { for(m=1; m<l; m++)
             { sqsum=sqsum+pow(scale[p][m]-scale[q][m],2); } 
             dist[p][q]=sqrt(sqsum); sqsum=0; } }

    for(p=1; p<i; p++) // Assign index of hypercube to zeroth col. of scalled array.
        { if(scale[p][0]==0) 
            { scale[p][0]=++cubeindex;  
            for(q=p; q<i; q++)
            { for(m=1; m<l; m++)
                { if(scale[q][m]<=1) { tmp1=0.5; }
                  else { tmp1=scale[q][m]; }
                  if(scale[p][m]<=1) { tmp2=0.5; }
                  else { tmp2=scale[p][m]; }
                if(ceil(tmp1)!=ceil(tmp2)) { include=0; } }
            if(include==1) { scale[q][0]=scale[p][0]; size++; } include=1; } }
        if(maxsize<size) { maxsize=size; mostdence=cubeindex; } size=0; } // maxsize=sizeof most dence hypercube.

    f=fopen(s1,"w"); // Print distance matrix in file.
    for(p=1; p<i; p++)
        { for(q=1; q<i; q++)
            { fprintf(f,"%f ",dist[p][q]); }
        fprintf(f,"\n"); }
    fclose(f);

    for(p=1; p<i; p++) // Print scaled matrix. 
    { for(q=0; q<=l; q++)
      { printf("%f ",scale[p][q]); }
      printf("\n"); } 
	
    k=3;
    for(p1=0,q1=1; p1<k; p1++,q1++)
    {
    int hypercube[cubeindex+1][maxsize+1]; float nbd[i],nearnbd[i],epsilon,mean[k+1][l+1],sum=0;     
    for(p=1; p<=cubeindex;  p++) { hypercube[p][0]=1; }
    for(p=1; p<i; p++) { nbd[p]=0; nearnbd[p]=0; }

    for(q=1; q<i; q++) // Fill the index of number to corr. hypercube.
    { if(scale[q][l]==-1) { row=scale[q][0]; col=hypercube[row][0]; hypercube[row][col]=q; hypercube[row][0]++; } }

    for(q=1,maxsize=hypercube[1][0],mostdence=1; q<=cubeindex; q++)
    {
        if(maxsize<hypercube[q][0]) { maxsize=hypercube[q][0]; mostdence=q; }
    }
    printf("Most dence cube-%d with %d ele(s)\n",mostdence,maxsize-1);

    for(q=1,t=0,m=0; q<hypercube[mostdence][0]; q++) // Find epsilon=max. of nearest nbd.
    { for(p=1; p<hypercube[mostdence][0]; p++)
      { if(p!=q) { nbd[++t]=dist[hypercube[mostdence][q]][hypercube[mostdence][p]]; } } nearnbd[++m]=Min(nbd,t); t=0; } epsilon=Max(nearnbd,m);

    printf("Epsilon=%f\n",epsilon);

    int new[i]; size=0;
    for(p=0; p<i; p++) // Take new ele to include in cluster corr. to mostdence hypercube.
    { new[p]=0; }
    for(p=1; p<hypercube[mostdence][0]; p++)
    { for(q=1; q<i; q++)
      { if(scale[q][l]==-1 && dist[hypercube[mostdence][p]][q]<=epsilon)
        { new[++size]=q; scale[q][l]=q1; } } }
    psize=size; 

    while(size!=0) // Take new ele to include in cluster corr. to new eles. of clusters iteratively while no furthur ele. includes.
    { psize=size; size=0;
      for(p=1; p<=psize; p++)
      for(q=1; q<i; q++)
      { if(scale[q][l]==-1 && dist[new[p]][q]<=epsilon)
        { new[++size]=q; scale[q][l]=q1; } } } 

    }

    for(p=1,d=1000; p<i; p++)
    { if(scale[p][l]==-1) 
      { for(q=1; q<i; q++)
        { if(scale[q][l]!=-1)
          { if(d>dist[p][q]) { d=dist[p][q]; scale[p][l]=scale[q][l]; } } }
          d=1000; } }

    // Print cluster.
    int points[k+1]; float avg_cluster_dist[k+1];
    for(p=0; p<=k; p++)
    { points[p]=0; avg_cluster_dist[p]=0; }
    for(m=1; m<=k; m++)
    { 
    printf("cluster-%d:\n",m);
    for(p=1; p<i; p++)
    { if(scale[p][l]==m)
      { points[m]++;
        printf("point-%d: ",p);
        for(q=1; q<l; q++)
        { printf("%f ",scale[p][q]); }
      printf("\n"); } }
    printf("Number of points=%d\n",points[m]); }

    float mean[k+1][l],A[i],B[i],S[i],M,inter_cluster=0,intra_cluster=0; int ci=0,cj=0;
    for(q1=1; q1<=k; q1++)
    { for(p=1; p<l; p++) // calc. mean of clusters.
    { for(q=1; q<i; q++)
    { if(scale[q][l]==q1) { sum+=scale[q][p]; } }
    mean[q1][p]=(float)(sum/points[q1]); sum=0; } }

    printf("Mean of cluster: \n");
    for(q1=1; q1<=k; q1++)
    { printf("Cluster-%d: ",q1);
    for(p=1; p<l; p++)
    { printf("%f ",mean[q1][p]); }
    printf("\n"); }

    float mean_mean_dist[k+1][k+1],mean_point_dist[i+1]; sqsum=0;
    for(p=1; p<=k; p++) // Fill mean_mean_distance matrix of mean.
        { for(q=1; q<=k; q++)
            { for(m=1; m<l; m++)
             { sqsum=sqsum+pow(mean[p][m]-mean[q][m],2); }
             mean_mean_dist[p][q]=sqrt(sqsum); sqsum=0; printf("%f\n",mean_mean_dist[p][q]); } }

    for(p=1; p<i; p++)
    { for(q=1; q<l; q++)
      { SCALE=scale[p][l];
        sqsum=sqsum+pow(mean[SCALE][q]-scale[p][q],2); }
    mean_point_dist[p]=sqrt(sqsum); sqsum=0; 
    printf("%d=%f\n",p,mean_point_dist[p]); }

    for(p=1; p<i; p++)
    { SCALE=scale[p][l]; 
      avg_cluster_dist[SCALE]+=mean_point_dist[p]; }
    for(p=1; p<=k; p++)
    { avg_cluster_dist[p]/=points[p]; 
      printf("%d=%f\n",p,avg_cluster_dist[p]); } 

    float R[k+1][k+1],D[k+1],DB=0;
    for(p=1; p<=k; p++)
    { for(q=1; q<=k; q++)
      { R[p][q]=(float)((avg_cluster_dist[p]+avg_cluster_dist[q])/mean_mean_dist[p][q]); printf("%f\n",R[p][q]); } }
     
    for(p=1; p<=k; p++)
    { D[p]=0;
    for(q=1; q<=k; q++)
    { if(D[p]<R[p][q] && p!=q) { D[p]=R[p][q]; } } 
    DB+=D[p]; }
    DB/=k;

    printf("DB_index=%f\n",DB);

    max_intra=0; min_inter=1000;
    for(m=1; m<=k; m++)
    { for(p=1; p<i; p++)
      { if(scale[p][l]==k) 
        { for(q=1; q<i; q++)
          { if(scale[q][l]==k && p!=q)  
            { if(max_intra<dist[p][q]) { max_intra=dist[p][q]; } }
            if(scale[q][l]!=k) 
            { if(min_inter>dist[p][q]) { min_inter=dist[p][q]; } } } } } } 

    DI=(float)(min_inter/max_intra);
    printf("Dunn Index=%f\n",DI);

    for(p=1; p<i; p++)
    { for(q=1; q<i; q++)
      { if(scale[p][l]==scale[q][l]) { inter_cluster+=dist[p][q]; ci++; }
        else { intra_cluster+=dist[p][q]; cj++; } }
        A[p]=(float)(inter_cluster/(ci-1));
        B[p]=(float)(intra_cluster/cj);
        M=(B[p]>A[p])?B[p]:A[p];
        S[p]=(float)((B[p]-A[p])/M); 
        S[0]+=S[p]; 
	ci=0,cj=0,inter_cluster=0,intra_cluster=0; }
    S[0]=(float)(S[0]/i);

    printf("Silhouette Index=%f\n",S[0]);

    float Result=(float)((DI)/DB);
    printf("Result=%f\n",Result);

}