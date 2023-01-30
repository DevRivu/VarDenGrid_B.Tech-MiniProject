#include<stdio.h> // standard I/O.
#include<ctype.h> // for isdigit() fn.
#include<stdlib.h> // for mrmory allocation.
#include<string.h> // to convert string to number. 
#include<math.h> // for pow() fn.

// 'num' is given number, median[] is set of k-means &
// fn call gives index of cluster to which num will be assign.
int Assign_cluster(float num,float median[],int k)
{
    float dist,min; int i,temp;
    for(i=1; i<=k; i++)
        { dist=(float)pow((num-median[i]),2);
        if(i==1) { min=dist; temp=i; }
        else if(i>1 && dist<min) { min=dist; temp=i; } }
    return temp;
}

typedef struct ll // Structure for linked list
{ 
    float data;
    struct ll *next;
}node;

node *addnode(node *ll,float data) // Fn. adds node(data) infront of linked list
{
    node *temp;
    temp=(node*)malloc(sizeof(node)); 
    temp->data=data;
    temp->next=ll; 
    ll=temp;
    return ll;
}

void print(node *ll) // print the linked list
{
    node *t=ll; // temperory pointer to traverse along linked list
    while(t!=NULL)
        { printf("%f ",t->data); t=t->next; }
    printf("\n");
}

void print_file(node *ll,FILE *f) // print linked list in file 
{
    node *t=ll; // temperory pointer to traverse along linked list
    while(t!=NULL)
        { fprintf(f,"%f ",t->data); t=t->next; }
    fprintf(f,"\n");
}

void BubbleSort(float a[],int n)
{
    int i,j; float t;
    for(i=1; i<n; ++i)
        { for(j=1; j<n; ++j)
            { if (a[j]>a[j+1])
            { t=a[j];
            a[j]=a[j+1];
            a[j+1]=t; } } }
}

float Med(node *ll)
{
    node *t=ll; int i,n=0; 
    while(t!=NULL)
        { t=t->next; n++; }
    t=ll;
    float a[n];
    for(i=1; i<=n; i++)
        { a[i]=t->data; t=t->next; }
    BubbleSort(a,n);
    if(n%2==0) { return ((float)((a[n/2]+a[n/2+1])/2)); }
    else { return (a[(n+1)/2]); }
}

float Max(float a[],int n) // Fn. returns maximum of array
{
    int i; float max=a[1];
    for(i=1; i<=n; i++)
        { if(a[i]>max) { max=a[i]; } }
    return max;
}

float Min(float a[],int n) // Fn. returns minimum of array
{
    int i; float min=a[1];
    for(i=1; i<=n; i++)
        { if(a[i]<min) { min=a[i]; } }
    return min;
}

float SSE(node *ll[],float mean[],int k) // Fn. returns sum_of_squared_error of final clustered data
{
    int i,n; float sse=0,temp=0; node *t; // temperory pointer to traverse along linked list
    for(i=1; i<=k; i++)
        { t=ll[i];
        while(t!=NULL)
            { temp=temp+pow((t->data-mean[i]),2); t=t->next; }
        sse=sse+temp; }
    return sse;
}

void main(int argc,char *argv[])
{
    int i=1,index,j=1,k,l,m=0,n=0,p,q,t=0,stop=1; float sse,p_med[10],c_med[10],max[20],min[20],a[800][20],b[800]; char c,*s,*s1,temp[10]; FILE *f; node *ll[10];
    s=(char*)malloc(sizeof(char)*(strlen(argv[1])+20));
    s1=(char*)malloc(sizeof(char)*(strlen(argv[1])+20));
    strcpy(s,argv[1]); strcat(s,".data");
    strcpy(s1,s); strcat(s1,"_cluster");
    f=fopen(s,"r"); // open file to read.
    fseek(f,0,2); // move to EOF.
    p=ftell(f); // 'p' gives number of character in file(f).
    fseek(f,0,0); // move to starting of file.
    for(q=0; q<p; q++) // read file and collect numerical data in 2-D array.
        { c=getc(f); // 'c' gives current character in file(f).
        if(isdigit(c) || c=='.') { temp[m++]=c; t=1; }
        else if(!isdigit(c) && c!='.' && t==1) {temp[m]='\0'; a[i][j]=atof(temp); m=0; j++; t=0; }
        if(c=='?') { a[i][j]=-1; j++; }
        if(c=='\n') { i++; l=j; j=1; }
    }

    for(q=1; q<l; q++) // Find max. and min. of all cols.
        {for(p=1; p<i; p++)
            { b[p]=a[p][q]; }
        max[q]=Max(b,i-1); 
        min[q]=Min(b,i-1); }

    for(q=1; q<l; q++) // Fill missing value
        {for(p=1; p<i; p++)
            { if(a[p][q]<0) { a[p][q]=max[q]+1; } } }

    printf("Number of clusters: "); scanf("%d",&k); // For 'k' clusters reads medians.
    printf("Median of clusters:\n");
    for(q=1; q<=k; q++)
        { printf("Median-%d) ",q); scanf("%f",&c_med[q]); p_med[q]=0; ll[q]=NULL; }

    for(q=1; q<l; q++) // assign all data to corr. cluster.
        {for(p=1; p<i; p++)
            { index=Assign_cluster(a[p][q],c_med,k); 
            ll[index]=addnode(ll[index],a[p][q]); } }

    for(q=1; q<=k; q++) // Find of new means of all clusters.
        { p_med[q]=c_med[q];
        c_med[q]=Med(ll[q]); 
        if(p_med[q]!=c_med[q]) { stop=0; } }

        while(stop!=1) // while iteration stops(current and prev. means are same) finish clustering. 
    {
        stop=1; // while stop=1 finish while loop.
        for(q=1; q<=k; q++) // refill all clusters.
            { ll[q]=NULL; }
        for(q=1; q<l; q++)
            {for(p=1; p<i; p++)
                { index=Assign_cluster(a[p][q],c_med,k);
                ll[index]=addnode(ll[index],a[p][q]); } }

        for(q=1; q<=k; q++) // Find of new means of all clusters.
            { p_med[q]=c_med[q];
            c_med[q]=Med(ll[q]);
            if(p_med[q]!=c_med[q]) { stop=0; } }
    }

    f=fopen(s1,"w"); // write clustered data in file f1.
    for(q=1; q<=k; q++)
        { printf("Cluster:%d) with median=%f\n",q,c_med[q]); print(ll[q]); print_file(ll[q],f); }

    sse=SSE(ll,c_med,k);
    printf("Sum of squared error=%f\n",sse); // prints sse of clustered data.

}