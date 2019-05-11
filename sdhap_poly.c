#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
//#include <atlas_enum.h>
#include <clapack.h>



/* DGESDD prototype */
extern void sgesdd_(char* jobz, int* m, int* n, float* a,
int* lda, float* s, float* u, int* ldu, float* vt, int* ldvt,
float* work, int* lwork, int* iwork, int* info);
/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, float* a, int lda );
#define Max_row 400000
#define Max_col 400000
#define Max_rank 19
#define Max_iter 1000
int max_array(int a[], int num_elements);
int  min_array(int arr[], int num_elements);
int  min_array_ind(int arr[], int num_elements);
int  min_array_indf(float arr[], int num_elements);
int MEC_calc_pos(int *list[], int *list2[], int *list_val[], int *list_val2[],int num_ct[],int num_ct2[], int *hapg1[], int pos, int K);
int MEC_calc(int *list[], int *list_val[], int num_ct[], int *hapg[],int max_row, int max_col, int K);
float obj_calc(float *conf_mat[],int *listc[],int readg[],int numc[], int max_row);
int low_rank_sdp(int row[],int col[],int val[],int max_row,int max_col,int basect,int blk,int num_ct[],int *list[],int *list_val[],int num_ct2[],int *list2[],int *list_val2[],int ploidy, int reord[],int *haplotype[]);
int low_rank_sdp2(int row[],int col[],int val[],int max_row,int max_col,int basect,int blk,int num_ct[],int *list[],int *list_val[],int num_ct2[],int *list2[],int *list_val2[]);
float func_ob_poly(float *Y[],float *conf_mat[],float *theta[], int *listc[],int numc[],int max_row,int r_rank);
float func_ob2(int X[],float *conf_mat[],int *listc[],int numc[],int max_row,int r_rank);
float func_ob_lag(float *X[],float *conf_mat[],int *listc[],int numc[],int max_row,int r_rank,float *theta[],float mu, int K);
float func_obS(int S[],float *conf_mat[],int *listc[],int numc[],int max_row);
float func_obSpos(int S[],float *conf_mat[],int *listc[],int numc[],int max_row,int pos);
int rank(float *X[],float th,int max_row,int r_rank);
float gaussian_random();
int cmpfunc (const void * a, const void * b);
int SWER(int *true_hap[],int *poly_list[],int reord[],int *hapg[], int K);
void swap (int *x, int *y);
/* void permute(int *a, int i, int n, int *hapg[], int *new_hap[], int K, int max_col, int *true_hap[], int reord[],int sw_best[]) ; */
void permute(int *a, int i, int n, int *new_hap[],int K, int max_col, int *true_hap[], int reord[], int flag,int pos) ;
int SWER2(int *true_hap[],int reord[],int *hapg[], int K, int len);
int compare (const void * a, const void * b);



int main(int argc, char *argv[])
{
	if (argc<4){ // ./exec_file frag_file output_hap number_haps
		printf("Error: The number of the arguments is %d which should be 3. \n",argc-1 );
		exit(-1);
	}

	FILE *fpt,*fpt1,*fpt2,*fpt5,*fpt7;
	int *row;
	int *col;
	int *val;
	int i,i2;
	int max_row, max_col;
	int iter,iter2;
	int cnt;
	clock_t t1;

	char tempo[10000];
	char qlty[10000];
	char str[10000];
	char *ru;
	int num_read;
	int num_col;
	int **snp_mat;
	int **snp_list;
	int *snp_ct;
	int MEC[10000];
	int num_base[10000];
	float SWER[10000];
	int block_len[10000];
	int block_rdlen[10000];
	int no_acc=0;

	int **arr_pos;
	arr_pos=malloc(10000*sizeof(int*));
	for (iter=0;iter<10000;iter++){
		arr_pos[iter]=malloc(10000*sizeof(int));
	}

	t1=clock();

	int K;
	fpt=fopen(argv[1],"r");
	fpt1=fopen(argv[1],"r");
	K=atoi(argv[3]);
	fscanf(fpt,"%d",&num_read);
	fscanf(fpt,"%d",&num_col);
	fscanf(fpt1,"%d",&num_read);
	fscanf(fpt1,"%d",&num_col);
	snp_mat=malloc(num_read*sizeof(int*));
	snp_list=malloc(num_read*sizeof(int*));
	snp_ct=malloc(num_read*sizeof(int));
	int **true_hap;
	int temp;
	true_hap=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		true_hap[iter]=malloc(num_col*sizeof(int));
	}
	int ct;
	int num_seg;
	int k;
	int pos;
	int d;
	int j;
	int *col_ct;
	int col_ct2;
	int *rd_array;
	int *rd_array2;
	int **list, **list2;
	int **list_val, **list_val2;
	int *num_ct, *num_ct2;
	int *read_ac;
	i=-1;i2=-1;
	num_ct=malloc(Max_row*sizeof(int));
	read_ac=malloc(Max_row*sizeof(int));
	list=malloc(Max_row*sizeof(int*));
	list_val=malloc(Max_row*sizeof(int*));
	for (iter=0;iter<Max_row;iter++){
		list[iter]=malloc(10000*sizeof(int));
		list_val[iter]=malloc(10000*sizeof(int));
	}
	num_ct2=malloc(Max_col*sizeof(int));
	list2=malloc(Max_col*sizeof(int*));
	list_val2=malloc(Max_col*sizeof(int*));
	for (iter=0;iter<Max_col;iter++){
		list2[iter]=malloc(10000*sizeof(int));
		list_val2[iter]=malloc(10000*sizeof(int));
	}
	rd_array=malloc(num_read*sizeof(int));
	rd_array2=malloc(num_read*sizeof(int));
	col_ct2=0;
	col_ct=malloc(num_col*sizeof(int));
	for(k=0;k<num_col;k++){
		col_ct[k]=0;
	}
	if(fpt==NULL)
	printf("\nError\n");
	else {
		while(i<num_read-1){ // while(!feof(fpt)){
			ct=-1;
			i++;
			fscanf(fpt,"%d",&num_seg);
			fscanf(fpt,"%s",tempo);
			fscanf(fpt1,"%d",&num_seg);
			fscanf(fpt1,"%s",tempo);
			for(k=0;k<num_seg;k++){
				fscanf(fpt,"%d",&pos);
				fscanf(fpt,"%s",str);
				d=strlen(str);
				for (j=0;j<d;j++){
					ct++;
				}
			}
			if((num_seg!=1)||(ct!=0)){
				i2++;
				snp_ct[i2]=ct+1;
				snp_mat[i2]=malloc(ct*sizeof(int*));
				snp_list[i2]=malloc(ct*sizeof(int*));
				ct=-1;
				for(k=0;k<num_seg;k++){
					fscanf(fpt1,"%d",&pos);
					fscanf(fpt1,"%s",str);
					d=strlen(str);
					for (j=0;j<d;j++){
						ct++;
						snp_list[i2][ct]=pos+j-1;
						col_ct[pos+j-1]++;
						col_ct2++;
						if(str[j]=='1'){
							snp_mat[i2][ct]=1;
						}
						if(str[j]=='2'){
							snp_mat[i2][ct]=2;
						}
						if(str[j]=='3'){
							snp_mat[i2][ct]=3;
						}
						if(str[j]=='4'){
							snp_mat[i2][ct]=4;
						}
					}
				}
				read_ac[i]=0;
				no_acc++;
			}
			else {
				for(k=0;k<num_seg;k++){
					fscanf(fpt1,"%d",&pos);
					fscanf(fpt1,"%s",str);
				}
				read_ac[i]=1;
			}
			fscanf(fpt,"%s",qlty);
			fscanf(fpt1,"%s",qlty);
		}
	}
	fclose(fpt);
	fclose(fpt1);
	int **haplotype;
	haplotype=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		haplotype[iter]=malloc(num_col*sizeof(int));
	}
	for(iter2=0;iter2<num_col;iter2++){
		for (iter=0;iter<K;iter++){
			haplotype[iter][iter2]=1;
		}
	}
	fpt=fopen(argv[1],"r");
	i=-1;i2=-1;
	fscanf(fpt,"%d",&num_read);
	fscanf(fpt,"%d",&num_col);
	int **snp_col;
	snp_col=malloc(num_col*sizeof(int*));
	for(k=0;k<num_col;k++){
		snp_col[k]=malloc(col_ct[k]*sizeof(int));
	}
	for(k=0;k<num_col;k++){
		col_ct[k]=0;
	}
	if(fpt==NULL)
	printf("\nError\n");
	else {
		while(i<num_read-1){
			ct=-1;
			i++;
			fscanf(fpt,"%d",&num_seg);
			fscanf(fpt,"%s",tempo);
			if(read_ac[i]==0){
				i2++;
				for(k=0;k<num_seg;k++){
					fscanf(fpt,"%d",&pos);
					fscanf(fpt,"%s",str);
					d=strlen(str);
					for (j=0;j<d;j++){
						ct++;
						snp_col[pos+j-1][col_ct[pos+j-1]]=i2;
						col_ct[pos+j-1]++;
					}
				}
			}
			else{
				for(k=0;k<num_seg;k++){
					fscanf(fpt,"%d",&pos);
					fscanf(fpt,"%s",str);
					d=strlen(str);
				}
			}
			fscanf(fpt,"%s",qlty);
		}
	}
	fclose(fpt);
	printf("Time to read data %f\n",(double)(clock()-t1)/CLOCKS_PER_SEC);
	num_read=no_acc;
	int *que_read;
	int *que_list;
	int *que_col;
	int iter_cnt;
	int que_ct;
	int *queue;
	int tu, iter3;

	t1=clock();
	row=malloc(10000000*sizeof(int));
	col=malloc(10000000*sizeof(int));
	val=malloc(10000000*sizeof(int));
	queue=malloc(10000000*sizeof(int));

	que_read=malloc(num_read*sizeof(int));
	que_list=malloc(num_read*sizeof(int));
	que_read=malloc(num_read*sizeof(int));
	que_col=malloc(num_col*sizeof(int));
	for(iter=0;iter<num_read;iter++){
		que_read[iter]=0;
	}
	for(iter=0;iter<num_read;iter++){
		que_list[iter]=0;
	}
	iter_cnt=1;
	que_ct=-1;

	int rd_cnt;
	int rd_cnt2;
	int *col_que;
	int colct;
	int *col_arr;
	int *col_arr2;
	int basect;
	int *reord;

	col_que=malloc(num_col*sizeof(int));
	reord=malloc(num_col*sizeof(int));
	col_arr=malloc(num_col*sizeof(int));
	col_arr2=malloc(num_col*sizeof(int));
	for(iter=0;iter<num_col;iter++){
		col_que[iter]=0;
	}
	rd_cnt2=0;
	for(iter=0;iter<num_read;iter++){
		rd_cnt=0;
		basect=0;
		if(que_read[iter]==0){
			que_read[iter]=iter_cnt;
			que_list[iter]=iter_cnt;
			rd_array[rd_cnt]=iter;
			rd_array2[iter]=rd_cnt;
			rd_cnt++;
			basect+=snp_ct[iter];
			for(iter2=0;iter2<snp_ct[iter];iter2++){
				for(iter3=0;iter3<col_ct[snp_list[iter][iter2]];iter3++){
					if(que_list[snp_col[snp_list[iter][iter2]][iter3]]==0){
						que_list[snp_col[snp_list[iter][iter2]][iter3]]=iter_cnt;
						que_ct++;
						queue[que_ct]=snp_col[snp_list[iter][iter2]][iter3];
					}
				}
			}
			while(que_ct!=-1){
				tu=queue[que_ct];
				que_read[tu]=iter_cnt;
				rd_array[rd_cnt]=tu;
				rd_array2[tu]=rd_cnt;
				rd_cnt++;
				basect+=snp_ct[tu];
				que_list[tu]=iter_cnt;
				que_ct--;
				for(iter2=0;iter2<snp_ct[tu];iter2++){
					for(iter3=0;iter3<col_ct[snp_list[tu][iter2]];iter3++){
						if(que_list[snp_col[snp_list[tu][iter2]][iter3]]==0){
							que_list[snp_col[snp_list[tu][iter2]][iter3]]=iter_cnt;
							que_ct++;
							queue[que_ct]=snp_col[snp_list[tu][iter2]][iter3];
						}
					}
				}
			}
			colct=0;
			for(iter2=0;iter2<rd_cnt;iter2++){
				for(iter3=0;iter3<snp_ct[rd_array[iter2]];iter3++){
					if(col_que[snp_list[rd_array[iter2]][iter3]]==0){
						col_que[snp_list[rd_array[iter2]][iter3]]=iter_cnt;
						colct++;
						col_arr[colct]=snp_list[rd_array[iter2]][iter3];
						col_arr2[snp_list[rd_array[iter2]][iter3]]=colct;
						reord[colct-1]=snp_list[rd_array[iter2]][iter3];
					}
				}
			}
			rd_cnt2+=rd_cnt;
			if(colct>1){
				i=0;
				for(iter2=0;iter2<rd_cnt;iter2++){
					for(iter3=0;iter3<snp_ct[rd_array[iter2]];iter3++){
						row[i]=iter2+1;
						col[i]=col_arr2[snp_list[rd_array[iter2]][iter3]];
						val[i]=snp_mat[rd_array[iter2]][iter3];
						i++;
					}
				}
				iter_cnt++;
				if(rd_cnt>1){
				}
				//printf("Base Count=%d\n",basect);
				num_base[iter_cnt-2]=basect;
				MEC[iter_cnt-2]=low_rank_sdp(row,col,val,rd_cnt,colct,basect,iter_cnt,num_ct,list,list_val,num_ct2,list2,list_val2,K,reord,haplotype);
				block_len[iter_cnt-2]=colct;
				block_rdlen[iter_cnt-2]=rd_cnt;
				qsort (reord,colct,sizeof(int), compare);
				for(iter=0;iter<colct;iter++){
					arr_pos[iter_cnt-2][iter]=reord[iter];
				}
			}
			else {
			}
		}
	}
	printf("Time to haplotype %f\n",(double)(clock()-t1)/CLOCKS_PER_SEC);
	int k2;
	int tot_MEC;
	int tot_base;
	float err_pc;
	tot_MEC=0;
	tot_base=0;
	fpt5=fopen(argv[2],"w");
	for(iter=0;iter<=iter_cnt-2;iter++){
		fprintf(fpt5,"Block %d\t Length of haplotype block %d\t Number of reads %d\t Total MEC %d\n",iter+1,block_len[iter],block_rdlen[iter],MEC[iter]);
		for(k=0;k<block_len[iter];k++){
			fprintf(fpt5,"%d\t",arr_pos[iter][k]);
			for (k2=0;k2<K;k2++){
				fprintf(fpt5,"%d\t",haplotype[k2][arr_pos[iter][k]]);
			}
			fprintf(fpt5,"\n");
		}
		tot_base+=num_base[iter];
		tot_MEC+=MEC[iter];
	}
	printf("Total MEC of all blocks= %d,Total num of bases=%d,error rate %f\n", tot_MEC, tot_base, (float)tot_MEC/(float)tot_base);
	fclose(fpt5);
	printf("Error Rate %f\n", (float)tot_MEC/(float)tot_base);
	/*  fpt5=fopen(argv[4],"w"); */
	/* fprintf(fpt5,"%f\n", (float)tot_MEC/(float)tot_base); */
	/*   fclose(fpt5); */
}




int low_rank_sdp(int row[],int col[],int val[],int max_row,int max_col,int basect,int blk,int num_ct[],int *list[],int *list_val[],int num_ct2[],int *list2[],int *list_val2[], int K, int reord[],int *haplotype[]){
	int i;
	int iter,iter2,iter3;
	int cnt;
	clock_t t1;
	int num_read;
	int num_col;
	for (iter=0;iter<max_row;iter++){
		num_ct[iter]=0;
	}
	for (iter=0;iter<basect;iter++){
		num_ct[row[iter]-1]++;
	}
	for (iter=0;iter<max_row;iter++){
		num_ct[iter]=0;
	}
	for (iter=0;iter<basect;iter++){
		num_ct[row[iter]-1]++;
		list[row[iter]-1][num_ct[row[iter]-1]-1]=col[iter]-1;
		list_val[row[iter]-1][num_ct[row[iter]-1]-1]=val[iter];
	}
	for (iter=0;iter<max_col;iter++){
		num_ct2[iter]=0;
	}
	for (iter=0;iter<basect;iter++){
		num_ct2[col[iter]-1]++;
	}
	for (iter=0;iter<max_col;iter++){
		num_ct2[iter]=0;
	}
	for (iter=0;iter<basect;iter++){
		num_ct2[col[iter]-1]++;
		list2[col[iter]-1][num_ct2[col[iter]-1]-1]=row[iter]-1;
		if(list2[col[iter]-1][num_ct2[col[iter]-1]-1]<0){
		}
		list_val2[col[iter]-1][num_ct2[col[iter]-1]-1]=val[iter];
	}
	for (iter=0;iter<max_col;iter++){
		// printf("%d %d \n",iter,num_ct2[iter]);
	}
	int **hapg;
	int **hapg5;
	int **hapg1;
	float temp;
	int k;
	hapg=malloc(K*sizeof(int*));
	hapg5=malloc(K*sizeof(int*));
	hapg1=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		hapg[iter]=malloc(max_col*sizeof(int));
		hapg5[iter]=malloc(max_col*sizeof(int));
		hapg1[iter]=malloc(max_col*sizeof(int));
	}

	/* Conflict matrix formation */
	int iter4, iter5;
	int *listn;
	int num,c1,c2;
	int **listc;
	float **conf_mat;
	int *numc;
	float u2;
	float conf_val=0.0;
	float conf_val2=0.0;
	int ltr;
	float **theta;
	float **theta2;
	int it;
	listc=malloc(max_row*sizeof(int*));
	conf_mat=malloc(max_row*sizeof(float*));
	theta=malloc(max_row*sizeof(float*));
	theta2=malloc(max_row*sizeof(float*));
	numc=malloc(max_row*sizeof(int));
	listn=malloc(max_row*sizeof(int));
	FILE *fpt3,*fpt6,*fpt8;
	for (iter=0;iter<max_row;iter++){
		for (iter5=0;iter5<max_row;iter5++){
			listn[iter5]=0;
		}
		num=0;
		for (iter2=0;iter2<num_ct[iter];iter2++){
			for (iter3=0;iter3<num_ct2[list[iter][iter2]];iter3++){
				if(list2[list[iter][iter2]][iter3]!=iter){
					if(listn[list2[list[iter][iter2]][iter3]]!=1){
						listn[list2[list[iter][iter2]][iter3]]=1;
						num++;

					}
				}
			}
		}
		numc[iter]=num;
		listc[iter]=malloc(num*sizeof(int));
		conf_mat[iter]=malloc(num*sizeof(float));
		theta[iter]=malloc(num*sizeof(float));
		theta2[iter]=malloc(num*sizeof(float));
		num=0;
		for (iter5=0;iter5<max_row;iter5++){
			listn[iter5]=0;
		}
		for (iter2=0;iter2<num_ct[iter];iter2++){
			for (iter3=0;iter3<num_ct2[list[iter][iter2]];iter3++){
				if(list2[list[iter][iter2]][iter3]!=iter){
					if(listn[list2[list[iter][iter2]][iter3]]!=1){
						listn[list2[list[iter][iter2]][iter3]]=1;
						c1=0;
						c2=0;
						for (iter4=0;iter4<num_ct[iter];iter4++){
							for (iter5=0;iter5<num_ct[list2[list[iter][iter2]][iter3]];iter5++){
								if(list[iter][iter4]==list[list2[list[iter][iter2]][iter3]][iter5]){
									c1++;
									if(list_val[iter][iter4]==list_val[list2[list[iter][iter2]][iter3]][iter5]){
										c2++;
									}
								}
							}
						}
						if(c1!=0){
							listc[iter][num]=list2[list[iter][iter2]][iter3];
							//u2=(2.0*(float)c2-(float)c1);
							u2=(2.0*(float)c2-(float)c1)/(float)c1;
							//u2=((float)c2-(float)c1)/(float)c1;
							u2=-u2;
							num++;
							conf_val+=fabs(u2);
							if(u2<0){
								conf_val2+=fabs(u2)*(-1.0/(1.0-(float)K));
							}
							else{
								conf_val2+=fabs(u2);
							}
							conf_mat[iter][num-1]=u2;
							if(u2>0){
								theta[iter][num-1]=fabs(u2)*(1.0-1.0/((float)K));
								theta2[iter][num-1]=theta[iter][num-1];
							}
							else {
								theta[iter][num-1]=0.0;
								theta2[iter][num-1]=0.0;
							}
						}
					}
				}
			}
		}
		numc[iter]=num;
	}

	int MEC;
	if(conf_val>0){
		/* Low rank solution */
		int r_rank;
		int iteration;
		int stop_iter, stop_iter2,stop_iter3,step_iter;
		float **X,**X_new, **Grad, **gr;
		float *X_n;
		float l,l2,l3;
		float eeta,eeta2;
		float obj, obj2;
		float *object;
		float *object2;
		float *conv_as;
		int k2,k3;
		int i2,j2,i3;
		int iter_rank,iter_rank2;
		float gr_norm;
		int lr;
		float en;
		float r_th=0.01*sqrt((float)max_row);
		int fl;
		int *S;
		float mu=1.0;
		float sc=0.8;
		float th=1.0;
		float cons[1000];
		float cons_gr[1000];
		float r_av;
		if(r_th<1.0){
			r_th=1.0;
		}
		//r_th=0.01;
		// printf("Into Low rank\n");
		object=malloc(Max_iter*sizeof(float));
		conv_as=malloc(Max_iter*sizeof(float));
		object2=malloc(Max_iter*sizeof(float));
		S=malloc(max_row*sizeof(int));
		X=malloc(max_row*sizeof(float*));
		for (k2=0;k2<max_row;k2++){
			X[k2]=malloc(Max_rank*sizeof(float));
		}
		X_n=malloc(max_row*sizeof(float));
		X_new=malloc(max_row*sizeof(float*));
		for (k2=0;k2<max_row;k2++){
			X_new[k2]=malloc(Max_rank*sizeof(float));
		}
		Grad=malloc(max_row*sizeof(float*));
		for (k2=0;k2<max_row;k2++){
			Grad[k2]=malloc(Max_rank*sizeof(float));
		}
		gr=malloc(max_row*sizeof(float*));
		for (k2=0;k2<max_row;k2++){
			gr[k2]=malloc(Max_rank*sizeof(float));
		}
		srand((unsigned int)time(NULL));
		for (k2=0;k2<max_row;k2++){
			S[k2]=(int)4.0*((float)(rand())/RAND_MAX);
			// printf("%d\t",S[k2]);
		}
		int stop;
		float val_max;
		float valalt[1000];
		float m2[K];
		float r;
		int out_it;
		int in_it, ctr, stop2, ty;
		float nr,sf;
		int cnt2;
		FILE *fpt5;
		float ratio;
		iter=-1;
		stop=0;
		val_max=func_obS(S,conf_mat,listc,numc,max_row);
		while(stop==0){
			iter++;
			for(k2=0;k2<max_row;k2++){
				for(k=0;k<K;k++){
					S[k2]=k;
					m2[k]=func_obSpos(S,conf_mat,listc,numc,max_row,k2);
					m2[k]=m2[k];
				}
				S[k2]=min_array_indf(m2,K);
			}
			valalt[iter]=func_obS(S,conf_mat,listc,numc,max_row);
			if(iter>1){
				if(valalt[iter]==valalt[iter-1]){
					stop=1;
				}
			}
		}
		r_rank=K;
		srand((unsigned int)time(NULL));
		for (k2=0;k2<max_row;k2++){
			for (k3=0;k3<r_rank;k3++){
				X[k2][k3]=-1.0/(float)K;
			}
			X[k2][S[k2]]+=1.0;
			l=1.0-1.0/(float)K;
			for (k3=0;k3<r_rank;k3++){
				X[k2][k3]=X[k2][k3]/sqrt(l);
			}
		}

		/* Gradient computation */
		object[0]=func_ob_poly(X,conf_mat,theta,listc,numc,max_row,r_rank);
		obj2=object[0];
		iteration=1;
		iter_rank2=r_rank;
		r_rank=iter_rank2;
		stop_iter2=0;
		stop_iter3=0;
		stop=0;
		out_it=0;
		in_it=0;
		ctr=0;
		iter=0;
		while(stop==0){
			ctr++;
			out_it++;
			iter++;
			in_it++;
			stop2=0;
			ty=0;
			while(stop2==0){
				ty++;
				for(k2=0;k2<max_row;k2++){
					for(j2=0;j2<r_rank;j2++){
						X_n[j2]=0.0;
					}
					for(i3=0;i3<numc[k2];i3++){
						for(j2=0;j2<r_rank;j2++){
							X_n[j2]-=(conf_mat[k2][i3]-theta[k2][i3])*X[listc[k2][i3]][j2];
						}
					}
					nr=0.0;
					for(j2=0;j2<r_rank;j2++){
						nr+=pow(X_n[j2],2);
					}
					nr=sqrt(nr);
					if(nr!=0.0){
						for(i3=0;i3<r_rank;i3++){
							X[k2][i3]=X_n[i3]/nr;
						}
					}
				}
				object[out_it]=func_ob_poly(X,conf_mat,theta,listc,numc,max_row,r_rank);
				if(iter>=5){
					iter=0;
					lr=rank(X,r_th,max_row,r_rank);
					//  printf("Rank %d\n",lr);
					if(lr<Max_rank){
						if(lr==iter_rank2){
							iter=0;
							iter_rank2++;
							r_rank++;
							for (k2=0;k2<max_row;k2++){
								l=0.0;
								for (k3=0;k3<r_rank;k3++){
									X[k2][r_rank-1]=0.01*(((float)(rand())/RAND_MAX)-0.5);
									l+=pow(X[k2][k3],2);
								}
								for (k3=0;k3<r_rank;k3++){
									X[k2][k3]=X[k2][k3]/sqrt(l);
								}
							}
						}
					}
				}
				if(ty==1){
					stop2=1;
				}
			}
			object[out_it]=func_ob_poly(X,conf_mat,theta,listc,numc,max_row,r_rank);
			//  printf("obj2 temp = %f\n",object[out_it]);
			/* Update mu, theta */
			//printf("Update\n");
			if(in_it>=1){
				//	printf("theta update\n");
				in_it=0;
				for(i2=0;i2<max_row;i2++){
					// printf("i2=%d\n",i2);
					for(i3=0;i3<numc[i2];i3++){
						//printf("hey2\n");
						r=0.0;
						for(j2=0;j2<r_rank;j2++){
							r+=-X[i2][j2]*X[listc[i2][i3]][j2];
						}
						r+=-(1.0/((float)K-1.0));
						theta2[i2][i3]=r;
						sf=r*fabs(conf_mat[i2][i3]);
						//printf("hey2\n");
						if(theta[i2][i3]>0.001){
							if(sf>0.0){
								//	theta[i2][i3]=(pow(2,r))*theta[i2][i3];
								theta[i2][i3]=(pow(2,r/(float)(iteration)))*theta[i2][i3];
							}
							else
							{
								// theta[i2][i3]=(pow(2,r))*theta[i2][i3];
								theta[i2][i3]=(pow(2,r/(float)(iteration)))*theta[i2][i3];
							}
						}
						else {
							theta[i2][i3]=theta[i2][i3]+mu*r/(float)(iteration);
						}
						if(theta[i2][i3]>5.0){
							theta[i2][i3]=5.0;
						}
						if(theta[i2][i3]<0.0){
							theta[i2][i3]=0.0;
						}
						if(r>0.0){
							theta2[i2][i3]=r;
						}
						else {
							theta2[i2][i3]=0.0;
						}
					}
				}
				cnt++;
				cnt2++;
			}
			if(out_it>50){
				r_av=0.0;
				for(it=0;it<51;it++){
					r_av+=object[out_it-it];
				}
				r_av=r_av/50.0;
				conv_as[out_it]=r_av;
				if(out_it>51){
					if(fabs(conv_as[out_it])>fabs(conv_as[out_it-1])){
						ratio=fabs(conv_as[out_it])/fabs(conv_as[out_it-1]);
						//   printf("ratio %f\n",ratio);
						if((ratio<1.00001)&&(out_it>100)&&(lr<r_rank)){
							stop=1;
						}
					}
				}
			}
			object[out_it]=func_ob_poly(X,conf_mat,theta,listc,numc,max_row,r_rank);
			// printf("obj2 = %f\n",object[out_it]);
			/* if((out_it==500)||(mu<0.01)||(cons[out_it]>-0.05)){ */
			/* 	stop=1; */
			/* } */
			if(out_it==Max_iter){
				stop=1;
			}
			/* all_obj(out_it,1:iteration-1)=object(1:iteration-1); */
			/* object_un(out_it)=func_ob(X,conf_mat,num_ar,numl); */
			stop_iter2=0;
			stop_iter3=0;
			iteration=1;
			//  printf("hi Cons %f out %d stop %d mu %f\n",cons[out_it],out_it,stop,mu);
		}
		// printf("col %d\n",max_col);
		/* FILE *fpt; */
		/* fpt=fopen("/home/polaris/shree/haplotyping/objective.txt","w"); */
		/* for(i2=0;i2<Max_iter;i2++){ */
		/*   fprintf(fpt,"%f\t",object[i2]); */
		/* } */
		/* //fprintf("\n"); */
		/* fclose(fpt); */
		//  printf("Time for conflict %f\n",(double)(clock()-t1)/CLOCKS_PER_SEC);
		t1=clock();
		float temp1;
		float **ran_proj;
		int *readg;
		int *readg_opt;
		int numran;
		float ran_max=50000.0;
		float ranval;
		float ran_run;
		int it;
		float *ran_arr;
		float *ran_mean;
		//fpt=fopen("/home/polaris/shree/haplotyping/rand_proj.txt","w");
		numran=max_row;
		readg=malloc(max_row*sizeof(int));
		readg_opt=malloc(max_row*sizeof(int));
		// printf("hi1\n");
		ran_proj=malloc(r_rank*sizeof(float*));
		for (k3=0;k3<r_rank;k3++){
			ran_proj[k3]=malloc(K*sizeof(float));
		}
		// printf("hi\n");
		/* if(max_row<500){ */
		/*   it=500*K;} */
		/* else{  */
		/*   it=(int)(max_row*K*log10(max_row)); */
		/* } */
		/* if(it>10000){ */
		/*   it=10000; */
		/* } */
		it=(int)(500*K*log10(max_row));
		// it=20;
		//printf("iter %d\n",it);
		ran_arr=malloc(it*sizeof(float));
		ran_mean=malloc(it*sizeof(float));
		ran_run=0.0;
		for(iteration=0;iteration<it;iteration++){
			// printf("iteration %d\n",iteration);
			for (k3=0;k3<K;k3++){
				l3=0.0;
				for (iter2=0;iter2<r_rank;iter2++){
					ran_proj[k3][iter2]=gaussian_random();
					l3+=ran_proj[k3][iter2];
				}
				for (iter2=0;iter2<r_rank;iter2++){
					ran_proj[k3][iter2]=ran_proj[k3][iter2]/sqrt(l3);
				}
			}
			//printf("hi\n");
			for (iter=0;iter<max_row;iter++){
				for (k3=0;k3<K;k3++){
					temp1=0.0;
					for (iter2=0;iter2<r_rank;iter2++){
						temp1+=X[iter][iter2]*ran_proj[k3][iter2];
					}
					m2[k3]=-temp1;
				}
				readg[iter]=min_array_indf(m2,K);
			}
			// printf("hi\n");
			ranval=func_obS(readg,conf_mat,listc,numc,max_row);
			ran_arr[iteration]=ranval;
			ran_run+=ranval;
			ran_mean[iteration]=ran_run/(float)(iteration+1);
			// fprintf(fpt,"%f\n",ran_arr[iteration]);
			//  printf("iteration  %d ranval %f ran_max %f \n",iteration,ranval,ran_max);
			if(ranval<ran_max){
				ran_max=ranval;
				for(iter3=0;iter3<numran;iter3++){
					readg_opt[iter3]=readg[iter3];
					S[iter3]=readg_opt[iter3];
				}
			}
		}
		// fclose(fpt);
		//   printf("Time for random projections  %f ran_max %f\n",(double)(clock()-t1)/CLOCKS_PER_SEC,ran_max);


		/* Greedy improvement */
		iter=-1;
		stop=0;
		val_max=func_obS(S,conf_mat,listc,numc,max_row);
		//printf("value %f\n",val_max);
		for(k2=0;k2<max_row;k2++){
			for(k=0;k<K;k++){
				X[k2][k]=-1.0/(float)K;
			}
			X[k2][readg_opt[k2]]+=1.0;
			l3=0.0;
			for(k=0;k<K;k++){
				l3+=pow(X[k2][k],2);
			}
			l3=sqrt(l3);
			for(k=0;k<K;k++){
				X[k2][k]=X[k2][k]/l3;
			}
		}
		lr=rank(X,r_th,max_row,K);
		//printf("Final projected value %f\n",func_ob_poly(X,conf_mat,theta,listc,numc,max_row,r_rank));
		// sleep(5);
		while(stop==0){
			iter++;
			for(k2=0;k2<max_row;k2++){
				for(k=0;k<K;k++){
					S[k2]=k;
					m2[k]=func_obSpos(S,conf_mat,listc,numc,max_row,k2);
				}
				/* for(k3=0;k3<K;k3++){ */
				/* 	printf("m[%d]=%f\t",k3,m2[k3]); */
				/* } */
				S[k2]=min_array_indf(m2,K);
				readg_opt[k2]=S[k2];
				/* printf("%d\n",S[k2]); */
			}
			valalt[iter]=func_obS(S,conf_mat,listc,numc,max_row);
			//   printf("value %f\n",valalt[iter]);
			if(iter>1){
				if(valalt[iter]==valalt[iter-1]){
					stop=1;
				}
			}
		}

		// printf("Cluster\n");
		/* for(iter=0;iter<max_row;iter++){ */
		/*   printf("%d\n",readg_opt[iter]); */
		/* } */
		float obj_val=0.0;
		for(iter=0;iter<max_row;iter++){
			for(iter2=0;iter2<numc[iter];iter2++){
				if(readg_opt[iter]==readg_opt[listc[iter][iter2]]){
					obj_val+=conf_mat[iter][iter2];
				}
				else{
					obj_val-=conf_mat[iter][iter2];
				}
			}
		}
		// printf("hello\n");
		// printf("Total %f\n",obj_val);
		float objective[10000];
		float s;
		iter=-1;
		stop=0;
		iter=-1;
		stop=0;
		// t1=clock();
		int **hap_ct;
		int m3[4];
		hap_ct=malloc(4*K*sizeof(int*));
		for(iter=0;iter<4*K;iter++){
			hap_ct[iter]=malloc(max_col*sizeof(int));
		}
		for(iter=0;iter<4*K;iter++){
			for(iter2=0;iter2<max_col;iter2++){
				hap_ct[iter][iter2]=0;
			}
		}
		for(iter=0;iter<max_row;iter++){
			for(iter2=0;iter2<num_ct[iter];iter2++){
				hap_ct[readg_opt[iter]*4+list_val[iter][iter2]-1][list[iter][iter2]]++;
			}
		}
		for(iter=0;iter<max_col;iter++){
			for(iter2=0;iter2<K;iter2++){
				for(iter3=0;iter3<4;iter3++){
					m3[iter3]=-hap_ct[iter2*4+iter3][iter];
				}
				hapg[iter2][iter]=min_array_ind(m3,4)+1;
				hapg1[iter2][iter]=hapg[iter2][iter];
			}
		}
		//  printf("%d\n",MEC_calc(list,list_val,num_ct,hapg,max_row,max_col,K));

		/* Greedy update of columns */
		int MEC2[100];
		int iter5;
		iter5=0;
		stop=0;
		while(stop==0){
			iter5++;
			for(iter=0;iter<max_col;iter++){
				for(iter2=0;iter2<K;iter2++){
					for(iter3=0;iter3<4;iter3++){
						hapg1[iter2][iter]=iter3+1;
						m3[iter3]=MEC_calc_pos(list,list2,list_val,list_val2,num_ct,num_ct2,hapg1,iter,K);
					}
					hapg[iter2][iter]=min_array_ind(m3,4)+1;
					hapg1[iter2][iter]=hapg[iter2][iter];
				}
			}
			MEC2[iter5]=MEC_calc(list,list_val,num_ct,hapg,max_row,max_col,K);
			//printf("%d %d\n",iter5,MEC2[iter5]);
			if(iter5>0){
				if(MEC2[iter5]==MEC2[iter5-1]){
					stop=1;
				}
			}
		}
		int k1;
		for(k2=0;k2<max_col;k2++){
			for(k1=0;k1<K;k1++){
				haplotype[k1][reord[k2]]=hapg[k1][k2];
			}
		}
		MEC2[0]=MEC_calc(list,list_val,num_ct,hapg,max_row,max_col,K);
		//printf("MEC=%d\n",MEC2[0]);
		MEC=MEC2[0];
	}
	else {
		MEC=0;
	}
	return(MEC);
}


//float SWER(int *poly_list[],int reord[],int *hapg[], int K){
int SWER2(int *true_hap[],int reord[],int *hapg[], int K,int max_col){
	int i;
	int *a;
	a=malloc(K*sizeof(int));
	for (i=0;i<K;i++){
		a[i]=i;
	}
	int iter,iter3;
	int j;
	int sw_best;
	int swch;
	int sw_ct;
	int **new_hap;
	int flag;
	new_hap=malloc(K*sizeof(int*));
	for (iter=0;iter<K;iter++){
		new_hap[iter]=malloc(max_col*sizeof(int));
	}
	for (iter=0;iter<max_col;iter++){
		for(iter3=0;iter3<K;iter3++){
			new_hap[iter3][iter]=hapg[iter3][iter];
		}
	}
	swch=0;
	for (iter=0;iter<max_col;iter++){
		sw_ct=0;
		for(iter3=0;iter3<K;iter3++){
			if(new_hap[iter3][iter]!=true_hap[iter3][reord[iter]]){
				sw_ct++;
			}
		}
		if(sw_ct!=0){
			if(iter!=0){
				swch++;
			}
			flag=0;
			for (i=0;i<K;i++){
				a[i]=i;
			}
			permute(a,0,K-1,new_hap,K,max_col,true_hap,reord,flag,iter);
		}
	}
	return swch;
	// sleep(10);
}


/* Function to swap values at two pointers */
void swap (int *x, int *y)
{
	int temp;
	temp = *x;
	*x = *y;
	*y = temp;
}

/* Function to print permutations of string
This function takes three parameters:
1. String
2. Starting index of the string
3. Ending index of the string. */
void permute(int *a, int i, int n, int *new_hap[],int K, int max_col, int *true_hap[], int reord[], int flag,int pos)
{
	int j;
	int i2;
	int iter;
	int iter3;
	int swch;
	int sw_ct;
	int temp[K];
	//printf("pos flag %d %d\n",pos,flag);
	/* sleep(5); */
	if(flag==0){
		if (i == n){
			/* for(i2=0;i2<=n;i2++){ */
			/* 	printf("%d\t", a[i2]); */
			/* } */
			/* printf("\n"); */
			/* for(i2=0;i2<K;i2++){ */
			/* 	  new_hap[i2][pos]=hapg[a[i2]][pos]; */
			/* 	} */
			/* printf("hello\n"); */
			/* for(iter3=0;iter3<K;iter3++){ */
			/*   for(iter=0;iter<max_col;iter++){ */
			/*     printf("%d\t",hapg[iter3][iter]); */
			/*   } */
			/*   printf("\n"); */
			/* } */
			/* for(iter3=0;iter3<K;iter3++){ */
			/*   for(iter=0;iter<max_col;iter++){ */
			/*     printf("%d\t",new_hap[iter3][iter]); */
			/*   } */
			/*   printf("\n"); */
			/* } */
			/* for(iter3=0;iter3<K;iter3++){ */
			/*   for(iter=0;iter<max_col;iter++){ */
			/*     printf("%d\t",true_hap[iter3][iter]); */
			/*   } */
			/*   printf("\n"); */
			/* } */
			sw_ct=0;
			for(iter3=0;iter3<K;iter3++){
				if(new_hap[a[iter3]][pos]!=true_hap[iter3][reord[pos]]){
					sw_ct++;
				}
			}
			if(sw_ct==0){
				flag=1;
				for(iter3=pos;iter3<max_col;iter3++){
					for(i2=0;i2<K;i2++){
						temp[i2]=new_hap[a[i2]][iter3];
					}
					for(i2=0;i2<K;i2++){
						new_hap[i2][iter3]=temp[i2];
					}
				}
			}
		}
		else
		{
			for (j = i; j <= n; j++)
			{
				swap((a+i), (a+j));
				permute(a,i+1,n,new_hap,K,max_col,true_hap,reord,flag,pos);
				swap((a+i), (a+j)); //backtrack
			}
		}
	}
}


int  max_array(int arr[], int num_elements)
{
	int i;
	int max=-32000;
	for (i=0; i<num_elements; i++)
	{
		if (arr[i]>max)
		{
			max=arr[i];
		}
	}
	return(max);
}

int  min_array(int arr[], int num_elements)
{
	int i;
	int min=32000;
	for (i=0; i<num_elements; i++)
	{
		if (arr[i]<min)
		{
			min=arr[i];
		}
	}
	return(min);
}

int  min_array_ind(int arr[], int num_elements)
{
	int i;
	int min=32000;
	int min_ind;
	for (i=0; i<num_elements; i++)
	{
		if (arr[i]<min)
		{
			min=arr[i];
			min_ind=i;
		}
	}
	return(min_ind);
}

int  min_array_indf(float arr[], int num_elements)
{
	int i;
	float min=32000.0;
	int min_ind;
	for (i=0; i<num_elements; i++)
	{
		if (arr[i]<min)
		{
			min=arr[i];
			min_ind=i;
		}
	}
	return(min_ind);
}

int MEC_calc_pos(int *list[], int *list2[], int *list_val[], int *list_val2[],int num_ct[],int num_ct2[], int *hapg1[], int pos, int K)
{
	int iter,iter2,iter3,k;
	int er2;
	int er[K];
	iter=0;
	iter2=0;
	er2=0;
	//printf("%d\n",pos);
	// printf("%d\n",num_ct2[pos]);
	for (iter=0;iter<num_ct2[pos];iter++){
		// printf("%d\n",list2[pos][iter]);
		//  er[0]=0;er[1]=0;
		for(iter3=0;iter3<K;iter3++){
			er[iter3]=0;
			for(iter2=0;iter2<num_ct[list2[pos][iter]];iter2++){
				//  printf("%d %d %d\n",iter,iter2,list_val[iter][iter2]);
				if(list_val[list2[pos][iter]][iter2]==hapg1[iter3][list[list2[pos][iter]][iter2]]){
					er[iter3]++;
				}
			}
		}
		//  printf("\n");
		// printf("%d %d\n",er[0],er[1]);
		er2+=num_ct[list2[pos][iter]]-max_array(er,K);
		// printf("%d\t",er2);
	}
	//printf("%d\t",er2);
	return(er2);
}

int MEC_calc(int *list[], int *list_val[], int num_ct[], int *hapg[],int max_row, int max_col,int K)
{
	int iter,iter2,iter3;
	int er2;
	int er[K];
	// printf("hi\n");
	er2=0;
	for(iter=0;iter<max_row;iter++){
		// printf("%d %d\t",iter,num_ct[iter]);
		for(iter3=0;iter3<K;iter3++){
			er[iter3]=0;
			for(iter2=0;iter2<num_ct[iter];iter2++){
				//  printf("%d %d %d\n",iter,iter2,list_val[iter][iter2]);
				if(list_val[iter][iter2]==hapg[iter3][list[iter][iter2]]){
					er[iter3]++;
				}
			}
		}
		er2+=num_ct[iter]-max_array(er,K);
	}
	return(er2);
}


float obj_calc(float *conf_mat[],int *listc[],int readg[],int numc[],int max_row){
	float s;
	int k2, iter2;
	s=0.0;
	for(k2=0;k2<max_row;k2++){
		for(iter2=0;iter2<numc[k2];iter2++){
			s+=conf_mat[k2][iter2]*readg[listc[k2][iter2]]*readg[k2];
		}
	}
	return(s);
}


float func_ob_poly(float *X[],float *conf_mat[],float *theta[], int *listc[],int numc[],int max_row,int r_rank){
	int i2,i3,k3;
	float obj;
	float l,l2,l3;
	obj=0.0;
	for(i2=0;i2<max_row;i2++){
		for(i3=0;i3<numc[i2];i3++){
			l=0.0;
			for (k3=0;k3<r_rank;k3++){
				l+=pow(X[i2][k3],2);
			}
			l=sqrt(l);
			l2=0.0;
			for (k3=0;k3<r_rank;k3++){
				l2+=pow(X[listc[i2][i3]][k3],2);
			}
			l2=sqrt(l2);
			l3=0.0;
			for (k3=0;k3<r_rank;k3++){
				l3+=X[listc[i2][i3]][k3]*X[i2][k3];
			}
			//printf("%f %f %f\n",l,l2,l3);
			obj+=(conf_mat[i2][i3]-theta[i2][i3])*(l3)/(l*l2);
			//printf("iter %d %d obj = %f\n",i2,i3,obj);
		}
	}
	// printf("%f\n",obj);
	if(isnan(obj)){
		printf("obj nan func %f %f %f %d\n",l,l2,l3,max_row);
		//sleep(10);
	}
	return(obj);
}


float func_ob_lag(float *X[],float *conf_mat[],int *listc[],int numc[],int max_row,int r_rank,float *theta[],float mu, int K){
	int i2,i3,k3;
	float obj;
	float l,l2,l3;
	float r;
	obj=0.0;
	for(i2=0;i2<max_row;i2++){
		for(i3=0;i3<numc[i2];i3++){
			l3=0.0;
			for (k3=0;k3<r_rank;k3++){
				l3+=X[listc[i2][i3]][k3]*X[i2][k3];
			}
			r=l3+(1.0/((float)K-1.0));
			if(r>theta[i2][i3]*mu){
				obj+=conf_mat[i2][i3]*r-(pow(theta[i2][i3],2)*mu/2.0);
			}
			else{
				obj+=conf_mat[i2][i3]*r-theta[i2][i3]*r+(1.0/(2.0*mu))*pow(fabs(r),2);
			}
		}
	}
	// printf("%f\n",obj);
	if(isnan(obj)){
		printf("obj nan func %f %f %f %d\n",l,l2,l3,max_row);
		//sleep(10);
	}
	return(obj);
}


float func_obS(int S[],float *conf_mat[],int *listc[],int numc[],int max_row){
	int i2,i3,k3;
	float obj;
	float l,l2,l3;
	obj=0.0;
	for(i2=0;i2<max_row;i2++){
		for(i3=0;i3<numc[i2];i3++){
			if(S[i2]==S[listc[i2][i3]]){
				obj+=conf_mat[i2][i3];
			}
			else{
				obj-=conf_mat[i2][i3];
			}
		}
	}
	//printf("%f\n",obj);
	return(obj);
}

float func_obSpos(int S[],float *conf_mat[],int *listc[],int numc[],int max_row,int pos){
	int i2,i3,k3;
	float obj;
	float l,l2,l3;
	obj=0.0;
	i2=pos;
	for(i3=0;i3<numc[i2];i3++){
		if(S[i2]==S[listc[i2][i3]]){
			obj+=conf_mat[i2][i3];
		}
		else{
			obj-=conf_mat[i2][i3];
		}
	}
	/* printf("%f\n",obj); */
	/* sleep(1); */
	return(obj);
}

float func_ob2(int X[],float *conf_mat[],int *listc[],int numc[],int max_row,int r_rank){
	int i2,i3,k3;
	float obj;
	float l,l2,l3;
	obj=0.0;
	for(i2=0;i2<max_row;i2++){
		for(i3=0;i3<numc[i2];i3++){
			obj+=conf_mat[i2][i3]*(float)X[i2]*(float)X[listc[i2][i3]];
		}
	}
	// printf("%f\n",obj);
	return(obj);
}

int rank(float *X[],float th,int max_row,int r_rank){
	int i,i2,i3;
	int l;
	int M=max_row;
	int N=r_rank;
	int LDU=5;
	int LDA=M;
	int LDVT=5;
	/* Locals */
	int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info, lwork;
	float wkopt;
	float* work;
	clock_t t1;
	int iwork[8*N];
	float s[N],u[LDU*M], vt[LDVT*N];
	srand((unsigned int)time(NULL));
	float a[LDA*N];
	int iter,iter2;
	for(iter=0;iter<N;iter++){
		for(iter2=0;iter2<LDA;iter2++){
			// a[iter*N+iter2]=((float)(rand())/RAND_MAX);
			a[iter*LDA+iter2]=X[iter2][iter];
			//printf("%f\t",a[iter*N+iter2]);
		}
	}

	lwork = -1;
	sgesdd_( "N", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt,
	&lwork, iwork, &info );
	lwork = (int)wkopt;
	work = (float*)malloc( lwork*sizeof(float) );
	/* Compute SVD */
	// printf("hi\n");
	sgesdd_( "N", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work,
	&lwork, iwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}
	// printf("Time for greedy %f\n",(float)(clock()-t1)/CLOCKS_PER_SEC);
	/* Print singular values */
	//  print_matrix( "Singular values", 1, n, s, 1 );
	l=r_rank;
	for(i=0;i<r_rank;i++){
		//  printf("%f\n",s[i]);
		if(s[i]>th){
			l=i+1;
			// printf("rank = %d %f\n",l,th);
		}
	}
	for(i=0;i<r_rank;i++){
		if(isnan(s[i])){
			//sleep(10);
		}
	}

	/* Print left singular vectors */
	//  print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu );
	/* Print right singular vectors */
	//   print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt );
	/* Free workspace */
	free( (void*)work );
	// printf("num sing %d\n",l);
	return(l);
}


/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, float* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}



float gaussian_random(){
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;
	if (iset == 0) { // We don't have an extra deviate handy, so
	do { v1=2.0*((float)(rand())/RAND_MAX)-1.0;// pick two uniform numbers in the square extending from -1 to +1 in each direction,
		v2=2.0*((float)(rand())/RAND_MAX)-1.0;
		rsq=v1*v1+v2*v2;// see if they are in the unit circle,
	}
	while (rsq >= 1.0 || rsq == 0.0);// and if they are not, try again.
	fac=sqrt(-2.0*log(rsq)/rsq);
	//  Now make the Box-Muller transformation to get two normal deviates. Return one and save the other for next time.
	gset=v1*fac;
	iset=1;// Set ag.
	return  v2*fac; }
	else {// We have an extra deviate handy,
		iset=0;// so unset the flag,
		return gset; //and return it.
	}
}
int cmpfunc (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}
int compare (const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}
