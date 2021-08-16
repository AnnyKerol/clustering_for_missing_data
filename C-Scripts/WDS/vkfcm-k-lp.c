/*#################################################################

VKFCM-K-LP - Whole data strategy - (WDS)
Version: 0.1
authors: Anny K G Rodrigues; Raydonal Ospina; Marcelo Ferreira

#################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//file pointers
FILE *e;
FILE *s;
FILE *par;
FILE *sig;
FILE *idx;
FILE *idx2;

struct partition //structure for the partition
{
    int n_clust; //cluster numbers
    int *hard;  //clusters for complete data
    int *hard2; 
    int *hard3; // clusters for missing data
    struct cluster *p; //pointer to cluster type structure
    double **pert;
};
struct cluster //Structure for the cluster
{
    int n_ind; //index number
    double *protot; //pointer that will store the address of the prototypes
    double *weight; // pointer for weights
};
struct dados //Structure for the data
{
    int *hard;
    int *hard2; 
    int *hard3; 
    int n_ind;
    int n_var; //number of variables
    double **x;
   // double **y; // matrix with no missing data
   // double **y_Miss; // matrix of missing
};

double **m2d_alloc(int n, int p, double **d);
double **aloca_mem_dados(int n, int p, double **d);
void leitura_dados(int n, int p, double **d); // function to read the data
void impressao_dados(int n, int p, double **d); //data printing function
double **aloca_mem_pert(int n, int p, double **d); // function for memory allocation of the membership matrix
void best_prototype(struct partition *p, struct dados *d, double m, double *sigma); // best prototype
void best_fuzzy_partition(struct partition *p, struct dados *d, double m, double *sigma); // best partition
void iniciar(struct partition *p, struct dados *d, double m, double *sigma); //function for initial partition
void inicial_best_fuzzy_partition(struct partition *p, struct dados *d, double m, int *sel, double *sigma);
void inicial_prototype(struct partition *p, struct dados *d, int *sel);
double criterio(struct partition *p, struct dados *d, double m, double *sigma); //calculates the J criterion
void particao_hard(struct partition *p, struct dados *d);
void inicial_best_weight(struct partition *p, struct dados *d);
void best_weight(struct partition *p, struct dados *d, double m, double *sigma);

double cr_index(struct partition *pp, struct dados *dd); 
void conf_mat(struct partition *p, struct dados *d, double **confmat); // confusion matrix
double misclassif_rate(struct partition *p, struct dados *d, double **confmat); 
double fmeasure(struct partition *p, struct dados *d, double **confmat);

double ** X_Completo(int n, int p, double **d); // complete Data 
double ** X_Missings(int n, int p, double **d); //data with missings
int verifica(struct dados *d, int k);
void distancia_parcial(struct partition *p,struct dados *d, double *sigma); // Partial Distance
int noMissing(double **d, int a,  int k);
void part_2(struct dados *d, struct partition *p, int **matriz, struct dados *f, int Num_Miss);

int main(int argc, char *argv[])
{
    int i,j,k,l,h,r,best_i,worst_i;
    int t, tmax, nrep,mconf;
    double m, eps, best_ant,best_atual, best_global,worst_global;
    double *cr, cr_mean, cr_sd, s_cr, ss_cr;
    double *mr, mr_mean, mr_sd, s_mr, ss_mr;
    double *fm, fm_mean, fm_sd, s_fm, ss_fm;
    double **confmat;
    double *sigma;
    struct dados *data, *data_Comp, *data_Miss;
    struct partition *particao;
    //char entrada[121];
    //char saida[121];
    //char param[121];
    char *entrada;
    char *saida;
    char *param;
    char *sigm;
    char *idxx;
    char *idxx2;

    //printf("Digite nome arquivo entrada\n");
    //scanf("%120s",entrada);
    //printf("Digite arquivo saida \n");
    //scanf("%120s",saida);
    //printf("Digite arquivo de parametros \n");
    //scanf("%120s",param);
    entrada = argv[1];
    saida = argv[2];
    param = argv[3];
    sigm = argv[4];
    idxx = argv[5];
      
 // input parameters
    e=fopen(entrada,"rt");
    s=fopen(saida,"wt");
    par=fopen(param,"rt");
    sig = fopen(sigm, "rt");
    idx = fopen(idxx, "wt");
    
    
    printf("Alocacao memoria dados 1\n");
    data = (struct dados *)malloc(sizeof(struct dados));
    data_Comp = (struct dados *)malloc(sizeof(struct dados));
    data_Miss = (struct dados *)malloc(sizeof(struct dados));
    //printf("num ind e num var:\n");
    //scanf("%d %d",&data->n_ind,&data->n_var);
    printf("Alocacao memoria particao 1\n");
    particao = (struct partition *)malloc(sizeof(struct partition));
    printf("Leitura dos parametros\n");
    fscanf(par, "%d %d %d %lf %d %lf %d", &data->n_ind, &data->n_var, &particao->n_clust, &m, &tmax, &eps, &nrep);
    printf("Alocacao memoria dados 2\n");
    data->x=aloca_mem_dados(data->n_ind,data->n_var,data->x); //allocating memory to the dataset
    data->hard3=(int *)malloc(data->n_ind*sizeof(int));
    printf("leitura dados\n");
    leitura_dados(data->n_ind,data->n_var,data->x);


  //counting the number of missings
    int Num_Miss=0;
	   for(j=0; j<data->n_ind; j++){
			  for(i=0; i<(data->n_var-1); i++){
				  if(data->x[j][i+1]==9999){
					   i=data->n_var+1;
					  Num_Miss=Num_Miss+1; 
					 
		  }
		 }  
	 }

 printf("Numero de missings %d\n", Num_Miss); 

	   int ii, b,**matriz;
         matriz = (int **)malloc(Num_Miss*sizeof(int *)); //allocating memory for an array of locations
         for(ii=0;ii<Num_Miss;ii++){
              matriz[ii]=(int*)malloc(2*sizeof(int));
           }   
           
           
// printf("Memoria alocada %d  %d\n", i,j); 
	     
	    int aux=0;   
	      for(j=0; j<data->n_ind; j++){
			  for(i=0; i<(data->n_var-1); i++){
				   if(data->x[j][i+1]==9999){
					        i=data->n_var+1; 
					   matriz[aux][0]=j;
					  // printf("Linhas dos missings %d\n", matriz[aux][0]);
					   aux++;
					   //j++;
					  // matriz[i][1]=i;//save the position of missings without the class column
					    
				  } 
			  }
	     } 
	     
	    
int c, N=Num_Miss;

data_Comp->n_ind = data->n_ind-Num_Miss;
data_Comp->n_var = data->n_var;

data_Comp->hard=(int *)malloc(data_Comp->n_ind*sizeof(int));

printf("data complete %d  %d\n", data_Comp->n_ind, data_Comp->n_var );
data_Comp->x =aloca_mem_dados(data_Comp->n_ind,data_Comp->n_var,data_Comp->x); //allocating memory to the complete database
data_Comp->x = X_Completo(data->n_ind,data->n_var,data->x); 



data_Miss->n_ind = Num_Miss;
data_Miss->n_var = data->n_var;
data_Miss->hard2=(int *)malloc(data_Miss->n_ind*sizeof(int));
data_Miss->x = aloca_mem_dados(data_Miss->n_ind,data_Miss->n_var,data_Miss->x);
data_Miss->x = X_Missings(data->n_ind,data->n_var, data->x);
 
             int linha, coluna;
             for(linha=0; linha<(Num_Miss); linha++){
			 for(coluna=0; coluna<(data->n_var); coluna++){
				printf("Meus dados completos: linha = %d  coluna= %d  %f\n ", linha, coluna, data_Miss->x[linha][coluna]);
				 
				 }
			
			}
 
    //printf("num grupos: \n");
    //scanf("%d",&particao->n_clust);
    printf("Alocacao memoria particao 2\n");
    particao->p=(struct cluster *)malloc(particao->n_clust*sizeof(struct cluster));
    printf("Alocacao memoria particao 3\n");
   
    particao->hard=(int *)malloc(data_Comp->n_ind*sizeof(int)); //banco completo
    particao->hard2=(int *)malloc(data_Miss->n_ind*sizeof(int)); //banco só com missings
    particao->hard3=(int *)malloc(data->n_ind*sizeof(int)); // banco original
    
    printf("Alocacao memoria particao 4\n");
    for(i=0;i<particao->n_clust;i++)
    {
        particao->p[i].protot=(double *)malloc((data->n_var-1)*sizeof(double));
    }
    for(i=0;i<particao->n_clust;i++)
    {
        particao->p[i].weight=(double *)malloc((data->n_var-1)*sizeof(double));
    }
    particao->pert=aloca_mem_pert(data_Comp->n_ind,particao->n_clust,particao->pert);
    //printf("Digite parametro m\n");
    //scanf("%lf",&m);
    //printf("Digite tmax\n");
    //scanf("%d",&tmax);
    //printf("Digite eps\n");
    //scanf("%lf",&eps);
    //printf("Digite numero repeticoes\n");
    //scanf("%d",&nrep);
    printf("Alocacao de memoria sigma\n");
    sigma = (double *)calloc(((data->n_var) - 1), sizeof(double));
    printf("Leitura sigma\n");
    for(j = 0; j < ((data->n_var) - 1); j++)
        fscanf(sig, "%lf", &sigma[j]);
    best_global = 9.9e+300;
    worst_global = 0.0;
    srand(time(NULL));

    /* Corrected Rand Index and Misclassification rate */
    cr = (double *)calloc(nrep, sizeof(double));
    mr = (double *)calloc(nrep, sizeof(double));
    fm = (double *)calloc(nrep, sizeof(double));

    for(i=0;i<nrep;i++)
    {
//        printf("Inicializacao\n");
        confmat = m2d_alloc(particao->n_clust, particao->n_clust, confmat);
        fprintf(s,"\nRepeticao %d",i+1);
        iniciar(particao,data_Comp,m,sigma);
        best_ant = criterio(particao,data_Comp,m,sigma);
        t=0;
        fprintf(s,"\nt=%d J=%f",t,best_ant);
       do
        {
            t=t+1;
            if(t> 1)
            {
                best_ant = best_atual;
            }
//            printf("Etapa 1: determinacao dos prototipos\n");
            best_prototype(particao,data_Comp,m,sigma);
//          printf("Etapa 2: determinacao da particao fuzzy\n");
            best_weight(particao,data_Comp,m,sigma);
            best_fuzzy_partition(particao,data_Comp,m,sigma);
            best_atual = criterio(particao,data_Comp,m,sigma);
            fprintf(s,"\nt=%d J=%f",t,best_atual);
          
        }
        while((fabs(best_ant-best_atual) > eps) && ( t < tmax));
        printf("\ni=%d J=%f\n",i,best_atual);
        
        particao_hard(particao,data_Comp);
        distancia_parcial(particao,data_Miss,sigma);
        part_2(data,particao,matriz, data_Comp, Num_Miss);
        
	 
	 
        if(best_global > best_atual)
        {
            best_i = i;
            best_global = best_atual;
        }
        if(worst_global < best_atual)
        {
            worst_i = i;
            worst_global = best_atual;
        } 
      
        
        cr[i] = cr_index(particao,data);
        s_cr += cr[i];
        ss_cr += cr[i] * cr[i];
        conf_mat(particao, data, confmat);
        mr[i] = misclassif_rate(particao, data, confmat);
        s_mr += mr[i];
        ss_mr += mr[i] * mr[i];
        fm[i] = fmeasure(particao, data, confmat);
        s_fm += fm[i];
        ss_fm += fm[i] * fm[i];
        fprintf(idx, "%d %.15lf %.15lf %.15lf %.15lf \n", i+1, best_atual, cr[i], mr[i], fm[i]);
        fprintf(s, "\n\nGrau de pertinencia\n");
        for(k = 0; k < data_Comp->n_ind; k++)
        {
            fprintf(s, "i=%d ", k+1);
            for(h = 0; h < particao->n_clust; h++)
            {
                fprintf(s, "%f ", particao->pert[k][h]);
            }
            fprintf(s, "%d %d \n", data_Comp->hard[k], particao->hard[k]);
        }
        fprintf(s,"\nPesos\n");
        for(k=0;k<particao->n_clust;k++)
        {
            fprintf(s,"Cluster %d: \n",k+1);
            for(j=0;j<(data_Comp->n_var-1);j++)
            {
				
                fprintf(s,"%f ",particao->p[k].weight[j]);
            }
            fprintf(s,"\n");
        }
        fprintf(s, "\nMatriz de Confusao\n");
        for(h = 0; h < particao->n_clust; h++)
        {
            for(l = 0; l < particao->n_clust; l++)
            {
                fprintf(s, "%d ", (int)confmat[h][l]);
            }
            fprintf(s, "\n");
        }
        fprintf(s, "\nClasses a priori\n");
        for(k = 0; k < particao->n_clust; k++)
        {
            fprintf(s, "Classe %d: \n", k+1);
            for(r = 0; r < data->n_ind; r++)
            {
                if(data->hard3[r] == (k+1))
                {
					
                    fprintf(s, "%d ", r+1);
                }
            }
            fprintf(s, "\n");
        }
        fprintf(s, "\nGrupos\n");
        for(k = 0; k < particao->n_clust; k++)
        {
            fprintf(s, "Grupo %d: \n", k+1);
            for(r = 0; r < data->n_ind; r++)
            {
                if(particao->hard3[r] == (k+1))
                {
                    fprintf(s, "%d ", r+1);
                }
            }
            fprintf(s, "\n");
        }
        free(confmat);
    }

    /* Compute the statistics */
    cr_mean = s_cr / nrep;
    cr_sd = pow(((ss_cr - nrep * cr_mean * cr_mean) / (nrep - 1)), 0.5);
    mr_mean = s_mr / nrep;
    mr_sd = pow(((ss_mr - nrep * mr_mean * mr_mean) / (nrep - 1)), 0.5);
    fm_mean = s_fm / nrep;
    fm_sd = pow(((ss_fm - nrep * fm_mean * fm_mean) / (nrep - 1)), 0.5);

    /* Print the results */
    fprintf(s, "\nFinal: best replication = %d, best J = %lf, worst replication = %d, worst J = %lf", best_i+1, best_global, worst_i+1, worst_global);
    fprintf(s, "\ncr_mean = %lf, cr_sd = %lf, best cr = %lf, worst cr = %lf, mr_mean = %lf, mr_sd = %lf, best mr = %lf, worst mr = %lf", cr_mean, cr_sd, cr[best_i], cr[worst_i], mr_mean, mr_sd, mr[best_i], mr[worst_i]);
    fprintf(s, "\nfm_mean = %lf, fm_sd = %lf, best_fm = %lf, worst_fm = %lf", fm_mean, fm_sd, fm[best_i], fm[worst_i]);

    free(matriz);
    free(cr);
    free(mr);
    free(fm);
    

    return 0;
}

double **m2d_alloc(int n, int p, double **d)
{
    int i;
    d = (double **)calloc(n, sizeof(double *));
    for(i = 0; i < n; i++)
        d[i] = (double *)calloc(p, sizeof(double));
    return d;
}
double **aloca_mem_dados(int n, int p, double **d)
{
    int i;
    d = (double **)malloc(n*sizeof(double *));
    for(i=0;i<n;i++)
    {
        d[i]=(double *)malloc(p*sizeof(double));
    }
    return d;
}
void leitura_dados(int n, int p, double **d)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            fscanf(e,"%lf",&d[i][j]);
        }
    }
}
void impressao_dados(int n, int p, double **d)
{
    int i, j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<p;j++)
        {
            fprintf(s,"%f ",d[i][j]);
        }
        printf("\n");
    }
}



double **aloca_mem_pert(int n, int k, double **d)
{
    int i;
    d = (double **)malloc(n*sizeof(double *));
    for(i=0;i<n;i++)
    {
        d[i]=(double *)malloc(k*sizeof(double));
    }
    return d;
}


int verifica(struct dados *d, int k){
	int j;
	int p=d->n_var;
	for(j=0; j<p; j++){
		   if(d->x[k][j]==9999){
			    
			    return k+1;
			   
			   } else{
				    
				    return k;
				   
				   }
		
		
	}

}
double ** X_Completo(int n, int p, double **d){ 
	  int i,j, k;
	   int Num_Miss=0;
	    for(j=0; j<n; j++){
			
			for(i=0; i<(p-1); i++){
				  if(d[j][i+1]==9999){
					  i=p+1;
					 Num_Miss=Num_Miss+1; 	 
		  }
		 }  
	 }
	 
	 
	 
	 // printf("numero de missing %d ", Num_Miss);
	 
	 int ii, b,**matriz;
         matriz = (int **)malloc(Num_Miss*sizeof(int *)); //allocating memory for an array of locations
         for(ii=0;ii<Num_Miss;ii++){
              matriz[ii]=(int*)malloc(2*sizeof(int));
           } 
           
          int aux=0;   
	      for(j=0; j<n; j++){
			  for(i=0; i<(p-1); i++){
				   if(d[j][i+1]==9999){
					        i=p+1; 
					   matriz[aux][0]=j;
					   aux++;
  
				  } 
			  }
	     } 
	     
	   
	 
	       double **y;
	       int m;       
	      
         y = (double **)malloc((n-Num_Miss)*sizeof(double *)); // alocando memória para a nova matriz 
             // printf(" passeii aquiiii  %d\n", n-Num_Miss);
         for(m=0;m<n-Num_Miss;m++){
              y[m]=(double*)malloc(p*sizeof(double));
             // printf(" Aloqueiiii \n");
              
           } 
	       
			/*for(i=0;i<Num_Miss;i++){
			   printf(" Verificaa  %d \n",matriz[i][0]);
	       }    
*/

			
	                            int linha, col, pp;          
			                    
			                    int aux2;
			                    
			                    aux=0;
					            for(linha=0; linha<n; linha++){
									
									aux2=0;
									for(pp=0; pp<Num_Miss; pp++) {
										   
									  if(linha!=matriz[pp][0]){
									
										for(col=0; col<p; col++){     
												 y[aux][col]=d[linha][col];
										}
										aux2++;
									  }
 								 	
								   }
									
											
									if(aux2==Num_Miss){
										aux++;
									}
 										
								 
							 }
	   
	       
	 
	// int N=n-Num_Miss;
 //printf("tamanho de y ==== %d\n ", N);
	//int lin, col2;
	//for(lin=0; lin<N; lin++){
	  // for(col2=0; col2<p; col2++){
		
	//}
	//}	 
return y;
free(matriz);
free(y);
}

double ** X_Missings(int n, int p, double **d){
	 
	  int i,j, k;
	   int Num_Miss=0;
	    for(j=0; j<n; j++){
			
			for(i=0; i<(p-1); i++){
				  if(d[j][i+1]==9999){
					  i=p+1;
					 Num_Miss=Num_Miss+1; //contando o número de missing		 
		  }
		 }  
	 }
	 
	 
	 
	 // printf("numero de missing %d ", Num_Miss);
	 
	   int ii, b,**matriz;
         matriz = (int **)malloc(Num_Miss*sizeof(int *)); 
         for(ii=0;ii<Num_Miss;ii++){
              matriz[ii]=(int*)malloc(2*sizeof(int));
           } 
           
          int aux=0;   
	      for(j=0; j<n; j++){
			  for(i=0; i<(p-1); i++){
				   if(d[j][i+1]==9999){
					        i=p+1; //gambiarra
					   matriz[aux][0]=j;
					   aux++;
					   
				  } 
			  }
	     } 
	     
	   
	 
	       double **y;
	       int m;       
	      
         y = (double **)malloc((Num_Miss)*sizeof(double *)); // alocando memória para a nova matriz 
             // printf(" passeii aquiiii  %d\n", n-Num_Miss);
         for(m=0;m<Num_Miss;m++){
              y[m]=(double*)malloc(p*sizeof(double));
             // printf(" Aloqueiiii \n");
              
           } 
	       
/*			for(i=0;i<Num_Miss;i++){
			   printf("%d \n",matriz[i][0]);
	       }
*/	       


			
	                            int linha, col, pp;          
			                    
			                    int aux2;
			                    
			                    aux=0;
					            for(linha=0; linha<n; linha++){
									
									aux2=0;
									for(pp=0; pp<Num_Miss; pp++) {
										   
									  if(linha==matriz[pp][0]){
										aux2++;
									  }
 								 	
								   }
									if(aux2==1){
										for(col=0; col<p; col++){
										y[aux][col]=d[linha][col];
										}
										aux++;
									}
									//printf("%d %d %d\n",aux,linha,aux2);
								 
							 }
	   
	       
	 
	// int N=Num_Miss;
//	 printf("tamanho de y ==== %d\n ", N);
	//int lin, col2;
	//for(lin=0; lin<N; lin++){
	 //  for(col2=0; col2<p; col2++){
		

		
	//}
//	}
	
	 
return y;
free(matriz);
free(y);
}


void iniciar(struct partition *p, struct dados *d, double m, double *sigma)
{
    int i,j,k,l,flag,*sel;
    sel=(int *)malloc(p->n_clust*sizeof(int));
    do
    {
        for(i=0; i<p->n_clust;i++)
        {
            sel[i] = rand() % (d->n_ind-1);
            fprintf(s,"sel[%d]=%d ",i,sel[i]);
        }
        fprintf(s,"\n");
        flag=0;
        for(i=0; i<(p->n_clust-1);i++)
        {
            k=sel[i];
            for(j=i+1; j<p->n_clust;j++)
            {
                l=sel[j];
                if(k==l)
                {
 //               printf("sel[%d]=%d sel[%d]=%d ",i, sel[i],j,sel[j]);
                    flag=1;
                    break;
                }
            }
        }
 //       printf("flag=%d ",flag);
    }
    while(flag==1);
    inicial_prototype(p,d,sel);
    inicial_best_weight(p,d);
    inicial_best_fuzzy_partition(p,d,m,sel,sigma);
    free(sel);
}
void inicial_best_weight(struct partition *p, struct dados *d)
{
    int j, k;
    for(k=0;k<p->n_clust;k++)
    {
        for(j=0;j<(d->n_var-1);j++)
        {
            p->p[k].weight[j] = 1.0;
        }
    }
}
void inicial_prototype(struct partition *p, struct dados *d, int *sel)
{
    int i,j;
    for(i=0; i<p->n_clust;i++)
    {
        for(j=0; j<(d->n_var-1);j++)
        {
            p->p[i].protot[j] = d->x[sel[i]][j+1];
        }
    }
/*    for(i=0; i<p->n_clust;i++)
    {
        for(j=0; j<(d->n_var-1);j++)
        {
            printf("%f ",p->p[i].protot[j]);
        }
        printf("\n");
    }
    system("pause");
    exit(1);*/
}
void best_prototype(struct partition *p, struct dados *d, double m, double *sigma)
{
    int i,j,k;
    double num, den;
    for(i=0; i<p->n_clust;i++)
    {
        for(j=0; j<(d->n_var-1);j++)
        {
            num = 0.0;
            den = 0.0;
            for(k=0;k<d->n_ind;k++)
            {
                num = num + pow(p->pert[k][i],m) * exp(-0.5* (d->x[k][j+1] - p->p[i].protot[j]) * (d->x[k][j+1] - p->p[i].protot[j])/sigma[j]) * d->x[k][j+1];
                den = den + pow(p->pert[k][i],m) * exp(-0.5* (d->x[k][j+1] - p->p[i].protot[j]) * (d->x[k][j+1] - p->p[i].protot[j])/sigma[j]);
            }
            p->p[i].protot[j] = num / den;
        }
    }
}
void best_fuzzy_partition(struct partition *p, struct dados *d, double m, double *sigma)
{
    int i,j,k,h,l;
    double num, den,frac, s, esp;
    esp = 1.0 / (m-1);
    for(i=0;i<d->n_ind;i++)
    {
        for(k=0;k<p->n_clust;k++)
        {
            num=0.0;
            for(j=0;j<(d->n_var-1);j++)
            {
                num=num + p->p[k].weight[j] * (1.0 - exp(-0.5* (d->x[i][j+1] - p->p[k].protot[j]) * (d->x[i][j+1] - p->p[k].protot[j])/sigma[j]));
            }
            s = 0.0;
            for(h=0;h<p->n_clust;h++)
            {
                den = 0.0;
                for(l=0;l<(d->n_var-1);l++)
                {
                    den = den + p->p[h].weight[l] * (1.0 - exp(-0.5* (d->x[i][l+1] - p->p[h].protot[l]) * (d->x[i][l+1] - p->p[h].protot[l])/sigma[l]));
                }
                frac = num / den;
                s = s + pow(frac, esp);
            }
            p->pert[i][k] = 1.0 / s;
        }
    }
}
void best_weight(struct partition *p, struct dados *d, double m, double *sigma)
{
    int i,j,k,h,l;
    double num, den,partial,eps;
    eps=1.0/(d->n_var-1);
    for(k=0;k<p->n_clust;k++)
    {
        for(j=0;j<(d->n_var-1);j++)
        {
            den = 0.0;
            for(i=0;i<d->n_ind;i++)
            {
                den = den + pow(p->pert[i][k],m)*(1.0 - exp(-0.5* (d->x[i][j+1] - p->p[k].protot[j]) * (d->x[i][j+1] - p->p[k].protot[j])/sigma[j]));
            }
            num=1.0;
            for(h=0;h<(d->n_var-1);h++)
            {
                partial=0.0;
                for(l=0;l<d->n_ind;l++)
                {
                    partial = partial + pow(p->pert[l][k],m)*(1.0 - exp(-0.5* (d->x[l][h+1] - p->p[k].protot[h]) * (d->x[l][h+1] - p->p[k].protot[h])/sigma[h]));
                }
                partial=pow(partial,eps);
                num = num *partial;
            }
            p->p[k].weight[j] = num/den;
/*            if(den > 0.0)
            {
                p->p[k].weight[j] = num/den;
            }
            else
            {
                p->p[k].weight[j]=1.0;
            }*/
//            printf("%d %d %f %f %f",k,j,num, den, p->p[k].weight[j]);
//            system("pause");
        }
    }
}

void inicial_best_fuzzy_partition(struct partition *p, struct dados *d, double m, int *sel, double *sigma)
{
    int i,j,k,h,l,count;
    double num, den,frac, s, esp;
    esp = 1.0 / (m-1);
    for(k=0;k<p->n_clust;k++)
    {
        for(i=0;i<d->n_ind;i++)
        {
            if(i==sel[k])
            {
                p->pert[i][k] = 1.0;
            }
            else
            {
                p->pert[i][k] = 0.0;
            }
        }
    }
   for(i=0;i<d->n_ind;i++)
    {
        for(k=0;k<p->n_clust;k++)
        {
            count = 0;
            if(i==sel[k])
            {
                ;
            }
            else
            {
                num=0.0;
                for(j=0;j<(d->n_var-1);j++)
                {
                    num=num + (1.0 - exp(-0.5* (d->x[i][j+1] - p->p[k].protot[j]) * (d->x[i][j+1] - p->p[k].protot[j])/sigma[j]));
                }
                s = 0.0;
                for(h=0;h<p->n_clust;h++)
                {
                    den = 0.0;
                    for(l=0;l<(d->n_var-1);l++)
                    {
                        den = den + (1.0 - exp(-0.5* (d->x[i][l+1] - p->p[h].protot[l]) * (d->x[i][l+1] - p->p[h].protot[l])/sigma[l]));
                    }
                    if(num > 0.0)
                    {
                        frac = num / den;
                        s = s + pow(frac, esp);
                    }
                    else
                    {
                        count = count + 1;
                        s=1.0;
                    }
                }
                p->pert[i][k] = 1.0 / s;
                if(count > 0)
                {
                    p->pert[i][k] =  p->pert[i][k] /count;
                }
            }
        }
    }
}
double criterio(struct partition *p, struct dados *d, double m, double *sigma)
{
    int i,j,k;
    double w=0.0;
    double s;
    for(k=0;k<p->n_clust;k++)
    {
        for(i=0;i<d->n_ind;i++)
        {
            s=0.0;
            for(j=0;j<(d->n_var-1);j++)
            {
                s = s + p->p[k].weight[j] * 2.0 * (1.0 - exp(-0.5* (d->x[i][j+1] - p->p[k].protot[j]) * (d->x[i][j+1] - p->p[k].protot[j])/sigma[j]));
            }
            w = w + pow(p->pert[i][k],m)*s;
       }
    }
    return w;
}


void distancia_parcial(struct partition *p, struct dados *d, double *sigma)
{

    int i,j,k;
    int  mindist;
    double s;
    int I, cont=0;
    double aux=9999;
    int grupo;
    	
        for(i=0;i<d->n_ind;i++)
        {
			aux=9999;
            for(k=0;k<p->n_clust;k++) 
			{
			
				s=0.0; 
           
            for(j=0;j<(d->n_var-1);j++)
            {
				
				
				if(d->x[i][j+1]!= 9999){
			
					    I=1;
                        cont=noMissing(d->x,d->n_var,i);
                        s = s + ((d->n_var-1)/cont)*(p->p[k].weight[j] * 2.0 * (1.0 - exp(-0.5* (d->x[i][j+1] - p->p[k].protot[j]) * (d->x[i][j+1] - p->p[k].protot[j])/sigma[j])))*I;
			           // printf(" Cluster === %d centros ===  %f  pesos=== %f  matriz de dados==== %f\n ", k, p->p[k].protot[j], p->p[k].weight[j], d->x[i][j+1]);
			              
			          }	           
                     
            }
            
             
         //printf( "row === %d  Cluster === %d  distance in each cluster  ===== %f\n ",i,k,s);
           //  printf(" \n ");

       
       
      
          if(s<aux){
			 grupo=k;
			 aux=s;
			 
			 }
	
			
       } 
      	
          p->hard2[i]=grupo+1; 
       //  printf("i === %d menor grupo  ===== %d\n ", i, p->hard2[i]); 

    }
    
    //return s;
}







void particao_hard(struct partition *p, struct dados *d)
{
    int i,k;
    double maxpert;
    int max;
    for(i=0;i<d->n_ind;i++)
    {
        maxpert=p->pert[i][0];
        max=0;
        for(k=0;k<p->n_clust;k++)
        {
            if(maxpert < p->pert[i][k])
            {
                maxpert = p->pert[i][k];
                max = k;
            }
//            printf("%d %d %d %f %f \n",i,k,max, p->pert[i][k], maxpert);
//            system("pause");
        }
        p->hard[i] = max+1;
       // printf("meu vetor segmentado %d %d\n ", i, p->hard[i]);
        d->hard[i] = d->x[i][0];
    }
/*    for(i=0;i<d->n_ind;i++)
    {
        for(k=0;k<p->n_clust;k++)
        {
            printf("%f ",p->pert[i][k]);
        }
        printf("%d %d \n", d->hard[i], p->hard[i]);
        system("pause");
    }*/
}

void part_2(struct dados *d, struct partition *p, int **matriz, struct dados *f, int Num_Miss){
	

	int i,aux,j, aux2, aux3;


//	printf("%d\n",Num_Miss);

	aux=0;
	aux2=0,aux3=0;
	
	for(i=0; i<d->n_ind; i++){
		if((aux < Num_Miss) && (i==matriz[aux][0])){
					
					
			p->hard3[i]=p->hard2[aux];
					
			//printf(" linha === %d aux === %d hard2 (just Missings) ==== %d\n ", i, aux, p->hard3[i]);
			aux++;

		}else{
				
			
			p->hard3[i]=p->hard[aux3];
			//printf(" linha == %d aux3 == %d hard (complete) == %d\n ", i, aux3, p->hard3[i]);
			aux3++;
				
			}

		//printf(" i== %d hard3 (Todos) === %d\n", i, p->hard3[i]);
		d->hard3[i]=d->x[i][0];
		}
}



double cr_index(struct partition *pp, struct dados *dd)
{
    int i, j;
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0, p, cr;
    for(i = 0; i < (dd->n_ind); i++)
    {
        for(j = (i+1); j < dd->n_ind; j++)
        {
			
            if((pp->hard3[i] == pp->hard3[j]) && (dd->hard3[i] == dd->hard3[j]))
                a += 1.0;
            if((pp->hard3[i] != pp->hard3[j]) && (dd->hard3[i] == dd->hard3[j]))
                b += 1.0;
            if((pp->hard3[i] == pp->hard3[j]) && (dd->hard3[i] != dd->hard3[j]))
                c += 1.0;
            if((pp->hard3[i] != pp->hard3[j]) && (dd->hard3[i] != dd->hard3[j]))
                d += 1.0;

        }
    }
    p = a + b + c + d;
    cr = ((a + d) - ((a + b) * (a + c) + (d + b) * (d + c)) * (1 / p)) / (p - ((a + b) * (a + c) + (d + b) * (d + c)) * (1 / p));
    //float h,t;
    //h=(a + d);//((a + b) * (a + c) + (d + b) * (d + c))* (1 / p);//((a + d) - ((a + b) * (a + c) + (d + b) * (d + c)) * (1 / p));
    //t=(p - ((a + b) * (a + c) + (d + b) * (d + c)) * (1 / p));
    //printf("%f\n",h);
    //printf("%f\n",t);
    return cr;
}

void conf_mat(struct partition *p, struct dados *d, double **confmat)
{
    int i, h, k;
    for(i = 0; i < p->n_clust; i++)
    {
        for(h = 0; h < p->n_clust; h++)
        {
            for(k = 0; k < d->n_ind; k++)
            {
                if(p->hard3[k] == (i + 1) && d->hard3[k] == (h + 1))
                    confmat[i][h] += 1;
                    //printf("%d %d %d %f\n",i,p->hard3[k],d->hard3[k] == (h + 1),confmat[i][h]);
            }
        }
    }
}

double misclassif_rate(struct partition *p, struct dados *d, double **confmat)
{
    int i, h;
    double max_mr, sum_mr = 0.0, tmp, mr;
    for(i = 0; i < p->n_clust; i++)
    {
        max_mr = confmat[i][0];
        for(h = 0; h < p->n_clust; h++)
        {
            if(max_mr < confmat[i][h])
                max_mr = confmat[i][h];
        }
        sum_mr += max_mr;
    }
    tmp = sum_mr / (d->n_ind);
    mr = 1.0 - tmp;
    return mr;
}
double fmeasure(struct partition *p, struct dados *d, double **confmat)
{
    int i, j;
    double *total_linha, *total_col, *maxfmeasure;
    double **precision, **recall, **fmeasure;
    double F = 0.0;

    total_linha = (double *)calloc(p->n_clust, sizeof(double));
    total_col = (double *)calloc(p->n_clust, sizeof(double));
    maxfmeasure = (double *)calloc(p->n_clust, sizeof(double));
    precision = m2d_alloc(p->n_clust, p->n_clust, precision);
    recall = m2d_alloc(p->n_clust, p->n_clust, recall);
    fmeasure = m2d_alloc(p->n_clust, p->n_clust, fmeasure);

    for(i = 0; i < p->n_clust; i++)
    {
        for(j = 0; j < p->n_clust; j++)
        {
            total_linha[i] += confmat[i][j];
            total_col[j] += confmat[i][j];
        }
    }
    for(i = 0; i < p->n_clust; i++)
    {
        for(j = 0; j < p->n_clust; j++)
        {
            precision[i][j] = confmat[i][j] / total_linha[i];
            recall[i][j] = confmat[i][j] / total_col[j];
        }
    }
    for(i = 0; i < p->n_clust; i++)
    {
        for(j = 0; j < p->n_clust; j++)
        {
            if((precision[i][j] + recall[i][j]) == 0)
                fmeasure[i][j] = 0.0;
            else
                fmeasure[i][j] = 2 * ((precision[i][j] * recall[i][j]) / (precision[i][j] + recall[i][j]));
        }
    }
    for(i = 0; i < p->n_clust; i++)
    {
        maxfmeasure[i] = fmeasure[0][i];
        for(j = 0; j < p->n_clust; j++)
        {
            if(maxfmeasure[i] < fmeasure[j][i])
                maxfmeasure[i] = fmeasure[j][i];
        }
    }
    for(i = 0; i < p->n_clust; i++)
    {
        F += (total_col[i] / d->n_ind) * maxfmeasure[i];
    }

    return F;

    free(total_linha);
    free(total_col);
    free(precision);
    free(recall);
    free(fmeasure);
    free(maxfmeasure);
}

int noMissing(double **d, int a,  int k){
	int i, cont=0;
	for(i=0; i<a-1; i++){
		   if(d[k][i+1]!=9999){
			   cont=cont+1;
			   }
		
		}
	return cont;
	}

