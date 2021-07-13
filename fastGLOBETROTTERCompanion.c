#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include "zlib.h"
 
// SCALE "chrom_ind_res" BY TEMP-WEIGHTS BEFORE STORING, TO REDUCE NUMBER OF DONORS (SLIGHTLY AT LEAST)

//typedef voidp gzFile;       /* opaque gzip file descriptor */

//ZEXTERN gzFile ZEXPORT gzopen OF((const char *path, const char *mode));
// This software is licenced under the "Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License" (http://creativecommons.org/licenses/by/4.0/)
void mutiply_temp_weight(double * tempweights,double *  weights_mat, int * num_surrogates, int * num_donors)
{    
    for (int k=0; k < *num_donors; k++)
		 for (int l=0; l < *num_donors; l++)
			  for (int m=0; m < *num_surrogates; m++)
				   for (int n=0; n < *num_surrogates; n++)
					    tempweights[(m*(*num_surrogates)+n)*(*num_donors)*(*num_donors)+k*(*num_donors)+l] = weights_mat[k*(*num_surrogates)+m]*weights_mat[l*(*num_surrogates)+n];

}



int reading(st, format, res)
    char **st, *format;
    void *res;
{
    int i;
    char *rs;
    rs = *st;
    for(i = 0; isspace(rs[i]); i++) ; 
    if (!rs[i]) return 0; 
    for(; !isspace(rs[i]); i++) ;
    if (rs[i]) rs[i++] = 0;  
    if (!sscanf(*st, format, res)) return 0; 
    *st += i;
    return 1;
}





void tabulate_chunks(double * geneticdistend_all,int nsamplestoconsider,int *ploidy,int ind_num, double *binwidth, int *num_bins, int *grid_start, int *num_donors, int num_chunksALL, int * chunk_donor_popALL, double * chunk_posALL, double * grid_sizeALL, double * gridweight_GENDIST_ALL, int *num_surrogates,double *** total_counts_mat_all)
{   
	
    int count_ans, times, chunk1,chunk2;
    int i,j,k,c;
    c = 0;
    int bin_spot_i_j;
    double dist_btwn_i_j,min_dist_btwn_i_j;
    
       
	 double chr_length = 0.0;
	
	   for (i=0; i <num_chunksALL; i++){
        if (geneticdistend_all[i] >= chr_length)
		{
			chr_length = geneticdistend_all[i];
		}
			
       }
   	  //double chr_length = geneticdistend_all[num_chunksALL-1];
	  int bin_size_const = 1; 
	  int num_bin = (int)chr_length/bin_size_const+1;
  
       int *bin = malloc(num_chunksALL*sizeof(int));
       int *numchunckinbin = malloc((num_bin)*sizeof(int));
       for (i=0; i <num_bin; i++)
       numchunckinbin[i] =0;
	   
   	   for (i=0; i <num_chunksALL; i++){
	   bin[i] = (int)chunk_posALL[i];
       }
  	   for (i=0; i <num_chunksALL; i++)
  	   { numchunckinbin[bin[i]] =  numchunckinbin[bin[i]]+1;
	   }
       int *start_bin = malloc((num_bin)*sizeof(int));
	   start_bin[0] = 0;
	   for (i=1; i <num_bin; i++)
	   {
				start_bin[i] = start_bin[i-1] + numchunckinbin[i-1];
		} 
	
   	double *tabulate_chunk_posALL = malloc((num_chunksALL)*sizeof(double));
   	double *tabulate_grid_sizeALL = malloc((num_chunksALL)*sizeof(double));
   	double *tabulate_gridweight_GENDIST_ALL = malloc((num_chunksALL)*sizeof(double));
   	int *tabulate_chunk_donor_popALL = malloc((num_chunksALL)*sizeof(int));
    int cc = 0;
		
	for ( j=0; j <num_bin; j++)
   	{
   		for ( k=0; k <num_chunksALL; k++){
   		if (bin[k] == j){
   		tabulate_chunk_posALL[cc] = chunk_posALL[k];
   		tabulate_grid_sizeALL[cc] = grid_sizeALL[k];
   		tabulate_gridweight_GENDIST_ALL[cc] = gridweight_GENDIST_ALL[k];
   		tabulate_chunk_donor_popALL[cc] = chunk_donor_popALL[k];
		if (tabulate_chunk_donor_popALL[cc] > *num_donors) 
			{
				printf("%i %i\n",tabulate_chunk_donor_popALL[cc],*num_donors);
				printf("The painting samples and ID file are not matched, please check number of individuals in the ID file. Exiting....\n");exit(1);
	
	         }
   		   	    cc = cc+1;
   		}
   	   }
   	}
    
   
       
	
		
	    int window_size = ceil(*binwidth*(*num_bins))+3;
      	count_ans = 0;
      	int t_1, t_2, j_lim;
		int Y_i;
      	for ( j=0; j <num_bin-1; j++)
      	{
      	j_lim =j+window_size;
      	if (j+window_size > num_bin){j_lim =num_bin;}
          for ( k=j+1; k <j_lim; k++)
      	{
   			 t_1 = numchunckinbin[j];
			 t_2 = numchunckinbin[k];
 			 Y_i = (int)t_1*t_2*exp(-0.05*(k-j))/10;
				if (t_1>0 && t_2>0 &&  Y_i == 0) {Y_i = 1;}
      		count_ans = count_ans+Y_i;
      		
      	    times = 0;		
      		while(times < Y_i) 
		 		{
         		  chunk1 = start_bin[j] + rand() %  (t_1); 
				  chunk2 = start_bin[k] + rand() %  (t_2); 
				  if (tabulate_chunk_donor_popALL[chunk2] == 0 ||tabulate_chunk_donor_popALL[chunk1] == 0 )
				  {
					  times = times+1;
					  continue;
				  }
	  		 dist_btwn_i_j=pow(pow((tabulate_chunk_posALL[chunk1]-tabulate_chunk_posALL[chunk2]),2.0),0.5);
      		 if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) <= ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)));
      		 if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) > ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)))+1;
	 
      		 bin_spot_i_j = bin_spot_i_j - (*grid_start);
	 
      		 min_dist_btwn_i_j=tabulate_grid_sizeALL[chunk1]/2.0+tabulate_grid_sizeALL[chunk2]/2.0;
      	     if (bin_spot_i_j >= 0 && bin_spot_i_j<(*num_bins) && dist_btwn_i_j >= min_dist_btwn_i_j) 
      			 {
      	
      			 if ((tabulate_chunk_donor_popALL[chunk1]-1) > (tabulate_chunk_donor_popALL[chunk2]-1))
      			 {
      			  total_counts_mat_all[ind_num][((tabulate_chunk_donor_popALL[chunk2]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk1]-1))][bin_spot_i_j]=total_counts_mat_all[ind_num][((tabulate_chunk_donor_popALL[chunk2]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk1]-1))][bin_spot_i_j]+tabulate_gridweight_GENDIST_ALL[chunk2]*tabulate_gridweight_GENDIST_ALL[chunk1];
      			
				 }	
      		     else 
      			 {
			 
      			  total_counts_mat_all[ind_num][((tabulate_chunk_donor_popALL[chunk1]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk2]-1))][bin_spot_i_j]=total_counts_mat_all[ind_num][((tabulate_chunk_donor_popALL[chunk1]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk2]-1))][bin_spot_i_j]+tabulate_gridweight_GENDIST_ALL[chunk1]*tabulate_gridweight_GENDIST_ALL[chunk2];      		     
				 }
				 }
	 
      	 	times = times+1;
      	    }
          }
          }

	  
	 
  free(numchunckinbin);
  free(bin);
  free(start_bin);
  free(tabulate_chunk_donor_popALL);
  free(tabulate_grid_sizeALL);
  free(tabulate_gridweight_GENDIST_ALL);
  free(tabulate_chunk_posALL);		  


}

void data_read(int *ploidy, int * ind_id_vec, char **filenameSAMPread, char **filenameRECOMread, double *binwidth, int *num_bins, int *grid_start, int *num_donors, int * donor_label_vec, int *num_inds_rec, char **recipient_label_vec, double *weight_min, double *gridweightMAX_cutoff, int * num_surrogates, double * temp_weights, double * means)
{
	
   
  	
  int h,i,j,p,s,n,k;
  int num_chrom, num_chrom2, nsamples, nsamplestoconsider, nsamplestoconsiderNULL, count, nhaps, ninds, nsites, numchunks, num_chunksALL, to_skip;
  int chunkcount, chunk_count_all, line_check, rec_count, hap_count;
  double totalgeneticdist, gridweight_val;
  FILE *fd, *fd3, *fd2;
  //gzFile *fd2;
  char *step;
  char * line = malloc(10000000 * sizeof(char));
  char waste[400];
  char waste2[400];
  char * filename = malloc(1000 * sizeof(char));
  char * filenameSAMP = malloc(1000 * sizeof(char));

  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  num_chrom=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom = num_chrom + 1;
    }
  fclose(fd);
  fd = fopen(*filenameRECOMread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
  num_chrom2=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom2 = num_chrom2 + 1;
    }
  fclose(fd);
  if (num_chrom != num_chrom2){printf("The number of files in %s and %s do not match. Exiting....\n",*filenameSAMPread,*filenameRECOMread);exit(1);}

	       // GET NUMBER OF RECIPIENT INDS AND SAMPLES:
  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  strcpy(filenameSAMP,waste);
  fclose(fd);
  fd2 = gzopen(filenameSAMP,"r");
  if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
  line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  while(strcmp(waste,"nsamples")!=0) 
    reading(&step,"%s",waste);
  reading(&step,"%s",waste);
  reading(&step,"%d",&nsamples);
  nsamplestoconsider=nsamples;
  if (nsamplestoconsider<=0) { printf("No painting samples in %s! Exiting....\n",filenameSAMP);exit(1);}
  //nsamplestoconsiderNULL=nsamples;
  nsamplestoconsiderNULL=1;
  count = 0;
  rec_count=0;
  hap_count=0;
  int * rec_ind_first_hap_index_vec = malloc((*num_inds_rec)*sizeof(int));
  while(!gzeof(fd2))
    {
      if(gzgets(fd2,line,10000000)!=NULL);
      step=line;
      reading(&step,"%s",waste);
      if(strcmp(waste,"HAP")==0) 
	{
	  reading(&step,"%s",waste2);
	  reading(&step,"%s",waste2);
	  for (i=0; i < *num_inds_rec; i++)
	    {
	      if (strcmp(waste2,recipient_label_vec[i])==0)
		{
		  if (rec_count>=((*ploidy)*(*num_inds_rec)))
		    {
		      printf("There are more than %d recipient haplotypes in %s! Exiting....\n",(*num_inds_rec)*(*ploidy),filenameSAMP);
		      exit(1);
		    }
		  rec_ind_first_hap_index_vec[(rec_count/(*ploidy))]=hap_count/(*ploidy);
		  rec_count=rec_count+1;
		  break;
		}
	    }
	  hap_count=hap_count+1;
	}
      count = count + 1;
    }
  gzclose(fd2);
  int ** ind_id_vec_new = malloc((*num_inds_rec)*sizeof(int*));
  for (i=0; i < (*num_inds_rec); i++)
    ind_id_vec_new[i]=malloc(num_chrom*sizeof(int));
  for (i=0; i < (*num_inds_rec); i++)
    {
      for (h=0; h < num_chrom; h++)
	ind_id_vec_new[i][h]=rec_ind_first_hap_index_vec[ind_id_vec[(i*num_chrom+h)]];
    }
  nhaps=count/(nsamples+1);
  if (nhaps != (((double)count)/(nsamples+1)) && nhaps != (((double)count-1)/(nsamples+1)))
    {
  
      printf("something wrong with file %s. Exiting....\n",filenameSAMP);
      exit(1);
    }
  ninds=nhaps/(*ploidy);
  if (ninds<=0){ printf("No individuals in %s! Exiting....\n",filenameSAMP);exit(1);}
  if (rec_count!=((*num_inds_rec)*(*ploidy)))
    {
      printf("Only %d recipient haplotypes in %s, while expecting %d! Exiting....\n",rec_count,filenameSAMP,(*num_inds_rec)*(*ploidy));
      exit(1);
    }
  
	
	  double *** total_counts_mat_all = malloc((*num_inds_rec)*sizeof(double **));
	  for (n=0; n < (*num_inds_rec); n++)
	    {
	      total_counts_mat_all[n]=malloc((*num_donors) * (*num_donors) * sizeof(double *));
	      for (k=0; k < ((*num_donors) * (*num_donors)); k++)
		{
		  total_counts_mat_all[n][k]=malloc((*num_bins)*sizeof(double));
		}
	    }
	  for (n=0; n < (*num_inds_rec); n++)
	    {
	      for (k=0; k < ((*num_donors) * (*num_donors)); k++)
		{
		  for (i=0; i < (*num_bins); i++)
		    total_counts_mat_all[n][k][i]=0.0;
		}
	    }
	
	
	  double *** results_all = malloc((*num_inds_rec)*sizeof(double **));
	  for (n=0; n < (*num_inds_rec); n++)
	    {
	      results_all[n]=malloc((*num_surrogates) * (*num_surrogates) * sizeof(double *));
	      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
		{
		  results_all[n][k]=malloc((*num_bins)*sizeof(double));
		}
	    }
	  for (n=0; n < (*num_inds_rec); n++)
	    {
	      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
		{
		  for (i=0; i < (*num_bins); i++)
		    results_all[n][k][i]=0.0;
		}
	    }
	
	//start call tabulate here
  for (h=0; h < num_chrom; h++)
    {
                       /* GET NUMBER OF SNPs: */
      fd = fopen(*filenameRECOMread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if(fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameRECOMread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filename,waste);
      fclose(fd);
      
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      nsites = -1;
      while(!feof(fd3))
	{
	  if(fgets(line,2047,fd3)!=NULL);
	  nsites = nsites + 1;
	}
      fclose(fd3);
      if (nsites<=0) {printf("Something wrong with %s -- no SNPs! Exiting....\n",filename); exit(1);}
	  
                       /* GET RECOMBINATION RATES AND SNP POSITIONS: */
      double * posvec = malloc(nsites * sizeof(double));
      double * recomrateperbp = malloc(nsites * sizeof(double));
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      for (j=0; j < nsites; j++)
	{
	  line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s at row %d. Exiting....\n",filename,j+1); exit(1);}
	  step=line;
	  reading(&step,"%lf",&posvec[j]);    // basepair position
	  reading(&step,"%lf",&recomrateperbp[j]);
	}
      fclose(fd3);
      for (i=0; i < nsites; i++)
	{
	  if (i > 0)
	    {
	      if (posvec[i]<posvec[(i-1)])
		{
		  printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filename,posvec[(i-1)],posvec[i]);
		  exit(1);
		}
	    }
	  if (recomrateperbp[i]<0)
	    {
	      printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",posvec[i]);
	      exit(1);
	    }
	}
      
      fd = fopen(*filenameSAMPread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filenameSAMP,waste);
      fclose(fd);

               /* IF NOT DOING NULL IND, CONSIDER COMPARISONS BETWEEN ALL CHUNKS: */
                       /* FIND NUMBER OF CHUNKS ACROSS ALL SAMPLES FOR IND: */
      fd2 = gzopen(filenameSAMP,"r");
      if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
      line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
      n=0;
      to_skip=(*ploidy)*(nsamples+1)*ind_id_vec_new[n][h];
      while (n < (*num_inds_rec))
	{
	  for (i=0; i < to_skip; i++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
	      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
	    }
	  
	  int ** ind_haps_mat=malloc((*ploidy*nsamplestoconsider)*sizeof(int *));
	  for (i=0; i < (*ploidy*nsamplestoconsider); i++)
	    ind_haps_mat[i]=malloc(nsites*sizeof(int));
	  int * numchunks_vec=malloc((*ploidy*nsamplestoconsider)*sizeof(int));
	  
	                 // FIND NUMBER OF CHUNKS AND STORE PAINTING:
	  num_chunksALL=0;
	  for (p=0; p < *ploidy; p++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;   // HAP label
	      if (line_check==0) {printf("Something wrong with %s -- perhaps does not contain %d individuals? Exiting....\n",filenameSAMP,ninds); exit(1);}
	      for (s=0; s < nsamplestoconsider; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,ind_id_vec_new[n][h]+1); exit(1);}
		  step=line;
		  reading(&step,"%s",waste);
		  for (j=0; j < nsites; j++)
		    reading(&step,"%d",&ind_haps_mat[(p*nsamplestoconsider+s)][j]);
		  numchunks = 1;
		  gridweight_val=1;
		  			 
		  for (j=1; j < nsites; j++)
		    {
			// printf("nsites %i %i\n",j, nsites);
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] == ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)] && gridweight_val>(*weight_min)) numchunks = numchunks + 1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=1;
		    }
			

		  num_chunksALL=num_chunksALL+numchunks;
		  numchunks_vec[(p*nsamplestoconsider+s)]=numchunks;
		}
	      for (s=nsamplestoconsider; s < nsamples; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,ind_id_vec_new[n][h]+1); exit(1);}
		}
	    }

	                 // FIND AND STORE ALL CHUNK INFORMATION:
	  double * chunk_posALL = malloc(num_chunksALL*sizeof(double));
	  double * grid_sizeALL = malloc(num_chunksALL*sizeof(double));
	  int * chunk_donor_popALL = malloc(num_chunksALL*sizeof(int));
	  double * gridweightGENDIST_ALL = malloc(num_chunksALL*sizeof(double));
	  double * geneticdiststart_all = malloc(num_chunksALL*sizeof(double));
	  double * geneticdistend_all = malloc(num_chunksALL*sizeof(double));
	
	  chunk_count_all=0;
	  for (p=0; p < *ploidy; p++)
	    {
			
	      for (s=0; s < nsamplestoconsider; s++)
		{
		              // FIND CHUNK LOCATIONS:
		  double * geneticdiststart = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(double));
		  double * geneticdistend = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(double));
		  int * hap_vecGRID = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(int));
		  
		  totalgeneticdist = 0;
		  geneticdiststart[0] = totalgeneticdist;
		  chunkcount = 1;
		  hap_vecGRID[0] = ind_haps_mat[(p*nsamplestoconsider+s)][0];
		  gridweight_val=1;
		  for (j=1; j < nsites; j++)
		    {
		      totalgeneticdist = totalgeneticdist + recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] == ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)] && gridweight_val>(*weight_min))
			{
				
				
			  geneticdiststart[chunkcount]=totalgeneticdist;
			  geneticdistend[(chunkcount-1)]=totalgeneticdist-recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
			  hap_vecGRID[chunkcount] = ind_haps_mat[(p*nsamplestoconsider+s)][j];
	  	        
			 // printf("%f %i\n",geneticdiststart[chunkcount],j);

			  chunkcount = chunkcount + 1;
			}
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=1;
		    }
		  geneticdistend[(chunkcount-1)]=totalgeneticdist;
		  for (i=0; i < numchunks_vec[(p*nsamplestoconsider+s)]; i++)
		    {
			  //printf("%i \n",numchunks_vec[(p*nsamplestoconsider+s)]);
		      chunk_posALL[chunk_count_all]=(geneticdiststart[i]+geneticdistend[i])/2.0;
			  geneticdiststart_all[chunk_count_all] = geneticdiststart[i];
			  geneticdistend_all[chunk_count_all]  = geneticdistend[i];
			  if (donor_label_vec[(hap_vecGRID[i]-1)]==-9)
			{
             
			  //printf("Individual has been painted by missing donor %i, Exiting....\n",hap_vecGRID[i]-1);
			  exit(1);
			} 
		      chunk_donor_popALL[chunk_count_all]=donor_label_vec[(hap_vecGRID[i]-1)];
		  	  
		      gridweightGENDIST_ALL[chunk_count_all]=geneticdistend[i]-geneticdiststart[i];
		      grid_sizeALL[chunk_count_all]=gridweightGENDIST_ALL[chunk_count_all];
		      if (gridweightGENDIST_ALL[chunk_count_all] > (*gridweightMAX_cutoff)) gridweightGENDIST_ALL[chunk_count_all]=(*gridweightMAX_cutoff);
 			  
			  //printf("%i %i %i %f\n",chunk_donor_popALL[chunk_count_all],hap_vecGRID[i], i, gridweightGENDIST_ALL[chunk_count_all]);

			  chunk_count_all=chunk_count_all+1;
			  

		    }
			 //printf("end\n");

			free(geneticdiststart);
			free(geneticdistend);
		  free(hap_vecGRID);
		}
	    }
	
    
		
	 //	printf("before1\n");

	  tabulate_chunks(geneticdistend_all,nsamplestoconsider,ploidy,n,binwidth,num_bins,grid_start,num_donors,num_chunksALL,chunk_donor_popALL,chunk_posALL,grid_sizeALL,gridweightGENDIST_ALL,num_surrogates,total_counts_mat_all);
	 	
	//	printf("end1\n");

	 
	  n=n+1;
	  if (n<(*num_inds_rec)) to_skip=(*ploidy)*(nsamples+1)*(ind_id_vec_new[n][h]-ind_id_vec_new[(n-1)][h]-1);
	  while(to_skip<0)
	    {
	  	  tabulate_chunks(geneticdistend_all,nsamplestoconsider,ploidy,n,binwidth,num_bins,grid_start,num_donors,num_chunksALL,chunk_donor_popALL,chunk_posALL,grid_sizeALL,gridweightGENDIST_ALL,num_surrogates,total_counts_mat_all);
	      n=n+1;
	      if (n==(*num_inds_rec)) break;
	      to_skip=(*ploidy)*(nsamples+1)*(ind_id_vec_new[n][h]-ind_id_vec_new[(n-1)][h]-1);
	    }
	  free(chunk_posALL);
	  free(grid_sizeALL);
	  free(chunk_donor_popALL);
	  free(gridweightGENDIST_ALL);
	  free(geneticdiststart_all);
	  free(geneticdistend_all);
	  for (i=0; i < (*ploidy*nsamplestoconsider); i++)
	    free(ind_haps_mat[i]);
	  free(ind_haps_mat);
	  free(numchunks_vec);
	}
      gzclose(fd2);
      free(posvec);
      free(recomrateperbp);
    }
	

for (int ind_num=0; ind_num<*num_inds_rec;ind_num++)
{		
 for (int  k=0; k < *num_donors; k++)
   {
   for (int h=k; h < *num_donors; h++)
   {			
	   for (int m=0; m < *num_surrogates; m++)
      {
        for ( int n=0; n < *num_surrogates; n++)
  	   {
   		  for (int i=0; i < *num_bins; i++)
			  
		results_all[ind_num][(m*(*num_surrogates)+n)][i]=results_all[ind_num][(m*(*num_surrogates)+n)][i]+temp_weights[(((*num_donors)*k+h)*(*num_surrogates)*(*num_surrogates)+(*num_surrogates)*m+n)]*total_counts_mat_all[ind_num][(k*(*num_donors)+h)][i];   	  }
		  //printf("%f ",temp_weights[(((*num_donors)*k+h)*(*num_surrogates)*(*num_surrogates)+(*num_surrogates)*m+n)]);
	 }
  	}
   } 
  }
	
	

            /* MAKE SYMMETRIC AND DIVIDE BY EXPECTATION, AND SUM ACROSS INDIVIDUALS: */
  double * sum_vec = malloc((*num_bins)*sizeof(double));
  double * exp_ind_res = malloc((*num_surrogates)*(*num_surrogates)*sizeof(double));
  double sum_tot;
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (i=0; i < *num_bins; i++)
	{
	  for (k=0; k < *num_surrogates; k++)
	    {  
                      // !!!!!! IS THIS RIGHT -- SHOULD BE RIGHT FOR DIAGONAL, BUT WHY DIVIDE BY 2 FOR OFF-DIAGONAL?? (THINK IT IS RIGHT, AS IT MEANS EACH OFF-DIAGONAL GETS HALF THE TOTAL -- OTHERWISE DONORS J,K (W/ J!=K) WOULD GET DOUBLE THE WEIGHT -- ONCE WITH J,K AND ONCE WITH K,J) !!!!!!!!!!!!!!!!!!!!!!!!!!
	      for (h=k; h < *num_surrogates; h++) results_all[n][(k*(*num_surrogates)+h)][i]=(results_all[n][(k*(*num_surrogates)+h)][i]+results_all[n][(h*(*num_surrogates)+k)][i])/2.0;
	    }
	  for (k=0; k < (*num_surrogates-1); k++)
	    {  
	      for (h=(k+1); h < (*num_surrogates); h++)
		results_all[n][(h*(*num_surrogates)+k)][i]=results_all[n][(k*(*num_surrogates)+h)][i];
	    }
	  sum_tot=0.0;
	  for (k=0; k < *num_surrogates; k++) sum_vec[k]=0.0;
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		sum_vec[k]=sum_vec[k]+results_all[n][(k*(*num_surrogates)+h)][i];
	      sum_tot=sum_tot+sum_vec[k];
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		exp_ind_res[(k*(*num_surrogates)+h)]=sum_vec[k]*sum_vec[h]/sum_tot;
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=k; h < *num_surrogates; h++)
		  exp_ind_res[(k*(*num_surrogates)+h)]=(exp_ind_res[(k*(*num_surrogates)+h)]+exp_ind_res[(h*(*num_surrogates)+k)])/2.0;
	    }
	  for (k=0; k < (*num_surrogates-1); k++)
	    {  
	      for (h=(k+1); h < (*num_surrogates); h++)
		exp_ind_res[(h*(*num_surrogates)+k)]=exp_ind_res[(k*(*num_surrogates)+h)];
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		{
		  means[(k*(*num_surrogates)*(*num_bins)+h*(*num_bins)+i)]=means[(k*(*num_surrogates)*(*num_bins)+h*(*num_bins)+i)]+(results_all[n][(k*(*num_surrogates)+h)][i]/exp_ind_res[(k*(*num_surrogates)+h)])/(*num_inds_rec);
		}
	    }
	}
    }

  free(line);
  free(filenameSAMP);
  free(filename);
  for (n=0; n < (*num_inds_rec); n++)
    free(ind_id_vec_new[n]);
  free(ind_id_vec_new);
  free(rec_ind_first_hap_index_vec);

  free(sum_vec);
  free(exp_ind_res);
  if (isnan(means[0])){
	  printf("NA values. Exiting....\n");   
	  exit(1);
  }
 
  
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (k=0; k < ((*num_donors) * (*num_donors)); k++)
	free(total_counts_mat_all[n][k]);
      free(total_counts_mat_all[n]);
    }
    free(total_counts_mat_all);
  
  
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
	free(results_all[n][k]);
      free(results_all[n]);
    }
    free(results_all);
	
}

void data_read_null(int *ploidy, char **filenameSAMPread, char **filenameRECOMread, double *binwidth, int *num_bins, int *grid_start, int *num_donors, int * donor_label_vec, int *num_inds_rec, char **recipient_label_vec, double *weight_min, double *gridweightMAX_cutoff, double * chrom_ind_res)
{
  int num_rec_indsMAX=100;   // only use first XX of recipients to generate NULL curves

  int h,i,j,p,s,n,g,m,k;
  int num_chrom, num_chrom2, nsamples, nsamplestoconsider, nsamplestoconsiderNULL, count, nhaps, ninds, nsites, numchunks, to_skip, num_inds_rec_toconsider;
  int chunkcount, line_check, rec_count, hap_count;
  double totalgeneticdist, gridweight_val, gridweight_val2;
  double dist_btwn_i_j,min_dist_btwn_i_j;
  int bin_spot_i_j;
  FILE *fd, *fd3, *fd2;
  char *step;
  char * line = malloc(10000000 * sizeof(char));
  char waste[400];
  char waste2[400];
  char * filename = malloc(1000 * sizeof(char));
  char * filenameSAMP = malloc(1000 * sizeof(char));

  num_inds_rec_toconsider=(*num_inds_rec);
  if (num_inds_rec_toconsider > num_rec_indsMAX) num_inds_rec_toconsider=num_rec_indsMAX;
  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  num_chrom=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom = num_chrom + 1;
    }
  fclose(fd);
  fd = fopen(*filenameRECOMread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
  num_chrom2=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom2 = num_chrom2 + 1;
    }
  fclose(fd);
  if (num_chrom != num_chrom2){printf("The number of files in %s and %s do not match. Exiting....\n",*filenameSAMPread,*filenameRECOMread);exit(1);}

	       // GET NUMBER OF RECIPIENT INDS AND SAMPLES:
  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  strcpy(filenameSAMP,waste);
  fclose(fd);
  fd2 = gzopen(filenameSAMP,"r");
  if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
  line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  while(strcmp(waste,"nsamples")!=0) 
    reading(&step,"%s",waste);
  reading(&step,"%s",waste);
  reading(&step,"%d",&nsamples);
  nsamplestoconsider=nsamples;
  if (nsamplestoconsider<=0) { printf("No painting samples in %s! Exiting....\n",filenameSAMP);exit(1);}
  //nsamplestoconsiderNULL=nsamples;
  nsamplestoconsiderNULL=1;
  count = 0;
  rec_count=0;
  hap_count=0;
  int * rec_ind_first_hap_index_vec = malloc((*num_inds_rec)*sizeof(int));
  while(!gzeof(fd2))
    {
      if(gzgets(fd2,line,10000000)!=NULL);
      step=line;
      reading(&step,"%s",waste);
      if(strcmp(waste,"HAP")==0) 
	{
	  reading(&step,"%s",waste2);
	  reading(&step,"%s",waste2);
	  for (i=0; i < *num_inds_rec; i++)
	    {
	      if (strcmp(waste2,recipient_label_vec[i])==0)
		{
		  if (rec_count>=((*ploidy)*(*num_inds_rec)))
		    {
		      printf("There are more than %d recipient haplotypes in %s! Exiting....\n",(*num_inds_rec)*(*ploidy),filenameSAMP);
		      exit(1);
		    }
		  rec_ind_first_hap_index_vec[(rec_count/(*ploidy))]=hap_count/(*ploidy);
		  rec_count=rec_count+1;
		  break;
		}
	    }
	  hap_count=hap_count+1;
	}
      count = count + 1;
    }
  gzclose(fd2);
  nhaps=count/(nsamples+1);
  if (nhaps != (((double)count)/(nsamples+1)) && nhaps != (((double)count-1)/(nsamples+1)))
    {
      printf("something wrong with file %s. Exiting....\n",filenameSAMP);
      exit(1);
    }
  ninds=nhaps/(*ploidy);
  if (ninds<=0){ printf("No individuals in %s! Exiting....\n",filenameSAMP);exit(1);}
  if (rec_count!=((*num_inds_rec)*(*ploidy)))
    {
      printf("Only %d recipient haplotypes in %s, while expecting %d! Exiting....\n",rec_count,filenameSAMP,(*num_inds_rec)*(*ploidy));
      exit(1);
    }
  for (h=0; h < num_chrom; h++)
    {
      fd = fopen(*filenameRECOMread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if(fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameRECOMread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filename,waste);
      fclose(fd);
      
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      nsites = -1;
      while(!feof(fd3))
	{
	  if(fgets(line,2047,fd3)!=NULL);
	  nsites = nsites + 1;
	}
      fclose(fd3);
      if (nsites<=0) {printf("Something wrong with %s -- no SNPs! Exiting....\n",filename); exit(1);}
      
      double * posvec = malloc(nsites * sizeof(double));
      double * recomrateperbp = malloc(nsites * sizeof(double));
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      for (j=0; j < nsites; j++)
	{
	  line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s at row %d. Exiting....\n",filename,j+1); exit(1);}
	  step=line;
	  reading(&step,"%lf",&posvec[j]);    // basepair position
	  reading(&step,"%lf",&recomrateperbp[j]);
	}
      fclose(fd3);
      for (i=0; i < nsites; i++)
	{
	  if (i > 0)
	    {
	      if (posvec[i]<posvec[(i-1)])
		{
		  printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filename,posvec[(i-1)],posvec[i]);
		  exit(1);
		}
	    }
	  if (recomrateperbp[i]<0)
	    {
	      printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",posvec[i]);
	      exit(1);
	    }
	}
      
      fd = fopen(*filenameSAMPread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filenameSAMP,waste);
      fclose(fd);

      int ** ind_haps_mat=malloc((*ploidy*nsamplestoconsiderNULL*(num_inds_rec_toconsider))*sizeof(int *));
      for (i=0; i < (*ploidy*nsamplestoconsiderNULL*(num_inds_rec_toconsider)); i++)
	ind_haps_mat[i]=malloc(nsites*sizeof(int));
      int * numchunks_vec=malloc((*ploidy*nsamplestoconsiderNULL*(num_inds_rec_toconsider))*sizeof(int));
      fd2 = gzopen(filenameSAMP,"r");
      if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
      line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
      n=0;
      to_skip=(*ploidy)*(nsamples+1)*rec_ind_first_hap_index_vec[n];
      for (n=0; n < (num_inds_rec_toconsider); n++)
	{
	  if (n>0) to_skip=(*ploidy)*(nsamples+1)*(rec_ind_first_hap_index_vec[n]-rec_ind_first_hap_index_vec[(n-1)]-1);
	  for (i=0; i < to_skip; i++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
	      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
	    }
	  
	                 // FIND NUMBER OF CHUNKS AND STORE PAINTING:
	  for (p=0; p < *ploidy; p++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;   // HAP label
	      if (line_check==0) {printf("Something wrong with %s -- perhaps does not contain %d individuals? Exiting....\n",filenameSAMP,ninds); exit(1);}
	      for (s=0; s < nsamplestoconsiderNULL; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,rec_ind_first_hap_index_vec[n]+1); exit(1);}
		  step=line;
		  reading(&step,"%s",waste);
		  for (j=0; j < nsites; j++)
		    reading(&step,"%d",&ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j]);
		  numchunks = 1;
		  gridweight_val=1;
		  for (j=1; j < nsites; j++)
		    {
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] == ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] != ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)] && gridweight_val>(*weight_min)) numchunks = numchunks + 1;
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] != ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)]) gridweight_val=1;
		    }
		  numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]=numchunks;
		}
	      for (s=nsamplestoconsiderNULL; s < nsamples; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,rec_ind_first_hap_index_vec[n]+1); exit(1);}
		}
	    }
	}
      gzclose(fd2);

      for (n=0; n < (num_inds_rec_toconsider-1); n++)
	{
	  for (p=0; p < *ploidy; p++)
	    {
	      for (s=0; s < nsamplestoconsiderNULL; s++)
		{
	                 // FIND CHUNK LOCATIONS:
		  double * geneticdiststart = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(double));
		  double * geneticdistend = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(double));
		  int * hap_vecGRID = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(int));
		  double * chunk_pos = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(double));
		  double * grid_size = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(double));
		  int * chunk_donor_pop = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(int));
		  double * gridweightGENDIST = malloc(numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]*sizeof(double));
		  totalgeneticdist = 0;
		  geneticdiststart[0] = totalgeneticdist;
		  chunkcount = 1;
		  hap_vecGRID[0] = ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][0];
		  gridweight_val=1;
		  for (j=1; j < nsites; j++)
		    {
		      totalgeneticdist = totalgeneticdist + recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] == ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] != ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)] && gridweight_val>(*weight_min))
			{
			  //geneticdiststart[chunkcount]=totalgeneticdist-recomrateperbp[(j-1)]*((posvec[j]-posvec[(j-1)])/2.0)*100;
			  geneticdiststart[chunkcount]=totalgeneticdist;
			  geneticdistend[(chunkcount-1)]=totalgeneticdist-recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
			  hap_vecGRID[chunkcount] = ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j];
			  chunkcount = chunkcount + 1;
			}
		      if (ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][j] != ind_haps_mat[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)][(j-1)]) gridweight_val=1;
		    }
		  geneticdistend[(chunkcount-1)]=totalgeneticdist;
		  for (i=0; i < numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]; i++)
		    {
		      chunk_pos[i]=(geneticdiststart[i]+geneticdistend[i])/2.0;
		      if (donor_label_vec[(hap_vecGRID[i]-1)]==-9)
			{
			  printf("Individual has been painted by missing donor %i, Exiting....\n",donor_label_vec[(hap_vecGRID[i]-1)]);
			  exit(1);
			}
		      chunk_donor_pop[i]=donor_label_vec[(hap_vecGRID[i]-1)];
		      gridweightGENDIST[i]=geneticdistend[i]-geneticdiststart[i];
		      grid_size[i]=gridweightGENDIST[i];
		      if (gridweightGENDIST[i] > (*gridweightMAX_cutoff)) gridweightGENDIST[i]=(*gridweightMAX_cutoff);
		    }
		  for (g=(n+1); g < num_inds_rec_toconsider; g++)
		    {
		      for (m=0; m < (*ploidy); m++)
			{
			  for (k=0; k < nsamplestoconsiderNULL; k++)
			    {
			           // FIND CHUNK LOCATIONS:
			      double * geneticdiststart2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(double));
			      double * geneticdistend2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(double));
			      int * hap2_vecGRID = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(int));
			      double * chunk_pos2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(double));
			      double * grid_size2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(double));
			      int * chunk_donor_pop2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(int));
			      double * gridweightGENDIST2 = malloc(numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]*sizeof(double));
			      totalgeneticdist = 0;
			      geneticdiststart2[0] = totalgeneticdist;
			      chunkcount = 1;
			      hap2_vecGRID[0] = ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][0];
			      gridweight_val2=1;
			      for (j=1; j < nsites; j++)
				{
				  totalgeneticdist = totalgeneticdist + recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
				  if (ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][j] == ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][(j-1)]) gridweight_val2=gridweight_val2+1;
				  if (ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][j] != ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][(j-1)] && gridweight_val2>(*weight_min))
				    {
				      //geneticdiststart[chunkcount]=totalgeneticdist-recomrateperbp[(j-1)]*((posvec[j]-posvec[(j-1)])/2.0)*100;
				      geneticdiststart2[chunkcount]=totalgeneticdist;
				      geneticdistend2[(chunkcount-1)]=totalgeneticdist-recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
				      hap2_vecGRID[chunkcount] = ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][j];
				      chunkcount = chunkcount + 1;
				    }
				  if (ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][j] != ind_haps_mat[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)][(j-1)]) gridweight_val2=1;
				}
			      geneticdistend2[(chunkcount-1)]=totalgeneticdist;
			      for (i=0; i < numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]; i++)
				{
				  chunk_pos2[i]=(geneticdiststart2[i]+geneticdistend2[i])/2.0;
				  if (donor_label_vec[(hap2_vecGRID[i]-1)]==-9)
				    {
			  printf("Individual has been painted by missing donor %i, Exiting....\n",donor_label_vec[(hap_vecGRID[i]-1)]);
				      exit(1);
				    }
				  chunk_donor_pop2[i]=donor_label_vec[(hap2_vecGRID[i]-1)];
				  gridweightGENDIST2[i]=geneticdistend2[i]-geneticdiststart2[i];
				  grid_size2[i]=gridweightGENDIST2[i];
				  if (gridweightGENDIST2[i] > (*gridweightMAX_cutoff)) gridweightGENDIST2[i]=(*gridweightMAX_cutoff);
				}
				       // COMPARE HAP1 CHUNKS TO HAP2 CHUNKS AND TALLY COUNTS:
			      for (i=0; i < numchunks_vec[(n*(*ploidy)*nsamplestoconsiderNULL+p*nsamplestoconsiderNULL+s)]; i++)
				{
				  for (j=0; j < numchunks_vec[(g*(*ploidy)*nsamplestoconsiderNULL+m*nsamplestoconsiderNULL+k)]; j++)
				    {
				      dist_btwn_i_j=pow(pow((chunk_pos[i]-chunk_pos2[j]),2.0),0.5);
				      if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) <= ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)));
				      if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) > ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)))+1;
				      bin_spot_i_j = bin_spot_i_j - (*grid_start);
				      min_dist_btwn_i_j=grid_size[i]/2.0+grid_size2[j]/2.0;
				      if (bin_spot_i_j >= 0 && bin_spot_i_j<(*num_bins) && dist_btwn_i_j >= min_dist_btwn_i_j && chunk_donor_pop[i]>0 && chunk_donor_pop2[j]>0) chrom_ind_res[((chunk_donor_pop[i]-1)*(*num_donors)*(*num_bins)+(chunk_donor_pop2[j]-1)*(*num_bins)+bin_spot_i_j)]=chrom_ind_res[((chunk_donor_pop[i]-1)*(*num_donors)*(*num_bins)+(chunk_donor_pop2[j]-1)*(*num_bins)+bin_spot_i_j)]+gridweightGENDIST[i]*gridweightGENDIST2[j];
				    }
				}
			      free(geneticdiststart2);
			      free(geneticdistend2);
			      free(hap2_vecGRID);
			      free(chunk_pos2);
			      free(grid_size2);
			      free(gridweightGENDIST2);
			      free(chunk_donor_pop2);
			    }
			}
		    }
		  free(geneticdiststart);
		  free(geneticdistend);
		  free(hap_vecGRID);
		  free(chunk_pos);
		  free(grid_size);
		  free(gridweightGENDIST);
		  free(chunk_donor_pop);
		}
	    }
	}
      for (i=0; i < (*ploidy*nsamplestoconsiderNULL*(num_inds_rec_toconsider)); i++)
	free(ind_haps_mat[i]);
      free(ind_haps_mat);
      free(numchunks_vec);
      free(posvec);
      free(recomrateperbp);
    }

  free(line);
  free(filenameSAMP);
  free(filename);
  free(rec_ind_first_hap_index_vec);
  if (isnan(chrom_ind_res[0])){printf("NA values. Exiting....\n"); exit(1);}
}


void tabulate_chunks_mode3(double * geneticdistend_all,int nsamplestoconsider,int *ploidy,int ind_num, double *binwidth, int *num_bins, int *grid_start, int *num_donors, int num_chunksALL, int * chunk_donor_popALL, double * chunk_posALL, double * grid_sizeALL, double * gridweight_GENDIST_ALL, int *num_surrogates, double * temp_weights, double *** results)
{   
    int count_ans, times, chunk1,chunk2;
    int h,i,j,k,m,n,c;
    c = 0;
    int  bin_spot_i_j;
    double dist_btwn_i_j,min_dist_btwn_i_j;
    double ** total_counts_mat = malloc((*num_donors)*(*num_donors)*sizeof(double *));
    for (k=0; k < (*num_donors * (*num_donors)); k++)
      total_counts_mat[k]=malloc((*num_bins)*sizeof(double));
      
    for (k=0; k < (*num_donors * (*num_donors)); k++)
      {
        for (i=0; i < *num_bins; i++)
    	total_counts_mat[k][i]=0.0;
      }
       
	 
      double chr_length = 0.0;
	
  	   for (i=0; i <num_chunksALL; i++){
          if (geneticdistend_all[i] >= chr_length)
  		{
  			chr_length = geneticdistend_all[i];
  		}
			
         }
	  int bin_size_const = 1; 
	  int num_bin = (int)chr_length/bin_size_const+1;
  
       int *bin = malloc(num_chunksALL*sizeof(int));
       int *numchunckinbin = malloc((num_bin)*sizeof(int));
       for (i=0; i <num_bin; i++)
       numchunckinbin[i] =0;
	   
   	   for (i=0; i <num_chunksALL; i++){
	   bin[i] = (int)chunk_posALL[i];
       }
  	   for (i=0; i <num_chunksALL; i++)
  	   { numchunckinbin[bin[i]] =  numchunckinbin[bin[i]]+1;
	   }
       int *start_bin = malloc((num_bin)*sizeof(int));
	   start_bin[0] = 0;
	   for (i=1; i <num_bin; i++)
	   {
				start_bin[i] = start_bin[i-1] + numchunckinbin[i-1];
		} 
	 // tabulate new chunk information  
   	double *tabulate_chunk_posALL = malloc((num_chunksALL)*sizeof(double));
   	double *tabulate_grid_sizeALL = malloc((num_chunksALL)*sizeof(double));
   	double *tabulate_gridweight_GENDIST_ALL = malloc((num_chunksALL)*sizeof(double));
   	int *tabulate_chunk_donor_popALL = malloc((num_chunksALL)*sizeof(int));
    int cc = 0;
   	for ( j=0; j <num_bin; j++)
   	{
   		for ( k=0; k <num_chunksALL; k++){
   		if (bin[k] == j){
   		tabulate_chunk_posALL[cc] = chunk_posALL[k];
   		tabulate_grid_sizeALL[cc] = grid_sizeALL[k];
   		tabulate_gridweight_GENDIST_ALL[cc] = gridweight_GENDIST_ALL[k];
   		tabulate_chunk_donor_popALL[cc] = chunk_donor_popALL[k];
    	    cc = cc+1;
   		}
   	   }
   	}

    int window_size = ceil(*binwidth*(*num_bins))+1;
      	count_ans = 0;
      	int t_1, t_2, j_lim;
      	for ( j=0; j <num_bin-1; j++)
      	{
      	j_lim =j+window_size;
      	if (j+window_size > num_bin){j_lim =num_bin;}
          for ( k=j+1; k <j_lim; k++)
      	{
      		int Y_i;
			 t_1 = numchunckinbin[j];
			 t_2 = numchunckinbin[k];
			Y_i = (int)t_1*t_2*exp(-0.05*(k-j))/8;
		    if (t_1>0 && t_2>0 &&  Y_i == 0) {Y_i = 1;}
		
      		count_ans = count_ans+Y_i;
      		
      	    times = 0;		
			
      		while(times < Y_i) 
			
      		{
         		  chunk1 = start_bin[j] + rand() %  (t_1); 
				  chunk2 = start_bin[k] + rand() %  (t_2); 
				 
	  		 dist_btwn_i_j=pow(pow((tabulate_chunk_posALL[chunk1]-tabulate_chunk_posALL[chunk2]),2.0),0.5);
      		 if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) <= ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)));
      		 if (((dist_btwn_i_j/(*binwidth)-floor(dist_btwn_i_j/(*binwidth)))*(*binwidth)) > ((*binwidth)/2.0)) bin_spot_i_j=((int)floor(dist_btwn_i_j/(*binwidth)))+1;
	 
      		 bin_spot_i_j = bin_spot_i_j - (*grid_start);
	 
      		 min_dist_btwn_i_j=tabulate_grid_sizeALL[chunk1]/2.0+tabulate_grid_sizeALL[chunk2]/2.0;
      	     if (bin_spot_i_j >= 0 && bin_spot_i_j<(*num_bins) && dist_btwn_i_j >= min_dist_btwn_i_j) 
      			 {
      			 if ((tabulate_chunk_donor_popALL[chunk1]-1) > (tabulate_chunk_donor_popALL[chunk2]-1))
      			 {
      			  total_counts_mat[((tabulate_chunk_donor_popALL[chunk2]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk1]-1))][bin_spot_i_j]=total_counts_mat[((tabulate_chunk_donor_popALL[chunk2]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk1]-1))][bin_spot_i_j]+tabulate_gridweight_GENDIST_ALL[chunk2]*tabulate_gridweight_GENDIST_ALL[chunk1];
      			 }	
      		     else 
      			 {
			 
      			  total_counts_mat[((tabulate_chunk_donor_popALL[chunk1]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk2]-1))][bin_spot_i_j]=total_counts_mat[((tabulate_chunk_donor_popALL[chunk1]-1)*(*num_donors)+(tabulate_chunk_donor_popALL[chunk2]-1))][bin_spot_i_j]+tabulate_gridweight_GENDIST_ALL[chunk1]*tabulate_gridweight_GENDIST_ALL[chunk2];
      		     }
      		 
          
      			  }
	 
      	 	times = times+1;
	    
      	    }
             
          }
          }

	 
	
 
  free(numchunckinbin);
  free(bin);
  free(start_bin);
  free(tabulate_chunk_donor_popALL);
  free(tabulate_grid_sizeALL);
  free(tabulate_gridweight_GENDIST_ALL);
  free(tabulate_chunk_posALL);		   
 
  for (m=0; m < *num_surrogates; m++)
    {
      for (n=0; n < *num_surrogates; n++)
	{
	  for (k=0; k < *num_donors; k++)
	    {
	      for (h=k; h < *num_donors; h++)
		{
		  for (i=0; i < *num_bins; i++)
		    results[ind_num][(m*(*num_surrogates)+n)][i]=results[ind_num][(m*(*num_surrogates)+n)][i]+temp_weights[(((*num_donors)*k+h)*(*num_surrogates)*(*num_surrogates)+(*num_surrogates)*m+n)]*total_counts_mat[(k*(*num_donors)+h)][i];
		}
	    }
	}
    }
   	
  for (k=0; k < (*num_donors * (*num_donors)); k++)
    free(total_counts_mat[k]);
  free(total_counts_mat);
    
}

void data_read_mode3(int *ploidy, int * ind_id_vec, char **filenameSAMPread, char **filenameRECOMread, double *binwidth, int *num_bins, int *grid_start, int *num_donors, int * donor_label_vec, int *num_inds_rec, char **recipient_label_vec, double *weight_min, double *gridweightMAX_cutoff, int * num_surrogates, double * temp_weights, double * means)
{
  int h,i,j,p,s,n,k;
  int num_chrom, num_chrom2, nsamples, nsamplestoconsider, nsamplestoconsiderNULL, count, nhaps, ninds, nsites, numchunks, num_chunksALL, to_skip;
  int chunkcount, chunk_count_all, line_check, rec_count, hap_count;
  double totalgeneticdist, gridweight_val;
  FILE *fd, *fd3, *fd2;
  //gzFile *fd2;
  char *step;
  char * line = malloc(10000000 * sizeof(char));
  char waste[400];
  char waste2[400];
  char * filename = malloc(1000 * sizeof(char));
  char * filenameSAMP = malloc(1000 * sizeof(char));

  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  num_chrom=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom = num_chrom + 1;
    }
  fclose(fd);
  fd = fopen(*filenameRECOMread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
  num_chrom2=-1;
  while(!feof(fd))
    {
      if (fgets(line,10000000,fd)!=NULL);
      num_chrom2 = num_chrom2 + 1;
    }
  fclose(fd);
  if (num_chrom != num_chrom2){printf("The number of files in %s and %s do not match. Exiting....\n",*filenameSAMPread,*filenameRECOMread);exit(1);}

	       // GET NUMBER OF RECIPIENT INDS AND SAMPLES:
  fd = fopen(*filenameSAMPread,"r");
  if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  strcpy(filenameSAMP,waste);
  fclose(fd);
  fd2 = gzopen(filenameSAMP,"r");
  if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
  line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
  step=line;
  reading(&step,"%s",waste);
  while(strcmp(waste,"nsamples")!=0) 
    reading(&step,"%s",waste);
  reading(&step,"%s",waste);
  reading(&step,"%d",&nsamples);
  nsamplestoconsider=nsamples;
  if (nsamplestoconsider<=0) { printf("No painting samples in %s! Exiting....\n",filenameSAMP);exit(1);}
  //nsamplestoconsiderNULL=nsamples;
  nsamplestoconsiderNULL=1;
  count = 0;
  rec_count=0;
  hap_count=0;
  int * rec_ind_first_hap_index_vec = malloc((*num_inds_rec)*sizeof(int));
  while(!gzeof(fd2))
    {
      if(gzgets(fd2,line,10000000)!=NULL);
      step=line;
      reading(&step,"%s",waste);
      if(strcmp(waste,"HAP")==0) 
	{
	  reading(&step,"%s",waste2);
	  reading(&step,"%s",waste2);
	  for (i=0; i < *num_inds_rec; i++)
	    {
	      if (strcmp(waste2,recipient_label_vec[i])==0)
		{
		  if (rec_count>=((*ploidy)*(*num_inds_rec)))
		    {
		      printf("There are more than %d recipient haplotypes in %s! Exiting....\n",(*num_inds_rec)*(*ploidy),filenameSAMP);
		      exit(1);
		    }
		  rec_ind_first_hap_index_vec[(rec_count/(*ploidy))]=hap_count/(*ploidy);
		  rec_count=rec_count+1;
		  break;
		}
	    }
	  hap_count=hap_count+1;
	}
      count = count + 1;
    }
  gzclose(fd2);
  int ** ind_id_vec_new = malloc((*num_inds_rec)*sizeof(int*));
  for (i=0; i < (*num_inds_rec); i++)
    ind_id_vec_new[i]=malloc(num_chrom*sizeof(int));
  for (i=0; i < (*num_inds_rec); i++)
    {
      for (h=0; h < num_chrom; h++)
	ind_id_vec_new[i][h]=rec_ind_first_hap_index_vec[ind_id_vec[(i*num_chrom+h)]];
    }
  nhaps=count/(nsamples+1);
  if (nhaps != (((double)count)/(nsamples+1)) && nhaps != (((double)count-1)/(nsamples+1)))
    { 
      printf("something wrong with file %s. Exiting....\n",filenameSAMP);
      exit(1);
    }
  ninds=nhaps/(*ploidy);
  if (ninds<=0){ printf("No individuals in %s! Exiting....\n",filenameSAMP);exit(1);}
  if (rec_count!=((*num_inds_rec)*(*ploidy)))
    {
      printf("Only %d recipient haplotypes in %s, while expecting %d! Exiting....\n",rec_count,filenameSAMP,(*num_inds_rec)*(*ploidy));
      exit(1);
    }
  double *** results = malloc((*num_inds_rec)*sizeof(double **));
  for (n=0; n < (*num_inds_rec); n++)
    {
      results[n]=malloc((*num_surrogates) * (*num_surrogates) * sizeof(double *));
      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
	{
	  results[n][k]=malloc((*num_bins)*sizeof(double));
	}
    }
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
	{
	  for (i=0; i < (*num_bins); i++)
	    results[n][k][i]=0.0;
	}
    }
  for (h=0; h < num_chrom; h++)
    {
                       /* GET NUMBER OF SNPs: */
      fd = fopen(*filenameRECOMread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameRECOMread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if(fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameRECOMread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filename,waste);
      fclose(fd);
      
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      nsites = -1;
      while(!feof(fd3))
	{
	  if(fgets(line,2047,fd3)!=NULL);
	  nsites = nsites + 1;
	}
      fclose(fd3);
      if (nsites<=0) {printf("Something wrong with %s -- no SNPs! Exiting....\n",filename); exit(1);}
	  
                       /* GET RECOMBINATION RATES AND SNP POSITIONS: */
      double * posvec = malloc(nsites * sizeof(double));
      double * recomrateperbp = malloc(nsites * sizeof(double));
      fd3 = fopen(filename,"r");
      if (fd3 == NULL) { printf("error opening %s\n",filename); exit(1);}
      line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;   // header
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filename); exit(1);}
      for (j=0; j < nsites; j++)
	{
	  line_check=0; if(fgets(line,2047,fd3)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s at row %d. Exiting....\n",filename,j+1); exit(1);}
	  step=line;
	  reading(&step,"%lf",&posvec[j]);    // basepair position
	  reading(&step,"%lf",&recomrateperbp[j]);
	}
      fclose(fd3);
      for (i=0; i < nsites; i++)
	{
	  if (i > 0)
	    {
	      if (posvec[i]<posvec[(i-1)])
		{
		  printf("positions in %s are not listed in increasing order at basepairs %lf and %lf. Exiting....\n",filename,posvec[(i-1)],posvec[i]);
		  exit(1);
		}
	    }
	  if (recomrateperbp[i]<0)
	    {
	      printf("recom rate must be > 0 (basepair %lf)!! Exiting....\n",posvec[i]);
	      exit(1);
	    }
	}
      
      fd = fopen(*filenameSAMPread,"r");
      if (fd == NULL) { printf("error opening %s\n",*filenameSAMPread); exit(1);}
      for (j=0; j <= h; j++) 
	{
	  line_check=0; if (fgets(line,1000000,fd)!=NULL) line_check=1;
	  if (line_check==0) {printf("Something wrong with %s. Exiting....\n",*filenameSAMPread); exit(1);}
	}
      step=line;
      reading(&step,"%s",waste);
      strcpy(filenameSAMP,waste);
      fclose(fd);

               /* IF NOT DOING NULL IND, CONSIDER COMPARISONS BETWEEN ALL CHUNKS: */
                       /* FIND NUMBER OF CHUNKS ACROSS ALL SAMPLES FOR IND: */
      fd2 = gzopen(filenameSAMP,"r");
      if (fd2 == NULL) { printf("error opening %s\n",filenameSAMP); exit(1);}
      line_check=0; if (gzgets(fd2,line,2047)!=NULL) line_check=1;
      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
      n=0;
      to_skip=(*ploidy)*(nsamples+1)*ind_id_vec_new[n][h];
      while (n < (*num_inds_rec))
	{
	  for (i=0; i < to_skip; i++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
	      if (line_check==0) {printf("Something wrong with %s. Exiting....\n",filenameSAMP); exit(1);}
	    }
	  
	  int ** ind_haps_mat=malloc((*ploidy*nsamplestoconsider)*sizeof(int *));
	  for (i=0; i < (*ploidy*nsamplestoconsider); i++)
	    ind_haps_mat[i]=malloc(nsites*sizeof(int));
	  int * numchunks_vec=malloc((*ploidy*nsamplestoconsider)*sizeof(int));
	  
	                 // FIND NUMBER OF CHUNKS AND STORE PAINTING:
	  num_chunksALL=0;
	  for (p=0; p < *ploidy; p++)
	    {
	      line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;   // HAP label
	      if (line_check==0) {printf("Something wrong with %s -- perhaps does not contain %d individuals? Exiting....\n",filenameSAMP,ninds); exit(1);}
	      for (s=0; s < nsamplestoconsider; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,ind_id_vec_new[n][h]+1); exit(1);}
		  step=line;
		  reading(&step,"%s",waste);
		  for (j=0; j < nsites; j++)
		    reading(&step,"%d",&ind_haps_mat[(p*nsamplestoconsider+s)][j]);
		  numchunks = 1;
		  gridweight_val=1;
		  for (j=1; j < nsites; j++)
		    {
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] == ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)] && gridweight_val>(*weight_min)) numchunks = numchunks + 1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=1;
		    }
		  num_chunksALL=num_chunksALL+numchunks;
		  numchunks_vec[(p*nsamplestoconsider+s)]=numchunks;
		}
	      for (s=nsamplestoconsider; s < nsamples; s++)
		{
		  line_check=0; if (gzgets(fd2,line,10000000)!=NULL) line_check=1;
		  if (line_check==0) {printf("Something wrong with %s at sample %d of individual %d. Exiting....\n",filenameSAMP,s+1,ind_id_vec_new[n][h]+1); exit(1);}
		}
	    }

	                 // FIND AND STORE ALL CHUNK INFORMATION:
	  double * chunk_posALL = malloc(num_chunksALL*sizeof(double));
	  double * grid_sizeALL = malloc(num_chunksALL*sizeof(double));
	  int * chunk_donor_popALL = malloc(num_chunksALL*sizeof(int));
	  double * gridweightGENDIST_ALL = malloc(num_chunksALL*sizeof(double));
	  double * geneticdiststart_all = malloc(num_chunksALL*sizeof(double));
	  double * geneticdistend_all = malloc(num_chunksALL*sizeof(double));
		  
	  chunk_count_all=0;
	  for (p=0; p < *ploidy; p++)
	    {
	      for (s=0; s < nsamplestoconsider; s++)
		{
		              // FIND CHUNK LOCATIONS:
		  double * geneticdiststart = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(double));
		  double * geneticdistend = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(double));
		  int * hap_vecGRID = malloc(numchunks_vec[(p*nsamplestoconsider+s)]*sizeof(int));
		  
		  totalgeneticdist = 0;
		  geneticdiststart[0] = totalgeneticdist;
		  chunkcount = 1;
		  hap_vecGRID[0] = ind_haps_mat[(p*nsamplestoconsider+s)][0];
		  gridweight_val=1;
		  for (j=1; j < nsites; j++)
		    {
		      totalgeneticdist = totalgeneticdist + recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] == ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=gridweight_val+1;
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)] && gridweight_val>(*weight_min))
			{
			  geneticdiststart[chunkcount]=totalgeneticdist;
			  geneticdistend[(chunkcount-1)]=totalgeneticdist-recomrateperbp[(j-1)]*(posvec[j]-posvec[(j-1)])*100;
			  hap_vecGRID[chunkcount] = ind_haps_mat[(p*nsamplestoconsider+s)][j];
	  	        
			  
			  chunkcount = chunkcount + 1;
			}
		      if (ind_haps_mat[(p*nsamplestoconsider+s)][j] != ind_haps_mat[(p*nsamplestoconsider+s)][(j-1)]) gridweight_val=1;
		    }
		  geneticdistend[(chunkcount-1)]=totalgeneticdist;
		  for (i=0; i < numchunks_vec[(p*nsamplestoconsider+s)]; i++)
		    {
		      chunk_posALL[chunk_count_all]=(geneticdiststart[i]+geneticdistend[i])/2.0;
			  geneticdiststart_all[chunk_count_all] = geneticdiststart[i];
			  geneticdistend_all[chunk_count_all]  = geneticdistend[i];
			  if (donor_label_vec[(hap_vecGRID[i]-1)]==-9)
			{
			  printf("Individual has been painted by missing donor %i, Exiting....\n",donor_label_vec[(hap_vecGRID[i]-1)]);
			  exit(1);
			}
		      chunk_donor_popALL[chunk_count_all]=donor_label_vec[(hap_vecGRID[i]-1)];
		      gridweightGENDIST_ALL[chunk_count_all]=geneticdistend[i]-geneticdiststart[i];
		      grid_sizeALL[chunk_count_all]=gridweightGENDIST_ALL[chunk_count_all];
		      if (gridweightGENDIST_ALL[chunk_count_all] > (*gridweightMAX_cutoff)) gridweightGENDIST_ALL[chunk_count_all]=(*gridweightMAX_cutoff);
 		      chunk_count_all=chunk_count_all+1;
		    }
			
			free(geneticdiststart);
			free(geneticdistend);
		  free(hap_vecGRID);
		}
	    }
	     
	  tabulate_chunks_mode3(geneticdistend_all,nsamplestoconsider,ploidy,n,binwidth,num_bins,grid_start,num_donors,num_chunksALL,chunk_donor_popALL,chunk_posALL,grid_sizeALL,gridweightGENDIST_ALL,num_surrogates,temp_weights,results);
	  n=n+1;
	  if (n<(*num_inds_rec)) to_skip=(*ploidy)*(nsamples+1)*(ind_id_vec_new[n][h]-ind_id_vec_new[(n-1)][h]-1);
	  while(to_skip<0)
	    {
	  	  tabulate_chunks_mode3(geneticdistend_all,nsamplestoconsider,ploidy,n,binwidth,num_bins,grid_start,num_donors,num_chunksALL,chunk_donor_popALL,chunk_posALL,grid_sizeALL,gridweightGENDIST_ALL,num_surrogates,temp_weights,results);
	      n=n+1;
	      if (n==(*num_inds_rec)) break;
	      to_skip=(*ploidy)*(nsamples+1)*(ind_id_vec_new[n][h]-ind_id_vec_new[(n-1)][h]-1);
	    }
	  free(chunk_posALL);
	  free(grid_sizeALL);
	  free(chunk_donor_popALL);
	  free(gridweightGENDIST_ALL);
	  free(geneticdiststart_all);
	  free(geneticdistend_all);
	  for (i=0; i < (*ploidy*nsamplestoconsider); i++)
	    free(ind_haps_mat[i]);
	  free(ind_haps_mat);
	  free(numchunks_vec);
	}
      gzclose(fd2);
      free(posvec);
      free(recomrateperbp);
    }

            /* MAKE SYMMETRIC AND DIVIDE BY EXPECTATION, AND SUM ACROSS INDIVIDUALS: */
  double * sum_vec = malloc((*num_bins)*sizeof(double));
  double * exp_ind_res = malloc((*num_surrogates)*(*num_surrogates)*sizeof(double));
  double sum_tot;
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (i=0; i < *num_bins; i++)
	{
	  for (k=0; k < *num_surrogates; k++)
	    {  
                      // !!!!!! IS THIS RIGHT -- SHOULD BE RIGHT FOR DIAGONAL, BUT WHY DIVIDE BY 2 FOR OFF-DIAGONAL?? (THINK IT IS RIGHT, AS IT MEANS EACH OFF-DIAGONAL GETS HALF THE TOTAL -- OTHERWISE DONORS J,K (W/ J!=K) WOULD GET DOUBLE THE WEIGHT -- ONCE WITH J,K AND ONCE WITH K,J) !!!!!!!!!!!!!!!!!!!!!!!!!!
	      for (h=k; h < *num_surrogates; h++) results[n][(k*(*num_surrogates)+h)][i]=(results[n][(k*(*num_surrogates)+h)][i]+results[n][(h*(*num_surrogates)+k)][i])/2.0;
	    }
	  for (k=0; k < (*num_surrogates-1); k++)
	    {  
	      for (h=(k+1); h < (*num_surrogates); h++)
		results[n][(h*(*num_surrogates)+k)][i]=results[n][(k*(*num_surrogates)+h)][i];
	    }
	  sum_tot=0.0;
	  for (k=0; k < *num_surrogates; k++) sum_vec[k]=0.0;
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		sum_vec[k]=sum_vec[k]+results[n][(k*(*num_surrogates)+h)][i];
	      sum_tot=sum_tot+sum_vec[k];
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		exp_ind_res[(k*(*num_surrogates)+h)]=sum_vec[k]*sum_vec[h]/sum_tot;
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=k; h < *num_surrogates; h++)
		  exp_ind_res[(k*(*num_surrogates)+h)]=(exp_ind_res[(k*(*num_surrogates)+h)]+exp_ind_res[(h*(*num_surrogates)+k)])/2.0;
	    }
	  for (k=0; k < (*num_surrogates-1); k++)
	    {  
	      for (h=(k+1); h < (*num_surrogates); h++)
		exp_ind_res[(h*(*num_surrogates)+k)]=exp_ind_res[(k*(*num_surrogates)+h)];
	    }
	  for (k=0; k < *num_surrogates; k++)
	    {  
	      for (h=0; h < *num_surrogates; h++) 
		{
		  means[(k*(*num_surrogates)*(*num_bins)+h*(*num_bins)+i)]=means[(k*(*num_surrogates)*(*num_bins)+h*(*num_bins)+i)]+(results[n][(k*(*num_surrogates)+h)][i]/exp_ind_res[(k*(*num_surrogates)+h)])/(*num_inds_rec);
		}
	    }
	}
    }

  free(line);
  free(filenameSAMP);
  free(filename);
  for (n=0; n < (*num_inds_rec); n++)
    free(ind_id_vec_new[n]);
  free(ind_id_vec_new);
  free(rec_ind_first_hap_index_vec);
  for (n=0; n < (*num_inds_rec); n++)
    {
      for (k=0; k < ((*num_surrogates) * (*num_surrogates)); k++)
	free(results[n][k]);
      free(results[n]);
    }
  free(results);
  free(sum_vec);
  free(exp_ind_res);
  if (isnan(means[0])){printf("NA values. Exiting....\n"); exit(1);}
}


