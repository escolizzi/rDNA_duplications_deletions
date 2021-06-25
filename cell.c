#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#include <ctype.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "cell.h"
#include "my.h"
#include "mersenne.h"
#include "constants.h"
#include "parameters.h"


// REPLICATION, 
// new function, everyone could be ready to replicate
// fitness inv. proportional to ic->lifespan
// of course keep constraints (although I wonder what they could do without them)

//|---   find some decent fitness function with both ic->health and ic->lifespan: 
//         lifespan*(2-Q): approaches lifespan when Q ~> 1, doubles it when Q is 0	
//| -- because you'll know lifespan, you can return max integration time 10x average for next time step    
    
double Replication(CELL *world,double **y,int cells_status,double max_growthtime)
{
  int i,j,howmany=0,howmanymay=0,howmanynot=0,rpos;
  double this_fitness,tot_fitness=0.,fitness[par_pop];
  int idcell[par_pop],idcellmay[par_pop],idcellnot[par_pop];
  
  double maxlifespan=-1.;
  double delta_step=1./par_init_interval_time;
  
  //CELL *slowest;
  //int slowguy;
  
  // ---
  // 1) -- We screen the population, to see who is fit and who's not, and we assign fitness
  //       'howmany' counts those that can replicate
  //       'howmanymay' counts those that can survive but not replicate
  //       'howmanynot' are those with zero fitness and can be replaced
  // ---
  for(i=0;i<par_pop;i++){
    
    int this_celldiv=world[i].celldiv;
    double this_health=world[i].health;
    
    if(this_celldiv==1 && this_health>0){
      
      idcell[howmany]=i;			//we get the position of those that can divide
      double this_lifespan=world[i].lifespan + delta_step;
      
      //suppose a guy got a big deletion, 
      // but is stil viable and half of his volume is larger than tarvol:
      // you REALLY want this guy to replicate, because he is VERY fit.
      this_fitness = world[i].health * max_growthtime/(this_lifespan);
      fitness[howmany] = this_fitness;
      
      tot_fitness+=this_fitness;
      howmany++;
      
      // this is just statistics, it is the guy that took longest to divide, but made it!
      if(this_lifespan>maxlifespan){
	maxlifespan=this_lifespan;
	//slowest=&world[i];
	//slowguy=i;
      }
      
    }else if(this_health>0 && world[i].cellvol > 0.5*world[i].tarvol ){
      idcellmay[howmanymay]=i;			// these are individuals that grew but not reach tar vol, so -MAY- survive, but not divide
      howmanymay++;
    }else{
      // these guys shrank during lifetime, and will be eliminated
      idcellnot[howmanynot]=i;
      howmanynot++;
    }
  }
  
  if(howmany==0) return -1.;	//because howmany counts only the fit guys, this should exclude also the case when tot_fitness==0
  
  //Now, because a shorter interdivision time (lifespan) should be rewarded with higher fitness
  // we loop through idcell[] and fitness[] to reflect this. 
  //In particular, if your lifespan was half, your p_repl should be double
  
  // ---
  // 2) -- We find positions in the population array that are going to be overwritten because they are not fit
  // ---
  
  //how many are going to be replaced? 
  // If more people have not reached target volume than par_pop*death
  // then, 'howmanynot' are the people who are going to be replaced
  // else it is going to be par_pop*death, which includes 
  //  those that have not reached tartget vol but also some that have
  //  proportional to fitness.
  
  int emptyspace =(int)(0.5+(double)(par_pop)*par_death);
  if(emptyspace==0) emptyspace++;
  
  if(howmanynot>=emptyspace) emptyspace=howmanynot;
  else{
    //we randomly draw from the rest of the population (i.e. idcell[] and idcellmay[])
    // enough individuals to reach par_pop*death
    for(i=howmanynot;i<emptyspace;i++){
      //take a random individual from idcell or idcellmay
      rpos=(int)((howmany+howmanymay)*genrand_real2());
      //pass it to idcellnot
      if(rpos<howmany){
	idcellnot[i]=idcell[rpos];
	//take the last position of idcell and put it in the position of the guy just drawn
	idcell[rpos]=idcell[howmany-1];
	//decrease size of idcell (howmany) by one, (increase howmanynot is not necessary, 
        // since we are going to have howmanynot=par_pop*death;
	howmany--;
      }else{
	rpos-=howmany; //we draw it from idcellmay, at position rpos-howmany
	idcellnot[i]=idcellmay[rpos];
	idcellmay[rpos]=idcellmay[howmanymay-1];
	howmanymay--;
      }
      //repeat
    }
  }
    
  // ---
  // 3) -- We duplicate and mutate fit individuals into positions of unfit ones
  // ---
  
  // Who will replicate in the positions specified by idcellnot[]? 
  // it depends on fitness, if your fitness is zero, then you won't
  
  //notice that this works only if tot_fitness>0
  // but
  for(i=0;i<emptyspace;i++){
    
    double where=tot_fitness*genrand_real2();
    double here=0.;
    for(j=0;j<howmany;j++){
      here += fitness[j];
      if(here>where){
        break; // now idcell[j] has the position of a fit guy that will replicate
      }
    }
    
    //Copy the guy at pos idcell[j] into the position specified by idcellnot[i], i.e into world[idcellnot[i]]!
    CopyCellHereToThere(world,y, idcell[j], idcellnot[i] );
    //fprintf(stderr,"here: %d, there: %d\n", idcell[j], idcellnot[i]);
    //Mutate him
    MutateCell( &world[idcellnot[i]] );
    
  }
  
  // ---
  // 4) -- Whoever is alive gets its volume halved
  // ---
  
  //Divide everyone's volume and cell content (except S), 
  // set volume and transcriptional load to zero
  for(i=0;i<par_pop;i++) HalfCellContent(world, y, i);
  
  return maxlifespan;
}

// makes an exact copy of a cell from parent_pos to child_pos
// re-set the pointers of child_pos to point to the correct parts of y
void CopyCellHereToThere(CELL *world,double **y, int parent_pos, int child_pos)
{
  int i;
  CELL *parent,*child;
  
  if(child_pos==parent_pos){
    //fprintf(stderr,"Warning, parent pos and child pos coincide\nNothing to copy, then.\n");
    return;
  }
  
  parent = &world[parent_pos];
  child  = &world[child_pos];
  
  //copies variables from parent to child
  for(i=0;i<N_VARIABLES;i++) y[child_pos][i] = y[parent_pos][i];
  
  //child takes values of parent
  *child = *parent;
  
  //but because I really thought this through in its inception, 
  // child pointers to variables are now pointing to parent's variables, 
  //so I have to reset all pointers to child_pos
  child->S=&(y[child_pos][_S]);
  child->A=&(y[child_pos][_A]);
  
  child->T=&(y[child_pos][_T]);
  child->Q=&(y[child_pos][_Q]);
  child->Rr=&(y[child_pos][_Rr]);
  child->Rp=&(y[child_pos][_Rp]);
  
  child->mQ=&(y[child_pos][_mQ]);
  child->mT=&(y[child_pos][_mT]);
  child->mRp=&(y[child_pos][_mRp]);
  
//   child->Crmt=&(y[child_pos][_Crmt]);
//   child->Crmq=&(y[child_pos][_Crmq]);
  
  strcpy(child->genome,parent->genome);
}

void HalfCellContent(CELL *world,double **y, int pos)
{
  int i;
  CELL *ic;
  
  //half volumes of all cell variables, except S
  for(i=1;i<N_VARIABLES;i++) y[pos][i] /= 2.;
  
  ic = &world[pos];
  
  //set transcr load to zero
  for(i=0;i<N_GENE_TYPES;i++){
    ic->trload[i]=0.;
  }
  
  //updates tags and stuff of parent and child
  CellVolume(ic);
  CellHealth(ic);
  ic->celldiv=0;
  
}

void MutateCell(CELL *child)
{
  //  MUTATIONS
  // MutationsMetabolism(child);	WE DO NOT MUTATE METABOLIC GENES ANY MORE
  
  if(MUT_REGULATION) MutationsRegulation(child);
  
  if(MUT_TRANSCR_LOAD) TranscrMutations_diffratesDupDelIn(child);	//these mutations occur during division, and only on the active strand, 
					                            // so it is ok to assign mutations only to child... maybe it should be to parent?
  //printf("Hello DivideCell 1, child->genomesize=%d\n", child->genomesize);
  
  if(MUT_BACKGROUND) BackgroundMutations(child); //This is background mutations, 
			      // break points duplicate and delete here
			      // along with all genes, with probability pMut
  if(MUT_MUTATION) MutMutations(child);	  // mutation of the mutational operators
			  // put here so they'll have effect next generation			     
  //GenomeStats(child);
  CellHealth(child);
}



// Background mutations involve single genes, 
// they can be duplications, deletions and inactivations.
// break points duplicate and delete here
void BackgroundMutations(CELL *ic)
{
  int i,j,k,nmut,mutpos;
    
  //fprintf(stderr,"before mut genome=%s\n",ic->genome);
  //with some small chance, B sites can be randomly inserted in the genome
  // we do this here because we assume that break points appear randomly in the genome
  
//   nmut=(int)bnldev( pow(par_mut,3.) , ic->genomesize);
//   for(i=0;i<nmut;i++){
//     if(genrand_real1() > 0.25) continue; //to balance with the frequency of other duplication elements
// 					 // only in 0.25 cases we will go through
//     mutpos=ic->genomesize*genrand_real2();
//     
//     if(ic->genomesize+1 > MAXGENOMESIZE) {
//       fprintf(stderr,"BackgroundMutations(): Error. genomesize exceeds MAXGENOMESIZE\n");
//       fprintf(stderr,"Attempt to introduce B failed\n");
//       fprintf(stderr,"current genome: %s\n",ic->genome);
//       fprintf(stderr,"The program will exit now.\n");
//       exit(1);
//     }
//       
//     for(j=ic->genomesize; j>mutpos;j--) ic->genome[j]=ic->genome[j-1]; //shift elements
//     
//     ic->genome[mutpos]='B';	//insert B
//     ic->genomesize+=1;
//     ic->genome[ic->genomesize]='\0';
//   }
  
  //All other genes can be duplicated,deleted or inactivated
  // 'B' can be deleted
  
  //In this model, if you hit 'B', another one is created randomly in the genome
  // this makes statistical sense... biologically I'm not sure, but I hope it works
  nmut=(int)bnldev( pow(par_mut,2.) , ic->genomesize);
  for(i=0;i<nmut;i++){
    mutpos=ic->genomesize*genrand_real2(); 
    SingleGeneMut(ic, mutpos,BGmut_scheme);
  }
  
  GenomeStats(ic);
  TargetVolume(ic); // cell vol does not change, but tarvol may
  
}

// this makes a single-gene mutation.
// It is assumed that 
// 1) that the mutation happens was already decided elsewhere
// 2) where it happens it is already decided elsewhere
//  in particular, this function is called by BackgroundMutations() and by TranscrMutations()
//  clearly, it is not possible to duplicate or delete 'B' if the function is called by the latter
int SingleGeneMut(CELL *ic,int mutpos, double *mut_scheme)
{
    int j,k,mut_type=10;
    char genes[N_GENE_TYPES]={'T','Q','R','P'}; // types of genes
    char delet_genes[N_GENE_TYPES]={'t','q','r','p'};
    double whateffect=genrand_real1();
    //fprintf(stderr,"nmut i=%d, mutation mutpos=%d which is %c, muteffect=%f\n",i,mutpos,ic->genome[mutpos],whateffect);
    if( whateffect < mut_scheme[DUPLIC] ){
      mut_type=1;
	if(ic->genomesize+1 > MAXGENOMESIZE){
	  fprintf(stderr,"BackgroundMutations(): Error. genomesize exceeds MAXGENOMESIZE\n");
          fprintf(stderr,"Attempt to duplicate a gene failed\n");
          fprintf(stderr,"current genome: %s\n",ic->genome);
          fprintf(stderr,"The program will exit now.\n");
	  exit(1);
	}
	if(ic->genome[mutpos] == 'B'){
	  int insBpos=ic->genomesize*genrand_real2();	//insert B randomly in the genome
	  for(j=ic->genomesize; j>insBpos;j--) ic->genome[j]=ic->genome[j-1]; //shift elements
          ic->genome[insBpos]='B';	//insert B
	}
	// else if it is any other kind of gene
        else{
          for(j=ic->genomesize; j>mutpos;j--) ic->genome[j]=ic->genome[j-1]; //duplication
	}
        // in both cases genome size increases
        ic->genomesize+=1;
	ic->genome[ic->genomesize]='\0';
	//fprintf(stderr,"hello duplication, genome %s\n", ic->genome);  
      }else if( whateffect < mut_scheme[DUPLIC]+ mut_scheme[DELET]){
	//deletion
	mut_type=-1;
	
	for(j=mutpos; j < ic->genomesize;j++) ic->genome[j]=ic->genome[j+1]; 
	ic->genomesize-=1;
	ic->genome[ic->genomesize]='\0';
	
      }else{
	mut_type= 0;
	//inactivating mutation, if the gene is TQR
	for(k=0;k<N_GENE_TYPES;k++)
	  if( ic->genome[mutpos]==genes[k] ){
	    ic->genome[mutpos] = delet_genes[k];
	    break;
	  }
      }
    //fprintf(stderr,"nmut=%d after mut genome=%s\n",i, ic->genome);  
  
  return mut_type;
  
}

// In this model max_ktrmut[] mutates, NOT mut_genes[]
//mutates the propensity for this or that type of mutation
void MutMutations(CELL *child)
{
  int i;
  //int nmut=(int)bnldev(par_mut,3);
  //for(i=0;i<nmut;i++){
  if(genrand_real1()<par_mut){
    int pos=(int)(3.*genrand_real2());
    double rmut=0.1*par_delta*(genrand_real1() - 0.5);
    //if a mutation happens its effect is more pronounced like this
    child->max_ktrmut[pos] += rmut;
    if(child->max_ktrmut[pos] >= 1.) child->max_ktrmut[pos] = 2. - child->max_ktrmut[pos];
    //if(child->max_ktrmut[pos] <  0.) child->max_ktrmut[pos] *= -1.;   LET THEM GO NEGATIVE
  }
  
}


//new version:
// REMOVE ALL CLUTTER ABOUT break points
// different rates for duplications, deletions, inactivations, independent from each other.
// --------
// 1- makes Mutations on child - ON ACTUAL GENOME !!!
// NOPE: 2- sets transcr load to zero  <-- NOT ANY MORE !!!
// 3- re-calculates target vol
void TranscrMutations_diffratesDupDelIn(CELL *child)
{
  int i,j,k,l,new_genomesize,singlegene_mut_flag,genetype;
  //int mut_count;
  double a,a_trload;
  double pmut;
  char genes[N_GENE_TYPES]={'T','Q','R','P'};	//this and next are used as conversion tables
  char delet_genes[N_GENE_TYPES]={'t','q','r','p'};
  //int howmanymuts_pergenetype[N_GENE_TYPES]={0};
  int mutations_array[child->genomesize];//contains the gene types on which we make mutations, 
					  // repeated as many times as we are mutating them, tot_mut marks its limit
  int mutpos_array[2*child->genomesize]; //marks genome positions where mutations are scheduled occur
					 // its twice genome size because duplications can happen that change location of pos
  int tot_mut=0;
  double delta_step=1./par_init_interval_time;
  //printf("Hello TranscrMutations0.1\n");
  
  // Calculate proportions of duplications, deletions and inactivations based on their rates
  double current_max_ktrmut_DUPLIC = ( child->max_ktrmut[DUPLIC] > 0. )? child->max_ktrmut[DUPLIC] : 0. ;   // to ensure that negative values are treated as zero.
  double current_max_ktrmut_DELET =  ( child->max_ktrmut[DELET]  > 0. )? child->max_ktrmut[DELET] : 0.;
  double current_max_ktrmut_INACT =  ( child->max_ktrmut[INACT]  > 0. )? child->max_ktrmut[INACT] : 0.;
  
  double current_max_ktrmut = current_max_ktrmut_DUPLIC + current_max_ktrmut_DELET + current_max_ktrmut_INACT;  //tot max rate transcr. mutations
  if(current_max_ktrmut>0.){ 
    child->mut_genes[DUPLIC] = current_max_ktrmut_DUPLIC / current_max_ktrmut;
    child->mut_genes[DELET] =  current_max_ktrmut_DELET  / current_max_ktrmut;
    child->mut_genes[INACT] =  current_max_ktrmut_INACT  / current_max_ktrmut;
  }
  
  // Here we see how many mutations occur per gene type, which depends on how much it was transcribed
  // when the transcr load is minimal, mutation rates converge to par_mut as follows
  // mut=par_mut+(1-par_mut)*(  trload^2 / (par_ktrmut^2 + trload^2)  )  
  
  //----- IMPORTANT -----
  // We calculate how many mutations occur per gene type first,
  //  then in which (random) order they occur for all gene types,
  //  and only *after* we find where they occur, 
  //  and we keep track of the genome with some flags that say where mutations can(not) occur
  
  double av_tr_mut=0.;	//average transcr. mut- used for plotting
  
  for(i=0;i<N_GENE_TYPES;i++){
    int tot_act_inact_genes = child->genesum[i]+child->delet_genesum[i];
    pmut=0.;
    
    if(child->trload[i] > 0){
      a_trload = child->trload[i]/((delta_step+child->lifespan) * tot_act_inact_genes); //average transcr load per gene divided by the number of ALL genes over
      
      // THIS IS WHAT WE USED THROUGHOUT 
      //a = (par_ktrmut*par_max_expression_per_gene)/a_trload; //0.9*0.2/<tr_load>
      //a = pow(a,10.);
      //pmut = current_max_ktrmut* ( 1./(1.+a)  );	//this is the per gene mutation rate
      
      //REVIEWER ASKED FOR LINEAR FUNCTION:
      // function should have same extremal behaviour:
      // should be 0 in 0 and pmut=1/(1+ (0.18/x)**10) -> ~0.7415 when x is 0.2
      // this is the line y=(0.74/0.2)*x, thus
      pmut=current_max_ktrmut*(0.7415/0.2)*a_trload; //current_max_ktrmut = par_max_ktrmut if you don't mutate it
      // Similarly the other functions listed in the paper can be intergrated here.
      
    }
    av_tr_mut += pmut/((double)N_GENE_TYPES);
    
    int nmut=(int)bnldev(pmut, tot_act_inact_genes);	// how many mutations on genes of this type
    
    //printf("Hello TranscrMutations0.1, i=%d, trload=%f, pmut=%f,nmut=%d\n",i, child->trload[i],pmut,nmut);
    
    for(j=tot_mut; j<tot_mut+nmut; j++) mutations_array[j]=i;
    tot_mut+=nmut;	//update tot mut count
    mutations_array[tot_mut]=-1; //just to be sure it ends somewhere
  }
  
  child->av_tr_mut = av_tr_mut;
  
  //Now we have the mutations that happen per gene type,
  // we generate the order at which they occur
  // by shuffling mutations_array[]
  
  // the loop goes 
  //   from tot_mut (int random pos=[0, tot_mut-1])
  //   to i=2, because when i=1, int random pos=[0,0] which would be pointless
  for(i=tot_mut; i>1; i--){
    int pos=(int)(i*genrand_real2()); // generate random position between 0 and i-1
    int tmp=mutations_array[pos]; // swap idcell[i] with idcell[random_pos]
    mutations_array[pos]=mutations_array[i-1]; //i-1, or it segfaults
    mutations_array[i-1]=tmp;
  }
  
  //now that we have which and how many mutations occur, 
  // we go on to find where they occur
  
  //first initialise mutpos_array: 
  // mutpos_array is initialised to -1 
  for(i=0;i<child->genomesize;i++) mutpos_array[i]=-1;
  mutpos_array[ child->genomesize ]=-1; //end of genome
  
  for(i=0;i<tot_mut;i++){
    //the mutation will involve a gene of type genes[ mutations_array[i] ]
    // we find which of all the copies on the genome will actually be mutated
    genetype=mutations_array[i];
    double howmany_act_inact=child->genesum[ genetype ] + child->delet_genesum[ genetype ]; // sum active and inactive of type genetype
    int whichgene= ((double)(howmany_act_inact)*genrand_real2());	//this is the gene where mutation occurs
    int countgenes=0;
    
    //locate it on genome
    for(j=0;j<=child->genomesize;j++)
      if(child->genome[j]== genes[ mutations_array[i]] || child->genome[j]== delet_genes[ mutations_array[i] ] ){
	if(countgenes==whichgene) break;
	countgenes++;
      }
    //now j cointains the position.. which we mark by numbering the order in which mutations will happen there ---let's try
    mutpos_array[j]=i+1; //the numbering starts from 1 because all non mutating positions are 0 and end of genome is -1
    if(j==child->genomesize){
      fprintf(stderr,"TranscrMutations(), Error. Gene that should be there was not found.\n");
      exit(1);
    }
  }
  
  //Now we have everything ready for applying mutations:
  // 1) total number is tot_mut
  // 2) the gene type on which mutation happens is in mutations_array[i];
  // 3) we loop over mutpos_array[], which marks the position where the mutation happens
  //    i.e. where mutpos_array[ ] == i+1.
  //    notice that mutpos_array also says the order in which mutations happen
  // when we apply the mutation we have to keep track of things, e.g. mutpos_array needs to be updated in multiple ways:
  // - duplication and deletions shift all values of mutpos_array[ ]
  // - if a duplication or a deletion influences the position of another mutation that mutation should not occur (for simplicity)
  //... hope it comes out ok
  
  //so, iterate over mutations
  for(i=0; i<tot_mut; i++){
    genetype=mutations_array[i]; // assign gene type-number, i.e. the gene type that will mutate
    
    //we find the gene on the genome, for which mutpos_array has index = i+1 (because 0 is where mutations do not happen)
    for(j=0; j<child->genomesize; j++){
      if( mutpos_array[j] == i+1 ){
	//printf("j=%d, mutpos_array[j]=%d\n",j,mutpos_array[j]);
        // j is the position on which a mutation happens
        int mutpos=j;
        
	// we go look for break points -- look left
        int leftB=-1;	//coordinate of B to the left
        int rightB=-1;	//coordinate of B to the right
        for(k=mutpos;k>=0;k--){
	  if(child->genome[k]=='B'){
	    leftB=k;
	    break;
	  }
        }
        // look right
        for(k=mutpos; k < child->genomesize; k++){
	  if(child->genome[k]=='B'){
	    rightB=k;
	    break;
	  }
        }
        
        //printf("genetype=%d, which is %c, mutpos=%d, leftB=%d, rightB=%d\n",genetype,genes[genetype] ,mutpos, leftB, rightB);
        
	//if there are break points left and right of mutpos, special mutations can happen
        //first we see if they are engaged, this depends on some "foldedness" which may also evolve... later?
        if(leftB!=-1 && rightB!=-1 && j!=child->genomesize && j!=child->genomesize){
	  double ldistance = mutpos - leftB;
	  double rdistance = rightB - mutpos; //makes it positive
	  //printf("%f %f\n",ldistance,rdistance);
	  if(ldistance<0 || rdistance<0) fprintf(stderr,"NOOOOOOO\n");
	  double pbreak=exp(-par_DNAbeta*(ldistance*rdistance));
	  if(genrand_real1() < pbreak){
	    // The break points are engaged
	    singlegene_mut_flag=0;
	    
	    //printf("Hello TranscrMutations0.5, pbreak=%f, breakpoint engaged\n",pbreak);
	    //exit(1);
	    //As soon as you decide what mutation you'll make, 
	    // update the part of mutpos inside the two break points.
	    //Later, if dupl or del happened, update mutpos_array as you update the genome
	    
	    double whateffect=genrand_real1();
	    
	    //DUPLICATION OF THE STRETCH
	    if(whateffect < child->mut_genes[0]){
	      //printf("Hello TranscrMutations0.6, duplication for mutation=%d\n",i);
	      //exit(1);
	      //first, update mutpos between break points
	      for(k=leftB;k<rightB;k++) mutpos_array[k]=0;
	      
	      //shift the right part of the genome right of an amount rightB-leftB
	      new_genomesize=child->genomesize+(rightB-leftB);
	      
	      if(new_genomesize > MAXGENOMESIZE){
		fprintf(stderr, "TranscrMutations(), Warning: reached genome limit. Reallocate memory for genome\n"); //maybe allocate more memory?
		fprintf(stderr,"Attempt to duplicate stretch failed\n");
		fprintf(stderr,"current genome: %s\n",child->genome);
		fprintf(stderr,"The program will exit now.\n");
		//ReallocateMemoryForGenome();
		exit(1);
	      }
	      
	      for(l=child->genomesize-1, k=0; l > rightB; l--, k++){ 
		child->genome[ new_genomesize -1 -k ] = child->genome[l]; //shift the part of genome right of rightB
		mutpos_array[ new_genomesize -1 -k ] = mutpos_array[l]; // shift mutpos_array
	      }
	      for(l=leftB+1,k=0; l<=rightB; l++,k++) {
		child->genome[rightB+1+k] = child->genome[l]; // copy the duplicated genome
		mutpos_array[rightB+1+k] = mutpos_array[l];   // copy mutpos_array
	      }
	      
	      child->genome[new_genomesize]='\0';
	      mutpos_array[new_genomesize]=-1;
	      child->genomesize=new_genomesize;
	      //printf("Duplication occurred, current genome: %s\n",child->genome);
	      //printf("Hello TranscrMutations0.51, new genomesize=%d new mutpos_array (line below is genome) is:\n",child->genomesize);
              //for(i=0;i<child->genomesize;i++) printf("%d ",mutpos_array[i]);
              //printf("\n");
	      
	    }else if(whateffect< child->mut_genes[0]+child->mut_genes[1]){
	      //DELETION
	      //printf("Hello TranscrMutations0.6, deletion for mutation %d\n",i);
	      //exit(1);
	      new_genomesize=child->genomesize - (rightB-leftB);
	      for(l=rightB,k=0;l<child->genomesize;l++,k++){
		child->genome[leftB+k]=child->genome[l];
		mutpos_array[leftB+k]=mutpos_array[l];
	      }
	      child->genome[new_genomesize]='\0';
	      mutpos_array[new_genomesize]=-1;
	      child->genomesize=new_genomesize;
	      //printf("Deletion occurred, current genome: %s\n",child->genome);
	      //exit(1);
	    }else{
	      //inactivating mutation despite involving break point
	      //it can be that despite break points exist and are close enough to work,  
	      // the mutation that actually occurs is simply inactivating
	      //(if the gene involved is already in deleterious count -> no effect)
	      
	      singlegene_mut_flag=1;
	    }
	  }else{
	    //singlegene_mut_flag=1;
	    
	    singlegene_mut_flag=2; // mutation because break point was not engaged
	    
	    //printf("Break point not engaged, singlegene_mut_flag activated");
	  }  
	    
	}// left or right breakpoint was not found
	else{
	  //single gene mutations (mostly inactivating ? )
	  
	  //singlegene_mut_flag=1;
	    
	  singlegene_mut_flag=2;	// mutation because no break point was found (or pair of B points)
        }
	
	//Because 
	//  1) breakpoint is engaged but no dupl/del happened, or
	//  2) break point exists but was not engaged, or
	//  3) break point pair does not exists
	// a single-gene mutation occurs:
        // which mutation is chosen depends on the value of singlegene_mut_flag
        // singlegene_mut_flag==1: inactivation (just as a possible (evolvable) outcome of break point activation)
        // singlegene_mut_flag==2: mutation with same scheme as background mutations
      if(singlegene_mut_flag==1){
	if(child->genome[mutpos] == genes[genetype]) child->genome[mutpos] = delet_genes[genetype]; //mutpos is where the gene is, genetype is the type number
	mutpos_array[mutpos]=0; 	// this is for safety, in practice it should not happen that we mutate this position again
	singlegene_mut_flag=0;
	//printf("mutflag activated, inact mut, current genome: %s\n",child->genome);
      }else if(singlegene_mut_flag==2){
	// Single gene mutation
	if(child->genome[mutpos]=='B'){
	  fprintf(stderr,"TranscrMutations(): Error. Got a 'B' character to mutate.\n");
	  exit(1);
	}
	int which_mut=SingleGeneMut(child,mutpos,child->mut_genes);	// In this model mut_genes[] is just the proportions of max_ktrmut[]	// if write child->mut_genes instead of TRmut_scheme mutations evolve!
	// here we have to update also other arrays, depending on the mutation
	switch(which_mut){
	  
	  case 1: 	// a duplication happened, set mutpos_array and shift mutpos_array to the right
	    mutpos_array[mutpos]=0;
	    for(k=child->genomesize; k>mutpos;k--) mutpos_array[k]=mutpos_array[k-1]; //duplication
	    mutpos_array[child->genomesize]=-1;
	    break;
	    
	  case -1:	// a deletion happened, shift mutpos_array to the left
	    for(k=mutpos; k < child->genomesize;k++) mutpos_array[k]=mutpos_array[k+1]; 
	    mutpos_array[child->genomesize]=-1;
	    break;
	    
	  case 0:	// inactivating mutations, 
	    break;
	  default:
	    fprintf(stderr,"TranscrMutations(): Error. SingleGeneMut() returned %d\n",which_mut);
	    exit(1);
	}
	
	singlegene_mut_flag=0;
	//printf("mutflag activated, inact mut, current genome: %s\n",child->genome);
      }
      //update genome stats and recalculate mutations (if one mutation deletes many genes,
      // the few remaining should not mutate as much as if only few mutations were made)
      GenomeStats(child);
	
	
      } //if statement was: if( mutpos_array[j] == i+1 )
    }//end of for loop to look for the mutation number i
     // we did not find the position on the genome, 
     // probably because it was deleted, or it hitch-hiked some other duplication
     
  }  // END OF MUTATIONAL DYNAMICS FOR mutation number 'i'
    
  //printf("Hello TranscrMutations2\n" );
  TargetVolume(child);	//recomputes target volume for child
}


//here we mutate regulation parameters, i.e.
// ktrscrE[4],ktrscrQ[4],ktrscrRp[4],ktrscrRr[4],ktrscrT[4]
void MutationsRegulation(CELL *ic)
{
  int i,gene_mut,prom_mut;
  int nmut=(int)bnldev(par_mut, N_GENE_TYPES*2);	// 5 genes, 2 signals each
  //double howmuch;
  
  if(nmut==0) return;
  
  for(i=0;i<nmut;i++)
  {
    gene_mut=(int)(N_GENE_TYPES*genrand_real2());
    prom_mut=(int)(2.*genrand_real2());
    
    switch(gene_mut){
      case 0:
	RegMutate( &(ic->ktrscrT[prom_mut]) );
	break;
      case 1:
	RegMutate( &(ic->ktrscrQ[prom_mut]) );
	break;
      case 2:
	RegMutate( &(ic->ktrscrRr[prom_mut]) );
	break;
       case 3:
 	RegMutate( &(ic->ktrscrRp[prom_mut]) );
 	break;
//       case 4:
// 	RegMutate( &(ic->ktrscrE[prom_mut]) );
// 	break;
      default:
	fprintf(stderr,"MutationsRegulation(). Error, something went wrong\n");
	exit(1);
    }
  }
}

void RegMutate(double *regpar)
{
  if(  (*regpar > -1.) || (*regpar < 1.)  ) *regpar += par_max_expression_per_gene*par_delta*(genrand_real1() - 0.5); 
  else *regpar += fabs(*regpar) * par_delta*(genrand_real1() - 0.5); 
  
  if(*regpar> MAX_PAR_EXPRESSION) *regpar= 2.*MAX_PAR_EXPRESSION - *regpar;
  if(*regpar<-MAX_PAR_EXPRESSION) *regpar= -2.*MAX_PAR_EXPRESSION - *regpar ;
}

//here we mutate the metabolic parameters
// i.e. the enzyme reaction constants:
// ks1,ks2,ks3,ksm1,ksm3,kc1,kc2,kc3,kcm1,kcm3;
// there is ten of them, also, we'll be using this order.
void MutationsMetabolism(CELL *ic)
{
  int i, who_mut;
  
  //we first calculate how many of these mutate,
  //I should be using a hypergeometric distribution, 
  //by the binomial distr (i.e. with replacement) included in the mersenne twister will do 
  // on a second thought, if clusters of genes are large, 
  // then drawing mutations with replacement actually makes sense... maybe?
  
  int nmut=(int)bnldev(par_mut, N_METAB_PAR);
  
  if(nmut==0) return;
  
  //here we decide who mutates 
  for(i=0;i<nmut;i++){
    
    who_mut= (int)( N_METAB_PAR*genrand_real2());
    switch(who_mut){
      case 0: 
	MetMutate(&(ic->kt));
	break;
      case 1:
	MetMutate(&(ic->alfat));
	break;
      case 2:
	MetMutate(&(ic->betat));
	break;
//       case 3:
// 	MetMutate(&(ic->ke));
// 	break;
//       case 4:
// 	MetMutate(&(ic->alfae));
// 	break;
//       case 5:
// 	MetMutate(&(ic->betae));
// 	break;
      default:
	fprintf(stderr,"MutationsMetabolism(). Error, something went wrong\n");
	exit(1);
    }
  }
  //fprintf(stderr,"Oh joy, Evolution!\n");
}


void MetMutate(double *metpar)
{
  if(*metpar<1.) *metpar += par_delta*(genrand_real1() - 0.5); //so that if metpar is small you still mutate
  else *metpar += (*metpar)*par_delta*(genrand_real1() - 0.5);
  //if(*metpar < 0.) *metpar *= -1;
  if(*metpar < par_min_met_par) *metpar = 2.*par_min_met_par - *metpar;
  if(*metpar > par_max_met_par) *metpar = 2.*par_max_met_par - *metpar;
}

// This updates the health of a cell health:[0,1] based on some constraints
// which in turn can be used as penalty for transcription regulation in the ODEs system (see main.c)
// 1) first and foremost, if you do not have any active gene of a certain type your health is 0.
// 2) your expressionof homeostatic protein Q has to be at the target
void CellHealth(CELL *ic)
{
  int i;
  
  if(ic->cellvol<=0.){
    ic->health=0.;
    return;
  }
  
  for(i=0;i<N_GENE_TYPES;i++)
    if(ic->genesum[i]<=0){
      ic->health=0.;
      return;
    }
  
  double Qdistance=((*(ic->Q))/ic->cellvol - par_tar_Q);
  // health can be extended to more stuff by summing other terms to Qdistance (e.g. one day toxicity)
  ic->health = exp( - transcr_penalty * Qdistance*Qdistance );
  
}

// CELL VOLUME: proteins + rRNA
// other mRNAs and metabolites do not count for volume
// ribosomes count as 2 molecules
void CellVolume(CELL *icel)
{
  double volume;
  volume =  *(icel->T) + *(icel->Q) + *(icel->Rr) + *(icel->Rp);
  //volume += *(icel->Crmt)*2. + *(icel->Crmq)*2.; // times 2 because there are two molecules per complex
  
  icel->cellvol=volume;
}

// CELL TARGET VOLUME, for now this is just dependent on genome size
// tarvol=par_k_tarvol*(n_genes)^0.9
void TargetVolume(CELL *ic)
{
  ic->tarvol=par_k_tarvol*pow((double)ic->genomesize,par_genome_to_volume_scale);
}

//takes genome and updates genesum, delet_genesum and genomesize
void GenomeStats(CELL *ic)
{
  int i;
  //zero the gene counts in genesum and delet_genesum
  for(i=0;i<N_GENE_TYPES;i++){
    ic->genesum[i]=0;
    ic->delet_genesum[i]=0;
  }
  //fprintf(stderr,"begin genome\n");
  //for(i=0;i<ic->genomesize+1;i++) fprintf(stderr,"i=%d -%c-\n",i,ic->genome[i]);
  //fprintf(stderr,"end genome\n");
  for(i=0;i<MAXGENOMESIZE;i++){
    //fprintf(stderr,"hello0.1\n");
    if(ic->genome[i]=='\0') break;	//end of genome
    
    //fprintf(stderr,"hello0.2\n");
    switch (ic->genome[i]){
      case 'T': ic->genesum[GENE_T]++;
                break;
      case 't': ic->delet_genesum[GENE_T]++;
                break;
      case 'Q': ic->genesum[GENE_Q]++;
                break;
      case 'q': ic->delet_genesum[GENE_Q]++;
                break;
      case 'R': ic->genesum[GENE_Rr]++;
                break;
      case 'r': ic->delet_genesum[GENE_Rr]++;
                break;
      case 'P': ic->genesum[GENE_Rp]++;
                break;
      case 'p': ic->delet_genesum[GENE_Rp]++;
                break;
      case 'B': break;
      default:
	fprintf(stderr,"GenomeStats(), Error. Invalid character in genome\n");
	exit(1);
    }
    //fprintf(stderr,"hello0.3\n");
  }
  //fprintf(stderr,"old genomesize=%d ",ic->genomesize);
  ic->genomesize=i;
  //fprintf(stderr,"new genomesize=%d\n",ic->genomesize);
  //fprintf(stderr,"hello2\n");
}

