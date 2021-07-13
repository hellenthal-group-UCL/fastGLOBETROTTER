## This software is licenced under the "Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International Public License" (http://creativecommons.org/licenses/by/4.0/)
## TAKE CHROMOPAINTER OUTPUT AND, FOR A SPECIFIED "RECIPIENT" POPULATION:
##      (1) INFERS DATES, PROPORTIONS AND SOURCES OF ADMIXTURE EVENT(S) IN ANCESTRAL HISTORY OF "RECIPIENT" POPULATION (IN PARTICULAR INFERRING THE ADMIXING SOURCES AS A LINEAR COMBINATION OF OTHER "DONOR" POPULATIONS)
##      (2) USING BOOTSTRAP RE-SAMPLING, INFERS CONFIDENCE INTERVAL AROUND ADMIXTURE DATES

## !!!! NOTE: MUST HAVE PACKAGE "nnls" INSTALLED !!!!


## usage:   (1) R CMD SHLIB -o fastGLOBETROTTERCompanion.so fastGLOBETROTTERCompanion.c -lz    (do only once)
## then:    (2) R < fastGLOBETROTTER.R parameter.infile.list samplefiles.infile.list recomratefiles.infile.list option --no-save > output.out
 
## example -- SOURCE/DATE/PROPORTION ESTIMATION:  R < fastGLOBETROTTER.R tutorial/paramfile.txt tutorial/samplefile.txt tutorial/recomfile.txt 1 --no-save > output.out

########################################
## COMMAND LINE INPUT:
ptm  = proc.time()
usage=function()
{
	print(noquote("run using: R < GLOBETROTTER.R parameter_infile painting_samples_filelist_infile recom_rate_filelist_infile --no-save > screen_output"))
	print(noquote("parameter input file format (NO defaults, all fields must be entered, even if not used):"))
	print(noquote("prop.ind: [0,1]"))
	print(noquote("bootstrap.date.ind: [0,1,2]"))
	print(noquote("null.ind: [0,1]"))
	#print(noquote("contrasts.ind: [0,1]"))
	print(noquote("input.file.ids: [input.filename1]"))
	print(noquote("input.file.copyvectors: [input.filename2]"))
	print(noquote("save.file.main: [output.filename1]"))
	print(noquote("save.file.bootstraps: [output.filename2]"))
	print(noquote("copyvector.popnames: [pop_1 pop_2 ... pop_j]"))
	print(noquote("surrogate.popnames: [pop_1 pop_2 ... pop_k]"))
	print(noquote("target.popname: [pop_rec]"))
	print(noquote("num.mixing.iterations: [0,1,...,5,...]"))
	print(noquote("props.cutoff: [0.0,...,1.0]"))
	print(noquote("bootstrap.num: [0,1,...]"))
	print(noquote("num.admixdates.bootstrap: [1,2]"))
	print(noquote("num.surrogatepops.perplot: [1,...]"))
	print(noquote("curve.range: [lower.lim upper.lim (in cM)]"))
	print(noquote("bin.width: [e.g. 0.1]"))
	print(noquote("xlim.plot: [lower.lim upper.lim (in cM)]"))
	print(noquote("prop.continue.ind: [0,1]"))
	print(noquote("haploid.ind: [0,1]"))
}

temp=commandArgs()

param.infile=as.character(temp[2])
if (param.infile=="help"){usage();q(save='no')}

######################################
## COMMAND LINE CHECK AND READ IN:

weightmin=0     ## a chunk must have > weightmin SNPs to be considered in the coancestry curves
gridweightMAXcutoff=1   ## any chunk's weighted contribution (i.e. the cM length of the chunk) to the coancestry curves cannot be bigger than this (in cM)
null.calc=0

error.input.message=function(file.name)
  {
    print(paste("Something wrong with input file ",file.name,". See below. Exiting...",sep=''))
    usage()
    q(save='no')
  }
line.check=function(file.name,skip.val,match.val)
  {
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) != as.character(match.val))
      {
        error.input.message(file.name)
      }
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) == as.character(match.val))
      return(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[2:length(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE))])
  }

          ## line read in and checks:
prop.ind=as.integer(line.check(param.infile,0,"prop.ind:"))
if ((length(prop.ind)!=1)||(prop.ind!=0 && prop.ind!=1)) error.input.message(param.infile)
admix.curves.ind=as.integer(line.check(param.infile,1,"bootstrap.date.ind:"))
if ((length(admix.curves.ind)!=1)||(admix.curves.ind!=0 && admix.curves.ind!=1 && admix.curves.ind!=2 )) error.input.message(param.infile)
null.ind=as.integer(line.check(param.infile,2,"null.ind:"))
if ((length(null.ind)!=1)||(null.ind!=0 && null.ind!=1)) error.input.message(param.infile)
#contrasts.ind=as.integer(line.check(param.infile,3,"contrasts.ind:"))
#if ((length(contrasts.ind)!=1)||(contrasts.ind!=0 && contrasts.ind!=1)) error.input.message(param.infile)
id.file=as.character(line.check(param.infile,3,"input.file.ids:"))
if (length(id.file)!=1) error.input.message(param.infile)
copyvector.file=as.character(line.check(param.infile,4,"input.file.copyvectors:"))
if (length(copyvector.file)!=1) error.input.message(param.infile)
save.file.props=as.character(line.check(param.infile,5,"save.file.main:"))
if (length(save.file.props)!=1) error.input.message(param.infile)
save.file.bootstraps=as.character(line.check(param.infile,6,"save.file.bootstraps:"))
if (length(save.file.props)!=1) error.input.message(param.infile)
donor.pops.all=unique(as.character(line.check(param.infile,7,"copyvector.popnames:")))
if (length(donor.pops.all)<2) error.input.message(param.infile)
newnames=as.character(line.check(param.infile,8,"surrogate.popnames:"))
if (length(unique(newnames))<2) error.input.message(param.infile)
recipient.pops=as.character(line.check(param.infile,9,"target.popname:"))
if (length(recipient.pops)!=1) error.input.message(param.infile)
popvec.admixcurves=recipient.pops
num.iterations = as.integer(line.check(param.infile,10,"num.mixing.iterations:"))
if ((length(num.iterations)>1) || (num.iterations < 0)) error.input.message(param.infile)
prop.cutoff = as.double(line.check(param.infile,11,"props.cutoff:"))
if ((length(prop.cutoff)>1) || (prop.cutoff < 0) || (prop.cutoff>=1)) error.input.message(param.infile)
bootstrap.samp = as.integer(line.check(param.infile,12,"bootstrap.num:"))
if ((length(bootstrap.samp)!=1) || (bootstrap.samp < 0)) error.input.message(param.infile)
nparam.toplot = as.integer(line.check(param.infile,13,"num.admixdates.bootstrap:"))      ## plot results from fitting one or two events?
if ((length(bootstrap.samp)!=1) || ((nparam.toplot!=1) && (nparam.toplot!=2))) error.input.message(param.infile)
numpops.toplot = as.integer(line.check(param.infile,14,"num.surrogatepops.perplot:"))    ## number of donor pops to include in plots
grid.range=as.integer(line.check(param.infile,15,"curve.range:"))  ## integer values for where to start and end fitting on the admixture copying grid
if ((length(grid.range)!=2)||(grid.range[1]<0)||(grid.range[2]<0)||(grid.range[1]>=grid.range[2])) error.input.message(param.infile)
bin.width=as.double(line.check(param.infile,16,"bin.width:"))  ## width of bins for tabulating chunk-pair counts when generating coancestry curves
if ((length(bin.width)!=1)||(bin.width[1]<=0)) error.input.message(param.infile)
xlim.plot=as.double(line.check(param.infile,17,"xlim.plot:"))   ## range of x-axis for admixture plots
if (((length(xlim.plot)==1) && (xlim.plot[1]>0)) || (xlim.plot[1]>=xlim.plot[2]) || (xlim.plot[2]<=0) || (length(xlim.plot)>2)) error.input.message(param.infile)
prop.continue.ind=as.integer(line.check(param.infile,18,"prop.continue.ind:"))
if ((length(prop.continue.ind)!=1)||(prop.continue.ind!=0 && prop.continue.ind!=1)) error.input.message(param.infile)
haploid.ind=as.integer(line.check(param.infile,19,"haploid.ind:"))
if ((length(haploid.ind)!=1)||(haploid.ind!=0 && haploid.ind!=1)) error.input.message(param.infile)
ploidy=2-haploid.ind
prop.save.post=".txt"
if (prop.continue.ind==0) prop.save.addon = ""
if (prop.continue.ind==1) prop.save.addon = "_continue"

	if (prop.continue.ind==1 && prop.ind==1)
	{
	    print(paste("YOU HAVE SPECIFIED BOTH 'prop.ind' AND 'prop.continue.ind' TO BE EQUAL TO 1 IN ",param.infile,". IF THERE IS AN EXISTING 'save.file.main', WHICH YOU WISH TO USE AS STARTING VALUES FOR YOUR ESTIMATION, THEN YOU WANT 'prop.continue.ind' TO BE 1. IF THERE IS NO EXISTING 'save.file.main', THEN YOU WANT 'prop.ind' TO BE 1. Exiting...",sep=''))
	    q(save='no')
	}
	if (max(c(prop.ind,prop.continue.ind,admix.curves.ind))==0)
	{
	    print(paste("YOU HAVE SPECIFIED NO ACTION IN ",param.infile,", WITH 'prop.ind', 'bootstrap.date.ind' and 'prop.continue.ind' ALL EQUAL TO 0. Exiting...",sep=''))
	    q(save='no')
	}
	if (admix.curves.ind==1 && bootstrap.samp==0)
	{
	    print(paste("YOU HAVE SPECIFIED TO BOOTSTRAP DATES IN ",param.infile," BUT HAVE REQUESTED 0 BOOTSTRAPS. Exiting...",sep=''))
	    q(save='no')
	}

#######################################
## PROGRAM:

dists=seq(grid.range[1],grid.range[2],bin.width)
num.bins=length(dists)
if (length(dists)<=1) {print(paste("Something wrong with bin.width and curve.range specification in ",param.infile,"! Exiting....",sep='')); q(save='no')}
#num.bins=ceiling((grid.range[2]-grid.range[1])/bin.width)

#############################################################
##############################################################
## (0) FUNCTIONS:
##############################################################
#############################################################

#install.packages("nnls")   ## IF THE FOLLOWING LINE DOESN'T WORK, USE THIS!
library(nnls)

options(scipen=999)
dyn.load("fastGLOBETROTTERCompanion.so")


findpeak  = function  (ourmix,means,bin.width) 
{
  
  dom = which(ourmix ==max(ourmix))
  new = means[length(ourmix)*(dom-1)+dom,]
 ##if (dim(means)[2] > 100) {wsize = 100  ## !!!!! pongsakorn orig !!!!!!!!
  ##} else  {wsize = round(dim(means)[2]/4)}  ## !!!!! pongsakorn orig !!!!!!!!
  wsize = round(dim(means)[2]/2) ## !!!!! garrett change (now can remove up to half the original -- fitted -- x-axis) !!!!!!!!
  lrange = c(3,5,7,9,11,13)
  mid = c(2,3,4,5,6,7)  
  keep_max = NULL
  c = 1
  for(i in lrange){
    slope = NULL
    j = 1
    ##while (j < wsize)  ## !!!!! pongsakorn orig !!!!!!!!
    while (j < wsize && ((j+i)<(dim(means)[2]-10)))  ## !!!!! garrett change (now must keep the last 10 gridpoints at least) !!!!!!!!
    {
      temp = lm(new[j:(j+i)]~c(j:(j+i)))
      temp = temp$coefficients[2]
      slope = append(slope,temp)
      j = j + 1
    }
    index.val=which(slope<0)  ## !!!!! garrett change !!!!!!!!
    if (sum(index.val)==0) keep_max = append(keep_max, length(slope))  ## !!!!! garrett change !!!!!!!!
    if (sum(index.val)>0) keep_max = append(keep_max, min(index.val))  ## !!!!! garrett change !!!!!!!!
    ##keep_max = append(keep_max, min(which(slope<0)))   ## !!!!! pongsakorn orig !!!!!!!!
    c = c+1
  }
  #print(new)
  #print(keep_max)
  #print(max(keep_max+mid))
  #plot(new[max(keep_max+mid):200])
   #print(keep_max[which(keep_max == max(keep_max))]-1+mid[which(keep_max == max(keep_max))])
  if (mean(keep_max) == 1) {return(0)}
    #else return(1)
     #else return(max(keep_max+mid))
  else return(max(keep_max[which(keep_max == max(keep_max))]-1+mid[which(keep_max == max(keep_max))]))
 
}


coancestry.curves.mode3=function(ploidy,ind.id.vec,infile.samples,infile.recom,bin.width,num.bins,grid.start,unique.pops,pop.label.vecII,rec.labels,weightmin,gridweightMAXcutoff,temp.weights,num.surrogates)
{
	means=matrix(0,nrow=num.surrogates^2,ncol=num.bins)
	out=.C("data_read_mode3",ploidy=as.integer(ploidy),ind_id_vec=as.integer(ind.id.vec),filenameSAMPread=as.character(infile.samples),filenameRECOMread=as.character(infile.recom),binwidth=as.double(bin.width),num_bins=as.integer(num.bins),grid_start=as.integer(grid.start),num_donors=as.integer(length(unique.pops)),donor_label_vec=as.integer(pop.label.vecII),num_inds_rec=as.integer(length(rec.labels)),recipient_label_vec=as.character(rec.labels),weight_min=as.double(weightmin),gridweightMAX_cutoff=as.double(gridweightMAXcutoff),num_surrogates=as.integer(num.surrogates),temp_weights=as.double(temp.weights),means=as.double(means))
	return(out$means)
}


coancestry.curves=function(ploidy,ind.id.vec,infile.samples,infile.recom,bin.width,num.bins,grid.start,unique.pops,pop.label.vecII,rec.labels,weightmin,gridweightMAXcutoff,temp.weights,num.surrogates)
{
	means=matrix(0,nrow=num.surrogates^2,ncol=num.bins)
	out=.C("data_read",ploidy=as.integer(ploidy),ind_id_vec=as.integer(ind.id.vec),filenameSAMPread=as.character(infile.samples),filenameRECOMread=as.character(infile.recom),binwidth=as.double(bin.width),num_bins=as.integer(num.bins),grid_start=as.integer(grid.start),num_donors=as.integer(length(unique.pops)),donor_label_vec=as.integer(pop.label.vecII),num_inds_rec=as.integer(length(rec.labels)),recipient_label_vec=as.character(rec.labels),weight_min=as.double(weightmin),gridweightMAX_cutoff=as.double(gridweightMAXcutoff),num_surrogates=as.integer(num.surrogates),temp_weights=as.double(temp.weights),means=as.double(means))
	return(out$means)
}

mutiply.temp.weight=function(unique.pops,num.surrogates,weights.mat)
{
	tempweights=matrix(0,nrow=num.surrogates*num.surrogates,ncol=unique.pops*unique.pops)
	out=.C("mutiply_temp_weight",tempweights = as.double(tempweights),weights_mat = as.double(weights.mat),num_surrogates = as.integer(num.surrogates),num_donors = as.integer(unique.pops))
	return(out$tempweights)
}

coancestry.curves.null=function(ploidy,infile.samples,infile.recom,bin.width,num.bins,grid.start,target.pop,unique.pops,pop.label.vecII,rec.labels,weightmin,gridweightMAXcutoff)
{
	chrom.ind.res=matrix(0,nrow=length(unique.pops)^2,ncol=num.bins)
	out=.C("data_read_null",ploidy=as.integer(ploidy),filenameSAMPread=as.character(infile.samples),filenameRECOMread=as.character(infile.recom),binwidth=as.double(bin.width),num_bins=as.integer(num.bins),grid_start=as.integer(grid.start),num_donors=as.integer(length(unique.pops)),donor_label_vec=as.integer(pop.label.vecII),num_inds_rec=as.integer(length(rec.labels)),recipient_label_vec=as.character(rec.labels),weight_min=as.double(weightmin),gridweightMAX_cutoff=as.double(gridweightMAXcutoff),chrom_ind_res=as.double(chrom.ind.res))
	if (is.element(target.pop,unique.pops))
	{
	   target.ind=(1:length(unique.pops))[unique.pops==target.pop]
       to.remove=sort(unique(c(seq(target.ind,length(unique.pops)^2,length(unique.pops)),((target.ind-1)*length(unique.pops)+1):(target.ind*length(unique.pops)))))
	   out$chrom_ind_res=c(t(matrix(out$chrom_ind_res,byrow=TRUE,ncol=num.bins)[-to.remove,]))
	}
	return(out$chrom_ind_res)
}

corr.find=function(x,vec.tocheck)
  {
    return(cor(x,vec.tocheck))
  }

#####function to more efficiently fit an initial mixture, essentially same as stepwise regression results but much faster.
getfit=function(predmat,fitdata,restrict=1){

  temp=matrix(predmat[-restrict,],ncol=dim(predmat)[2])
  for(i in 1:nrow(temp)) temp[i,]=temp[i,]-predmat[restrict,]

  fitdata2=fitdata-predmat[restrict,]

  v=nnls(t(temp),fitdata2)
  x=v$x
  newx=1:nrow(predmat)
  newx[!((1:nrow(predmat))==restrict)]=x
  newx[restrict]=1-sum(x)
  v$x=newx
  names(v$x)=rownames(predmat)
  
  return(v)
}

#####the following does an initial linear model fit quicker than before...
getoverallfit=function(predmat,fitdata){
  restrict=1
  rep=1
  i=1
  while(rep==1){
    q=getfit(predmat,fitdata,restrict=i)
    
    if(q$x[i]>0) rep=0
    i=i+1
  }
  
  return(q)
}

#####need to require sums to equal 1
getfitfull=function(newpredmat,target,restrict){

  temp=newpredmat[,-restrict]
  for(i in 1:(ncol(temp)/2)) temp[,i]=temp[,i]-newpredmat[,restrict[1]]
  for(i in (ncol(temp)/2+1):(ncol(temp))) temp[,i]=temp[,i]-newpredmat[,restrict[2]]
  target2=target-newpredmat[,restrict[1]]-newpredmat[,restrict[2]]

  v=nnls(temp,target2)
  x=v$x
  newx=1:ncol(newpredmat)
  newx[!((1:ncol(newpredmat))%in%restrict)]=x
  newx[restrict[1]]=1-sum(x[1:(length(x)/2)])
  newx[restrict[2]]=1-sum(x[(length(x)/2+1):(length(x))])
  v$x=newx

  return(v)
}

####two populations, need 
getoverallfitfull=function(newpredmat,target){
  ourguess=nnls(newpredmat,target)$x
     ####use as a ranking

  guess=matrix(ourguess,nrow=2,byrow=T)
  orderi=order(guess[1,],decreasing=T)
  orderj=order(guess[2,],decreasing=T)+length(orderi)

  maxi=ncol(newpredmat)/2
  minj=ncol(newpredmat)/2+1
  maxj=ncol(newpredmat)
  
  finish=0;
  for(i in orderi){
    for(j in orderj){
      
      q=getfitfull(newpredmat,target,restrict=c(i,j))
      if(q$x[i]>0 & q$x[j]>0) finish=1;
      if(finish==1) break;
    }
    if(finish==1) break;
  }
  return(q)
}

getafit=function(alpha,mu,weights.mat, predmat, intercepts, fitdata){

  test=weights.mat%*%t(predmat)

  temppredmat=predmat
  vals=colSums(weights.mat)
            ######deal with lack of self copying
  for(i in 1:nrow(predmat)) temppredmat[i,]=temppredmat[i,]*vals
  for(i in 1:nrow(predmat)) temppredmat[i,]=temppredmat[i,]/sum(temppredmat[i,])

  newpredmat=cbind(mu*sqrt(alpha)*sqrt(1-alpha)*test,-mu*sqrt(alpha)*sqrt(1-alpha)*test)
  newpredmat=rbind(newpredmat,cbind((1-mu)*alpha*t(predmat), (1-mu)*(1-alpha)*t(predmat)))

  target=c(mu*intercepts,(1-mu)*fitdata)

  r=getoverallfitfull(newpredmat,target)
  names(r$x)=colnames(newpredmat)
  
  return(r)
}

########function that uses the above functions to obtain information on fitted values and likelihoods
fitwhoadmixed=function(weights.mat,predmat,intercepts,fitdata){
	alpha.vec=seq(0.01,0.99,0.01)
	mu=0.5
	tol=1e-5
	results=vector(length=0)
	for(alpha in alpha.vec){
		  mu=sqrt(length(fitdata)*mean(fitdata^2))/(sqrt(length(fitdata)*mean(fitdata^2))+sqrt(length(intercepts)*mean(intercepts^2)))
		  ####this is close to 1 if fitdata is long compared to intercepts
		  ####so intercepts get upweighted in resids if there are few of them, which seems right
		  ####also small if intercepts are large, in which case accept reasonable error level given large values

		  q=getafit(alpha,mu,weights.mat,predmat,intercepts,fitdata)
		  munew1=mean(q$res[1:length(intercepts)]^2)/mu^2
		  munew2=mean(q$res[(length(intercepts)+1):length(q$res)]^2)/(1-mu)^2
		  lhood=-nrow(weights.mat)*log(munew1)/2-nrow(predmat)*log(munew2)/2

	               ####q is all we need
		   mufitted=mu
		   lhood=lhood
		   errors1=q$res[1:length(intercepts)]
		   errors2=q$res[(length(intercepts)+1):length(q$res)]

                      ###think about model
		   errors1=errors1/mu
		   errors2=errors2/(1-mu)
		   resids=c(errors1,errors2)
		   curverelmse=mean(errors1^2)/mean(intercepts^2)
		   proprelmse=mean(errors2^2)/mean(fitdata^2)
		   avrelmse=(curverelmse+proprelmse)/2
		   fittedintercepts=intercepts-errors1
		   fittedprops=fitdata-errors2
		   curvecor=cor(fittedintercepts,intercepts)
		   propcor=cor(fittedprops,fitdata)
		   avcor=(curvecor+propcor)/2
		   pop1=q$x[1:(length(q$x)/2)]
		   pop2=q$x[(length(q$x)/2+1):length(q$x)]
		   pop1fit=pop1[pop1>0]
		   names(pop1fit)=names(pop1)[pop1>0]
		   pop2fit=pop2[pop2>0]
		   names(pop2fit)=names(pop2)[pop2>0]
		   mixturefit=alpha*pop1+(1-alpha)*pop2
		   mixturefit=mixturefit[mixturefit>0]

                      ####explore correlations
        	   z=q$x[1:(length(q$x)/2)]
		   if(length(pop1fit)){
			that=t(matrix(predmat[names(pop1fit),],ncol=ncol(predmat))) %*%matrix(pop1fit,ncol=1)
			msevec=1:nrow(predmat);
			for(i in 1:nrow(predmat)) msevec[i]=mean((predmat[i,]-that)^2)
			bestpop1mse=min(msevec)
			optimalpop1group=rownames(predmat)[which(msevec==bestpop1mse)]
                   }
	 	   else{ bestpop1mse=1;optimalpop1group=NA;}

	 	   if(length(pop2fit)){
			that=t(matrix(predmat[names(pop2fit),],ncol=ncol(predmat))) %*%matrix(pop2fit,ncol=1)
	 		msevec=1:nrow(predmat);
	 		for(i in 1:nrow(predmat)) msevec[i]=mean((predmat[i,]-that)^2)
	 		bestpop2mse=min(msevec)
	 		optimalpop2group=rownames(predmat)[which(msevec==bestpop2mse)]
	           }
	 	   else{ bestpop2mse=1;optimalpop2group=NA;}

	 	   avbestmse=(bestpop1mse+bestpop2mse)/2
	 	   curresults=list(alpha=alpha,pop1fit=pop1fit,pop2fit=pop2fit,mixturefit=mixturefit,lhood=lhood,curverelmse=curverelmse,proprelmse=proprelmse,avrelmse=avrelmse,curvecor=curvecor,propcor=propcor,avcor=avcor,bestpop1mse=bestpop1mse,optimalpop1group=optimalpop1group,bestpop2mse=bestpop2mse, optimalpop2group=optimalpop2group,avbestmse=avbestmse,resids=resids,mufitted=mufitted)

	 	   results=c(results,curresults)
         }
	 return(results)
}

######new approach to get intercepts - need a matrix "weights.mat" to tell us dimensions, and fitvec to get rid of frequencies
getfitted=function(weights.mat,fitvec,interceptmat,minval=0,pcnum=1){
  d=matrix(1:nrow(weights.mat),nrow=1)
  colnames(d)=rownames(weights.mat)

        #####make a matrix
  interceptmat[is.na(interceptmat)]=0
  interceptmat[fitvec<minval,]=0
  interceptmat[,fitvec<minval]=0

         #####start with raw matrix
  propvals=interceptmat*2
  ourmix=matrix(fitvec,nrow=length(fitvec))
  propvals=(propvals)*(ourmix %*% t(ourmix))

         ###make sums all zero
  test=propvals
  for(i in 1:ncol(test)){
    for(j in 1:ncol(test)){
      propvals[i,j]=test[i,j]-colMeans(test)[j]-rowMeans(test)[i]+mean(test);
    }
  }

  qfitted=eigen(propvals)$vectors[,pcnum]*sqrt(eigen(propvals)$values[pcnum])

  quality=1-sum((qfitted %*% t(qfitted)-propvals)^2)/sum(propvals^2)
  return(c(quality,qfitted))
}

                             ## for inferring dates:
getlhood=function(means,times,admixtimes,popfracs,likpenalty=likpenalty){
	###treat each row of means separately
	###filtering of times done already
	pred=matrix(nrow=length(times),ncol=length(admixtimes))
	for(i in 1:length(admixtimes)){
	vec=exp(-times*admixtimes[i]/100)
	pred[,i]=vec
	}
	
	lhood=0.0

                ## fits means (i.e. admix curves) as a function of 1 or two times (i.e. pred is either a length(times) vector or a length(times)x2 matrix)
	lhoodnegind=0
	to.sum=rep(0,length(admixtimes))
	seq.to.capture=seq(1,nrow(means),by=sqrt(nrow(means)))+0:(sqrt(nrow(means))-1)
	for(i in 1:nrow(means)){

		temp=lm(means[i,]~pred)
		curlhood=-length(temp$res)/2*log(mean((temp$res)^2))
		if (sum(is.element(i,seq.to.capture))==1)
		{
			mean.plus.3sd=3*sd(abs(temp$res[floor(length(temp$res)/2):length(temp$res)]))
			if (length(admixtimes)==2)
			{
				if (admixtimes[1] != admixtimes[2] && !is.na(temp$coeff[3]))
				{
					if ((sign(temp$coeff[2]) != sign(temp$coeff[3])) && (min(temp$coeff[2:3]) < -1*mean.plus.3sd)) lhoodnegind=1
					for (j in 1:length(admixtimes))
			    		    to.sum[j]=to.sum[j]+popfracs[((i-1)/sqrt(nrow(means))+1)]*abs(temp$coeff[(j+1)])
				}
				if (is.na(temp$coeff[3])) lhoodnegind=1
			}
			if (length(admixtimes)==1)
			{
				for (j in 1:length(admixtimes))
				    to.sum[j]=to.sum[j]+popfracs[((i-1)/sqrt(nrow(means))+1)]*abs(temp$coeff[(j+1)])
			}
		}
		lhood=lhood+curlhood
	}
	if (is.na(max(to.sum)))
	{
		to.sum=0
		lhoodnegind=1
	}
	if (lhoodnegind==1 || max(to.sum)>2.0) return(lhood/likpenalty)
	else return(lhood)
	return(lhood)
}

getmax=function(means,times,ourpops="",nparam,popfracs,returnall=F,method="Nelder-Mead"){
	###start with 5 but transform
	x=10^seq(0,3,.25)
	param.mat=matrix(x,ncol=1)
        to.remove=NULL
        for (i in 1:length(x))
          to.remove=c(to.remove,(length(x)*i+1):(length(x)*i+i))
	if (nparam == 2)
	   param.mat=cbind(rep(x,length(x):1),rep(x,length(x))[-to.remove])
	maxlik=0
	initnew=sqrt(param.mat[1,])	
	likpenalty=getlhood(means,times,param.mat[1,],popfracs,1)
	# Longest time
	for (i in 1:dim(param.mat)[1])
	{
		likval=getlhood(means,times,param.mat[i,],popfracs,likpenalty)
		if (likval>maxlik)
		{
			initnew=sqrt(param.mat[i,])
			maxlik=likval			
		}
    }
	
	init=initnew
	ourlhood=function(params){
		newparams=1+params^2
		return(-getlhood(means,times,newparams,popfracs,likpenalty))  
	}
	this=optim(init,ourlhood,method=method)
	likval=this$value
	##transform
	this$par=1+this$par^2    ## maximum parameter value (i.e. date)
	if(returnall==T){       ## (just returns coefficients of linear regression of observed admix times versus expected admix times for date=this$par, in additon to this$par and the "residudal fit" of this linear regression)
		admixtimes=this$par
		pred=matrix(nrow=length(times),ncol=length(admixtimes))
		for(i in 1:length(admixtimes)){
		vec=exp(-times*admixtimes[i]/100)
		pred[,i]=vec
		}
	
		lhood=0.0
		lhoodFULL=0.0
		nulllhood=0.0
		lhood1=1:nrow(means)*0
		lhood2=lhood1
		lhood3=lhood1
		news=matrix(nrow=nrow(means),ncol=(nparam+1))
		y.intercept=matrix(NA,nrow=nrow(means),ncol=nparam)
		for(i in 1:nrow(means)){

			temp=lm(means[i,]~pred)
			y.intercept[i,]=round(temp$coeff[2:length(temp$coeff)],10)
			news[i,]=temp$coeff
			
			ourpar=admixtimes
			lhood2[i]=getlhood(matrix(means[i,],nrow=1),times,ourpar,popfracs[((i-1)/sqrt(nrow(means))+1)],likpenalty)
			quan.val=.1
			fac.val=10
			if ((temp$fitted[length(temp$fitted)] == min(temp$fitted)) && (temp$fitted[1] == max(temp$fitted))) cutoff.val=max((1:length(temp$fitted))[temp$fitted>(temp$fitted[length(temp$fitted)]+(max(temp$fitted)-min(temp$fitted))/fac.val)])
			if ((temp$fitted[length(temp$fitted)] == max(temp$fitted)) && (temp$fitted[1] == min(temp$fitted))) cutoff.val=max((1:length(temp$fitted))[temp$fitted<(temp$fitted[length(temp$fitted)]-(max(temp$fitted)-min(temp$fitted))/fac.val)])
			if (((temp$fitted[length(temp$fitted)] == min(temp$fitted)) && (temp$fitted[1] != max(temp$fitted))) || ((temp$fitted[length(temp$fitted)] != min(temp$fitted)) && (temp$fitted[1] == max(temp$fitted))) || ((temp$fitted[length(temp$fitted)] == max(temp$fitted)) && (temp$fitted[1] != min(temp$fitted))) || ((temp$fitted[length(temp$fitted)] != max(temp$fitted)) && (temp$fitted[1] == min(temp$fitted)))) cutoff.val=length(temp$fitted)
		        lhood1[i]=summary(lm(means[i,]~pred))$r.squared/nrow(means)
			lhood=lhood+lhood1[i];
			lhoodFULL=lhoodFULL+lhood2[i];
			lhood3[i]=getlhood(matrix(means[i,],nrow=1),times,rep(0,length(ourpar)),popfracs,likpenalty)
			nulllhood=nulllhood+lhood3[i];
		}
		if(ourpops[1] !=""){
			indices=(1:length(means))[order(-lhood1+lhood2)][1:min(25,nrow(news))]
			
			par(mfrow=c(5,5))
			for(i in 1:length(indices)){
				predline=news[indices[i],1]
				for(j in 1:nparam) predline=predline+pred[,j]*news[indices[i],(j+1)]
			}
		}

		return(list(par=admixtimes,y.intercept=y.intercept,coeffs=news,lhood=lhood,lhoodFULL=lhoodFULL,lhoodall=nrow(means)*lhood1))
		
	}
	else return(this)
}
prerun=function(temp,donor.pops.all2,num.bins,id.file,copyvector.file){
	id.mat=read.table(id.file,as.is=TRUE)
	num_target =  length(intersect(which(id.mat[,2]==recipient.pops),which(id.mat[,3]=='1')))
	 
	                       ## assign labels (in painting) to each donor group:
	donor.label.vec=rep(-9,ploidy*length(id.mat[,1]))
	donor.pops.all.temp=donor.pops.all
	if (length((1:length(donor.pops.all))[donor.pops.all==simname])>0)
	donor.pops.all.temp=donor.pops.all.temp[-(1:length(donor.pops.all))[donor.pops.all==simname]]
	for (i in 1:length(donor.pops.all.temp))
	{
		id.labels.i=(1:length(id.mat[,1]))[id.mat[,2]==donor.pops.all.temp[i] & id.mat[,3]!=0]
		if (ploidy==2) id.labels.i=sort(c(2*id.labels.i,2*id.labels.i-1))
		donor.label.vec[id.labels.i]=i
	}
	id.labels.i=(1:length(id.mat[,1]))[id.mat[,2]==simname & id.mat[,3]!=0]
	if (ploidy==2) id.labels.i=sort(c(2*id.labels.i,2*id.labels.i-1))
	donor.label.vec[id.labels.i]=0
	to.remove.id=(1:dim(id.mat)[1])[as.character(id.mat[,3])=='0']
	if (length(to.remove.id)>0) id.mat=id.mat[-to.remove.id,]
	if (length(dim(id.mat))==0)
	{
		print(paste("SOMETHING WRONG WITH ",id.file," -- NO NON-EXCLUDED INDS? Exiting....",sep=''))
		q(save='no')
	}
	recipient.label.vec=id.mat[,1][id.mat[,2]==simname]

	                      ## GET RAW COPY PROPS FOR ALL DONORS:
	probs=read.table(copyvector.file,header=TRUE)
	rownames.copyvector=as.character(probs[,1])
	colnames.copyvector=as.character(read.table(copyvector.file,as.is=TRUE,nrows=1)[-1])
	for (i in 1:length(donor.pops.all))
	{
		if (!is.element(donor.pops.all[i],as.character(id.mat[,2])) && !is.element(donor.pops.all[i],colnames.copyvector))
		{
			print(paste("COPY VECTOR COLUMN LABEL ",donor.pops.all[i]," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
			q(save='no')
		}
	}
	for (i in 1:length(donor.pops.all2))
	{
		if (!is.element(donor.pops.all2[i],as.character(id.mat[,2])) && !is.element(donor.pops.all2[i],rownames.copyvector))
		{
			print(paste("DONOR POPULATION ",donor.pops.all2[i]," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
			q(save='no')
		}
	}
	if (!is.element(simname,as.character(id.mat[,2])) && !is.element(simname,rownames.copyvector))
	{
		print(paste("RECIPIENT POPULATION ",simname," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}


	                                   ## combine columns across copy-vector pops:
	predmat.orig=NULL
	for (i in 1:length(donor.pops.all))
	  {
	    id.labels.i=c(donor.pops.all[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all[i]])
	    match.i=NULL
	    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(colnames.copyvector)==id.labels.i[j]))
	    #match.i=match(id.labels.i,as.character(colnames.copyvector))
	    #match.i=match.i[!is.na(match.i)]
	    if (length(match.i)==0)
	    {
		print(paste("NO INDS OF ",donor.pops.all[i]," FOUND AMONG COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	    }
	    predmat.orig=cbind(predmat.orig,apply(matrix(as.matrix(probs[,2:dim(probs)[2]][,match.i]),nrow=dim(probs)[1]),1,sum))
	  }
	rownames(predmat.orig)=rownames.copyvector
	colnames(predmat.orig)=donor.pops.all

	                                   ## combine rows across donor pops:
	predmat=NULL
	for (i in 1:length(donor.pops.all2))
	  {
	    id.labels.i=c(donor.pops.all2[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all2[i]])
	    match.i=NULL
	    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(rownames.copyvector)==id.labels.i[j]))
	    #match.i=match(id.labels.i,as.character(rownames.copyvector))
	    #match.i=match.i[!is.na(match.i)]
	    if (length(match.i)==0)
	    {
	#print(id.labels.i[j])

		print(paste("NO INDS OF ",donor.pops.all2[i]," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	    }
	    predmat=rbind(predmat,apply(matrix(as.matrix(predmat.orig[match.i,]),ncol=dim(predmat.orig)[2]),2,mean))
	  }
	rownames(predmat)=donor.pops.all2
	colnames(predmat)=donor.pops.all

	                       ## GET RAW COPY PROPS FOR RECIPIENT:
	id.labels.rec=c(simname,as.character(id.mat[,1])[as.character(id.mat[,2])==simname])
	match.rec=NULL
	for (j in 1:length(id.labels.rec)) match.rec=c(match.rec,which(as.character(rownames.copyvector)==id.labels.rec[j]))
	#match.rec=match(id.labels.rec,as.character(rownames.copyvector))
	#match.rec=match.rec[!is.na(match.rec)]
	if (length(match.rec)==0)
	{
		print(paste("NO INDS OF ",simname," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
	raw.copyprops.all=apply(matrix(as.matrix(predmat.orig[match.rec,]),ncol=dim(predmat.orig)[2]),2,mean)
	names(raw.copyprops.all)=donor.pops.all

	                        ## remove self and standardize:
	donor.pops.all.orig=donor.pops.all
	donor.pops.all=donor.pops.all[!(colnames(predmat) %in% simname)]
	donor.pops.all2=donor.pops.all2[!(rownames(predmat) %in% simname)]
	predmat=predmat[,!(colnames(predmat) %in% simname)]
	predmat=predmat[!(rownames(predmat) %in% simname),]
	for(i in 1:nrow(predmat)) predmat[i,]=predmat[i,]/sum(predmat[i,])
	raw.copyprops.all=raw.copyprops.all[!(names(raw.copyprops.all) %in% simname)]
	raw.copyprops.all=raw.copyprops.all/sum(raw.copyprops.all)
 ndonor = length(predmat[1,])
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

cor_col=matrix(0,nrow=dim(predmat)[2],ncol=dim(predmat)[2])

for (i in 1:dim(predmat)[2])
{
	for (j in 1:dim(predmat)[2])
	{
	    cor_col[i,j] = euc.dist(predmat[,i],predmat[,j])
    }
}
temp_keep = NULL
for (i in 1:dim(predmat)[2])
{
		temp_keep = append(temp_keep, list(donor.pops.all[which(abs(scale(cor_col[,i])- min(scale(cor_col[,i]))) < 0.2)]))
}

to_skip = NULL
final_col = NULL
for (i in 1:dim(predmat)[2])
{   
	if (i %in% to_skip){next}
    if (length(unlist(temp_keep[i]))>=2)
	{
		final_col = append(final_col,list(unlist(temp_keep[i])[which(!(which(donor.pops.all %in% unlist(temp_keep[i])) %in% to_skip) == TRUE)]))
		to_skip = append(to_skip,which(donor.pops.all %in% unlist(temp_keep[i]))[2:length(unlist(temp_keep[i]))])
	}
	else
	{
		final_col = append(final_col,temp_keep[i])
	}
	
}

final_col = unique(final_col)
keep = final_col
col_temp_keep = keep
donor.label.vec_old = donor.label.vec

ref = unique(keep)
for (i in 1:length(ref))
{
	pop_temp = unlist(ref[i])
	temp =  which(donor.pops.all %in% pop_temp)
	for (j in 1:length(temp))
	{
		donor.label.vec[donor.label.vec_old == temp[j]] = rep(i, length(donor.label.vec_old[donor.label.vec_old == temp[j]]))
	}
}
for  (i in 1:length(ref)){
	ref[i] =  paste(unlist(ref[i]), collapse = '_')
}


                   ## combine rows across donor pops:
				  predmat_newcol = NULL
				  raw.copyprops.all_newcol =NULL
				for (i in 1:length(keep))
			   	  {
			   		  col_temp = NULL
			   		   for (j in 1:length(unlist(keep[i])))
					   {
			   			   col_temp = c(col_temp,which(donor.pops.all == unlist(keep[i])[j]))
			   		   }
			   		   predmat_newcol = cbind(predmat_newcol, rowMeans(as.matrix(predmat[,col_temp])))
			   		   raw.copyprops.all_newcol = c(raw.copyprops.all_newcol, colMeans(as.matrix(raw.copyprops.all[col_temp])))
			   	  }	
	
raw.copyprops.all = raw.copyprops.all_newcol
predmat = predmat_newcol
colnames(predmat)=ref
names(raw.copyprops.all)=ref
print(ref)

mem1 = num_target*ndonor*(ndonor*num.bins*8/10^9)
ndonor2 = length(predmat[1,])
mem2 = num_target*ndonor2*(ndonor2*num.bins*8/10^9)
print(paste("number of targets:",num_target))
if (num.bins<= 300) {print(paste("number of bins:",num.bins))}
if (num.bins>300) {print(paste0("number of bins:",num.bins,". ***In the parameter file, we suggest to use 'curve.range: 1 30' which is enough for making inference."))}
print(paste("number of donors in mode1:",ndonor))
print(paste("number of donors after merging in mode2:",ndonor2))
print(paste('Mode 1 requires memory of',round(mem1+1,1), 'G.'))
print(paste('Mode 2 requires memory of',round(mem2+1,1), 'G (may provide less precise inference).'))
print(paste('Mode 3 generally requires much less memory but can take longer time.'))
print(paste('(Note that, for larger estimated memory usage, these estimates may be off by 1G or so.)'))
quit()
}




##############################################################
#############################################################
### (I)-(a) START READING IN DATA:
##############################################################
#############################################################

donor.pops.all2=unique(newnames)
simname=recipient.pops
if (temp[5] == 'mem'){prerun(temp,donor.pops.all2,num.bins,id.file,copyvector.file)}		
		if (length(temp) < 5) 
		   {
		   	print('Please specify running mode whether 1,2 or 3. Now using default as mode 1.')
			temp[5] = '1'
		   }
		if (temp[5] != '1' && temp[5] != '2' && temp[5] != '3') 
		    {
			print('Please specify running mode whether 1,2 or 3. Now using default as mode 1.')
			temp[5] = '1'
		    }
			mode = temp[5]
	
                        ## GET IDS FILE:
id.mat=read.table(id.file,as.is=TRUE)
#to.remove.id=(1:dim(id.mat)[1])[as.character(id.mat[,3])!='R' | as.character(id.mat[,3])!='r' | as.character(id.mat[,3])!='D' | as.character(id.mat[,3])!='d']
                       ## assign labels (in painting) to each donor group:
donor.label.vec=rep(-9,ploidy*length(id.mat[,1]))
donor.pops.all.temp=donor.pops.all
if (length((1:length(donor.pops.all))[donor.pops.all==simname])>0)
donor.pops.all.temp=donor.pops.all.temp[-(1:length(donor.pops.all))[donor.pops.all==simname]]
for (i in 1:length(donor.pops.all.temp))
{
	id.labels.i=(1:length(id.mat[,1]))[id.mat[,2]==donor.pops.all.temp[i] & id.mat[,3]!=0]
	if (ploidy==2) id.labels.i=sort(c(2*id.labels.i,2*id.labels.i-1))
	donor.label.vec[id.labels.i]=i
}
id.labels.i=(1:length(id.mat[,1]))[id.mat[,2]==simname & id.mat[,3]!=0]
if (ploidy==2) id.labels.i=sort(c(2*id.labels.i,2*id.labels.i-1))
donor.label.vec[id.labels.i]=0
to.remove.id=(1:dim(id.mat)[1])[as.character(id.mat[,3])=='0']
if (length(to.remove.id)>0) id.mat=id.mat[-to.remove.id,]
if (length(dim(id.mat))==0)
{
	print(paste("SOMETHING WRONG WITH ",id.file," -- NO NON-EXCLUDED INDS? Exiting....",sep=''))
	q(save='no')
}
recipient.label.vec=id.mat[,1][id.mat[,2]==simname]

                      ## GET RAW COPY PROPS FOR ALL DONORS:
probs=read.table(copyvector.file,header=TRUE)
rownames.copyvector=as.character(probs[,1])
colnames.copyvector=as.character(read.table(copyvector.file,as.is=TRUE,nrows=1)[-1])
for (i in 1:length(donor.pops.all))
{
	if (!is.element(donor.pops.all[i],as.character(id.mat[,2])) && !is.element(donor.pops.all[i],colnames.copyvector))
	{
		print(paste("COPY VECTOR COLUMN LABEL ",donor.pops.all[i]," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
}
for (i in 1:length(donor.pops.all2))
{
	if (!is.element(donor.pops.all2[i],as.character(id.mat[,2])) && !is.element(donor.pops.all2[i],rownames.copyvector))
	{
		print(paste("DONOR POPULATION ",donor.pops.all2[i]," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
		q(save='no')
	}
}
if (!is.element(simname,as.character(id.mat[,2])) && !is.element(simname,rownames.copyvector))
{
	print(paste("RECIPIENT POPULATION ",simname," NOT FOUND IN DONORS/RECIPIENTS OF ",id.file," OR ROWS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
}
#match.colnames=match(as.character(colnames.copyvector),as.character(id.mat[,1]))
#match.colnamesII=match(as.character(colnames.copyvector),as.character(id.mat[,2]))
#for (i in 1:length(colnames.copyvector))
#{
#	if (is.na(match.colnames[i]) && is.na(match.colnamesII[i]))
#	{
#		print(paste("COLUMN LABEL ",colnames.copyvector[i]," IN ",copyvector.file," IS NOT FOUND IN EITHER COLUMN OF ",id.file,"! Exiting....",sep=''))
#		q(save='no')
#	}
#}
#match.rownames=match(as.character(rownames.copyvector),as.character(id.mat[,1]))
#if (sum(is.na(match.rownames))!=0)
#{
#	print(paste("SOME ROW LABELS IN ",copyvector.file," ARE NOT FOUND IN THE FIRST COLUMN OF ",id.file,"! Exiting....",sep=''))
#	q(save='no')
#}
#for (i in 1:length(donor.pops.all2))
#{
#	if (!is.element(donor.pops.all2[i],as.character(id.mat[,2][match.rownames])))
#	{
#		print(paste("DONOR POPULATION LABEL ",donor.pops.all[i]," NOT FOUND IN ROWS OF",copyvector.file,"! Exiting....",sep=''))
#		q(save='no')
#	}
#}

                                   ## combine columns across copy-vector pops:
predmat.orig=NULL
for (i in 1:length(donor.pops.all))
  {
    id.labels.i=c(donor.pops.all[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all[i]])
    match.i=NULL
    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(colnames.copyvector)==id.labels.i[j]))
    #match.i=match(id.labels.i,as.character(colnames.copyvector))
    #match.i=match.i[!is.na(match.i)]
    if (length(match.i)==0)
    {
	print(paste("NO INDS OF ",donor.pops.all[i]," FOUND AMONG COLUMNS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
    }
    predmat.orig=cbind(predmat.orig,apply(matrix(as.matrix(probs[,2:dim(probs)[2]][,match.i]),nrow=dim(probs)[1]),1,sum))
  }
rownames(predmat.orig)=rownames.copyvector
colnames(predmat.orig)=donor.pops.all

                                   ## combine rows across donor pops:
predmat=NULL
for (i in 1:length(donor.pops.all2))
  {
    id.labels.i=c(donor.pops.all2[i],as.character(id.mat[,1])[as.character(id.mat[,2])==donor.pops.all2[i]])
    match.i=NULL
    for (j in 1:length(id.labels.i)) match.i=c(match.i,which(as.character(rownames.copyvector)==id.labels.i[j]))
    #match.i=match(id.labels.i,as.character(rownames.copyvector))
    #match.i=match.i[!is.na(match.i)]
    if (length(match.i)==0)
    {
#print(id.labels.i[j])

	print(paste("NO INDS OF ",donor.pops.all2[i]," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
    }
    predmat=rbind(predmat,apply(matrix(as.matrix(predmat.orig[match.i,]),ncol=dim(predmat.orig)[2]),2,mean))
  }
rownames(predmat)=donor.pops.all2
colnames(predmat)=donor.pops.all

                       ## GET RAW COPY PROPS FOR RECIPIENT:
id.labels.rec=c(simname,as.character(id.mat[,1])[as.character(id.mat[,2])==simname])
match.rec=NULL
for (j in 1:length(id.labels.rec)) match.rec=c(match.rec,which(as.character(rownames.copyvector)==id.labels.rec[j]))
#match.rec=match(id.labels.rec,as.character(rownames.copyvector))
#match.rec=match.rec[!is.na(match.rec)]
if (length(match.rec)==0)
{
	print(paste("NO INDS OF ",simname," FOUND AMONG ROWS OF ",copyvector.file,"! Exiting....",sep=''))
	q(save='no')
}
raw.copyprops.all=apply(matrix(as.matrix(predmat.orig[match.rec,]),ncol=dim(predmat.orig)[2]),2,mean)
names(raw.copyprops.all)=donor.pops.all

                        ## remove self and standardize:
donor.pops.all.orig=donor.pops.all
donor.pops.all=donor.pops.all[!(colnames(predmat) %in% simname)]
donor.pops.all2=donor.pops.all2[!(rownames(predmat) %in% simname)]
predmat=predmat[,!(colnames(predmat) %in% simname)]
predmat=predmat[!(rownames(predmat) %in% simname),]
for(i in 1:nrow(predmat)) predmat[i,]=predmat[i,]/sum(predmat[i,])
raw.copyprops.all=raw.copyprops.all[!(names(raw.copyprops.all) %in% simname)]
raw.copyprops.all=raw.copyprops.all/sum(raw.copyprops.all)

##############################################################
#############################################################
### (I)-(b) IF REQUESTED, ONLY INFER PROPORTIONS OF ANCESTRY: (I.E. AS IN LESLIE ET AL 2015, NATURE)
##############################################################
#############################################################

if ((prop.ind==1 || prop.continue.ind==1) && num.iterations==0)
{
        ourmix=getoverallfit(predmat,raw.copyprops.all)$x
     	ourmix=ourmix[ourmix>0]
     	ourmix=ourmix[ourmix>prop.cutoff]
     	#ourmix=ourmix/sum(ourmix)
     	donor.pops=names(ourmix)
	donor.pops=donor.pops[sort(ourmix,index.return=TRUE)$ix]
	ourmix=ourmix[sort(ourmix,index.return=TRUE)$ix]
	write("### INFERRED ANCESTRY PROPORTIONS (MIXING COEFFICIENTS), WITHOUT INFERRING ANY ADMIXTURE EVENTS:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
     	write.table(rbind(donor.pops,ourmix),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
	print("Finished!")
	q(save='no')
}

##############################################################
#############################################################
### (I)-(c) OTHERWISE CONTINUE READING IN:
##############################################################
#############################################################

samples.filein=as.character(temp[3])
recomrates.filein=as.character(temp[4])


                        ## find number of samples files:
num.samplefiles=0
readfile=file(samples.filein,open="r")
line2=1
while(!is.na(line2[1]))
{
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
    	if (!is.na(line2[1])) num.samplefiles=num.samplefiles+1
}
if (num.samplefiles==0)
{
	print(paste("SOMETHING WRONG WITH INPUT FILE ",samples.filein," (PERHAPS EMPTY?). Exiting....",sep=''))
	q(save='no')
}
close(readfile)

                        ## find number of recom files:
num.recomfiles=0
readfile=file(recomrates.filein,open="r")
line2=1
while(!is.na(line2[1]))
{
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
    	if (!is.na(line2[1])) num.recomfiles=num.recomfiles+1
}
if (num.recomfiles==0)
{
	print(paste("SOMETHING WRONG WITH INPUT FILE ",recomrates.filein," (PERHAPS EMPTY?). Exiting....",sep=''))
	q(save='no')
}
close(readfile)

                         ## get number of individuals:
readfile=file(as.character(read.table(samples.filein,nrow=1,as.is=TRUE)),open="r")
line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
nsamples=as.integer(strsplit(as.character(line2[21]),split=",")[[1]][1])
if (is.na(nsamples) || nsamples<=0)
{
	print(paste("No painting samples in ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE)),"! Exiting....",sep=''))
	q(save='no')
}
count=0
rec.hap.count=0
while(!is.na(line2[1]))
{
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
    	if (!is.na(line2[1]))
	{
		if (as.character(line2[1])=="HAP")
		{
			if (is.element(line2[3],recipient.label.vec)) rec.hap.count=rec.hap.count+1
		}
		count=count+1
	}
}
close(readfile)
ninds.tot=count/(ploidy*(nsamples+1))

if (ninds.tot<=0)
{
	print(paste("No individuals in ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE)),"! Exiting....",sep=''))
	q(save='no')
}
ninds=rec.hap.count/ploidy
if (ninds != floor(rec.hap.count/ploidy))
{
	print(paste("Number of target individuals in ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE))," does not match ploidy! Exiting....",sep=''))
	q(save='no')
}
if (ninds<=0)
{
	print(paste("No target individuals in ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE)),"! Exiting....",sep=''))
	q(save='no')
}
if (ninds<length(recipient.label.vec))
{
	print(paste("Number of target individuals in",id.file," and ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE))," do not match. Exiting....",sep=''))
	q(save='no')
}
if (length(recipient.label.vec) != ninds) 
{
	print(paste("Warning!!! Number of target individuals in",id.file," and ",as.character(read.table(samples.filein,nrow=1,as.is=TRUE))," do not match.",sep=''))
}

ind.id.vec=matrix(rep(1:ninds-1,each=num.samplefiles),byrow=FALSE,nrow=num.samplefiles)

##############################################################
#############################################################
### (II) SOURCE/PROPORTION/DATE ESTIMATION
##############################################################
#############################################################
#Start surrogate grouping here
#save.image("test.RData")
if (mode == 2){


euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

cor_col=matrix(0,nrow=dim(predmat)[2],ncol=dim(predmat)[2])

for (i in 1:dim(predmat)[2])
{
	for (j in 1:dim(predmat)[2])
	{
	    cor_col[i,j] = euc.dist(predmat[,i],predmat[,j])
    }
}
temp_keep = NULL
for (i in 1:dim(predmat)[2])
{
		temp_keep = append(temp_keep, list(donor.pops.all[which(abs(scale(cor_col[,i])- min(scale(cor_col[,i]))) < 0.2)]))
}

to_skip = NULL
final_col = NULL
for (i in 1:dim(predmat)[2])
{   
	if (i %in% to_skip){next}
    if (length(unlist(temp_keep[i]))>=2)
	{
		final_col = append(final_col,list(unlist(temp_keep[i])[which(!(which(donor.pops.all %in% unlist(temp_keep[i])) %in% to_skip) == TRUE)]))
		to_skip = append(to_skip,which(donor.pops.all %in% unlist(temp_keep[i]))[2:length(unlist(temp_keep[i]))])
	}
	else
	{
		final_col = append(final_col,temp_keep[i])
	}
	
}

final_col = unique(final_col)
keep = final_col
#print(final_col)
col_temp_keep = keep
donor.label.vec_old = donor.label.vec

ref = unique(keep)
for (i in 1:length(ref))
{
	pop_temp = unlist(ref[i])
	temp =  which(donor.pops.all %in% pop_temp)
	for (j in 1:length(temp))
	{
		donor.label.vec[donor.label.vec_old == temp[j]] = rep(i, length(donor.label.vec_old[donor.label.vec_old == temp[j]]))
	}
}
for  (i in 1:length(ref)){
	ref[i] =  paste(unlist(ref[i]), collapse = '_')
}


                   ## combine rows across donor pops:
				  predmat_newcol = NULL
				  raw.copyprops.all_newcol =NULL
				for (i in 1:length(keep))
			   	  {
			   		  col_temp = NULL
			   		   for (j in 1:length(unlist(keep[i])))
					   {
			   			   col_temp = c(col_temp,which(donor.pops.all == unlist(keep[i])[j]))
			   		   }
			   		   predmat_newcol = cbind(predmat_newcol, rowMeans(as.matrix(predmat[,col_temp])))
			   		   raw.copyprops.all_newcol = c(raw.copyprops.all_newcol, colMeans(as.matrix(raw.copyprops.all[col_temp])))
			   	  }	
	
raw.copyprops.all = raw.copyprops.all_newcol
predmat = predmat_newcol
colnames(predmat)=ref
donor.pops.all = ref
names(raw.copyprops.all)=ref
print(ref)
### End of grouping 
}

if (prop.ind==1)
{
                       ## GET INITIAL FIT:
        ourmix=getoverallfit(predmat,raw.copyprops.all)$x
     	ourmix=ourmix[ourmix>0]
     	ourmix=ourmix[ourmix>prop.cutoff]
     	ourmix=ourmix/sum(ourmix)
     	donor.pops=names(ourmix)
	if (length(donor.pops)==1)
	{
		print(donor.pops)
		print("A SINGLE DONOR POP FITS MIXTURE AFTER 0 ITERATIONS OF SOURCE/DATE ESTIMATION! Exiting early....")
		write("### INFERRED SOURCES AND DATES ('best-guess' conclusion: Unknown",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT EVIDENCE, DATE ESTIMATE, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.1date","proportion.source1","maxR2fit.1date","fit.quality.1event","fit.quality.2events","bestmatch.event1.source1","bestmatch.event1.source2","proportion.event2.source1","bestmatch.event2.source1","bestmatch.event2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
		write(rep(NA,10),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT EVIDENCE, DATE ESTIMATES, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.2dates.date1","gen.2dates.date2","maxScore.2events","proportion.date1.source1","bestmatch.date1.source1","bestmatch.date1.source2","proportion.date2.source1","bestmatch.date2.source1","bestmatch.date2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
		write(rep(NA,9),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(donor.pops,ourmix),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		write.table(rbind(NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC2:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE1-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE2-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
 		q(save='no')
	}
			## make weights:
	weights.mat=matrix(predmat[match(donor.pops,donor.pops.all2),],ncol=dim(predmat)[2])
	rownames(weights.mat)=donor.pops
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]*ourmix
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]/sum(weights.mat[,i])
        #weights.mat2=cbind(weights.mat,0)   ## for self-copying (only necessary for curves)
        #rownames(weights.mat2)=donor.pops
        #weights.mat2[is.nan(weights.mat2)]=0.0;
        newdim=nrow(weights.mat)
        olddim=length(donor.pops.all)
	    tempweights = matrix(mutiply.temp.weight(olddim,newdim,weights.mat),byrow=TRUE,ncol=length(donor.pops.all)*length(donor.pops.all))
		
			## generate (weighted) coancestry curves:
	dists.shift=abs(seq(0,grid.range[2],bin.width)-grid.range[1])
	ptm1  = proc.time()
	if (mode == 1 || mode== 2)
	{
	means=matrix(coancestry.curves(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
	}
	if (mode == 3)
	{
	means=matrix(coancestry.curves.mode3(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
	}
	num_target =  length(intersect(which(id.mat[,2]==recipient.pops),which(id.mat[,3]=='1')))
	
	if (null.ind==1 && num_target == 1)  null.ind = 0 # avoid Null analysis if there is only one target individual
	if (null.ind==1 && null.calc==0)
	{
		null.calc=1
		tempchromindres.null=matrix(coancestry.curves.null(ploidy,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,simname,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff),byrow=TRUE,ncol=num.bins)
    }
	if (null.ind==1)
	{
		expindres=matrix(0,nrow=nrow(weights.mat)^2,ncol=length(dists))
		#save.image("test.RData")
		#print(tempweights)
		#print(tempchromindres.null)
	        results=tempweights %*% tempchromindres.null
		for (k in 1:length(dists)) results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))

	               #####expectations
		sums=rep(0,length(dists))
		sums2=matrix(0, nrow=nrow(weights.mat),ncol=length(dists))
		sums3=sums2
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat))
	      		{
			    sums=sums+results[(k-1)*nrow(weights.mat)+l,]
			    sums2[k,]=sums2[k,]+results[(k-1)*nrow(weights.mat)+l,]
			    sums3[l,]=sums3[l,]+results[(k-1)*nrow(weights.mat)+l,]		
	      		}
		}
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat)) expindres[(k-1)*nrow(weights.mat)+l,]=sums2[k,]*sums3[l,]/sums
		}
		for (k in 1:length(dists))
  		{
			#results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
			expindres[,k]=c(t((matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
 		}
		results=results/expindres
		means=means/results
	}

                 ## infer new date and intercept.mat:
      	#names3=rownames(weights.mat2)
      	names3=rownames(weights.mat)
		#times=dists
		#new.grid.range = findpeak(times,ourmix,means,bin.width) 
        #dists=seq(new.grid.range,grid.range[2],bin.width)
      	times=dists
        pop.comb=paste(rep(names3,each=length(names3)),rep(names3,length(names3)),sep=' ')
        rownames(means)=pop.comb

        check=getmax(means,times,nparam=1,popfracs=ourmix,returnall=T,ourpops=rownames(means))
        admixtimes=check$par
	    pred=matrix(exp(-times*admixtimes/100),ncol=1)
        r.squared.est=r.squared.est.2event=intercept.est=rep(NA,nrow(means))
        for(i in 1:nrow(means))
 	{
            temp=lm(means[i,]~pred)
            r.squared.est[i]=summary(temp)$r.squared
	    intercept.est[i]=round(temp$coeff[2],10)
        }
        if (num.iterations==0) ## NEVER GETS CALLED!
	{
		check2=getmax(means,times,nparam=2,popfracs=ourmix,returnall=T,ourpops=rownames(means))
		admixtimes2=check2$par
        	pred2=matrix(nrow=length(times),ncol=2)
        	for(k in 1:2)
        	{
			vec=exp(-times*admixtimes2[k]/100)
            		pred2[,k]=vec
		}
		intercept.est.2date=matrix(NA,nrow=nrow(means),ncol=2)
        	for(i in 1:nrow(means))
 		{
            		temp2=lm(means[i,]~pred2)
            		r.squared.est.2event[i]=summary(temp2)$r.squared
			intercept.est.2date[i,]=round(temp2$coeff[2:length(temp2$coeff)],10)
        	}
        	intercept.2date.mat1=matrix(intercept.est.2date[,1],ncol=length(donor.pops),byrow=TRUE)
        	intercept.2date.mat2=matrix(intercept.est.2date[,2],ncol=length(donor.pops),byrow=TRUE)	
	}
	r.squared.est.2=r.squared.est
	r.squared.est.2event.2=r.squared.est.2event
        intercept.mat=matrix(intercept.est,ncol=length(donor.pops),byrow=TRUE)
        donor.pops.bestrep=strsplit(pop.comb[intercept.est<0][r.squared.est[intercept.est<0]==max(r.squared.est[intercept.est<0])][1],split=" ")[[1]]
}
if (prop.ind==0)
{
                       ## GET INITIAL FIT:
     props.read=read.table(paste(save.file.props,prop.save.post,sep=''),skip=11,nrows=2,as.is=TRUE)
     mixture.props.1=as.double(props.read[2,2:dim(props.read)[2]])
     names(mixture.props.1)=as.character(props.read[1,2:dim(props.read)[2]])
     prop.val.est=as.double(read.table(paste(save.file.props,prop.save.post,sep=''),skip=3,nrows=1,header=TRUE)$proportion.source1)
     props.read=read.table(paste(save.file.props,prop.save.post,sep=''),skip=13,nrows=2,as.is=TRUE)
     mixture.props.2=as.double(props.read[2,2:dim(props.read)[2]])
     names(mixture.props.2)=as.character(props.read[1,2:dim(props.read)[2]])

     ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops.all2))
     ourmix.pop1[match(names(mixture.props.1),donor.pops.all2)]=mixture.props.1
     ourmix.pop2[match(names(mixture.props.2),donor.pops.all2)]=mixture.props.2
     ourmix=prop.val.est*ourmix.pop1+(1-prop.val.est)*ourmix.pop2
     names(ourmix)=donor.pops.all2
     ourmix=ourmix[ourmix>prop.cutoff]
     ourmix=ourmix/sum(ourmix)
     donor.pops=names(ourmix)
	 #print(donor.pops)
#     if (length(donor.pops)==1)
#     {
#	write("### INFERRED SOURCES AND DATES",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
#	write("### 1-DATE FIT EVIDENCE, DATE ESTIMATE, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write(c("gen.1date","maxR2fit.1date","fit.quality.1event","fit.quality.2events","bestmatch.event1.source1","bestmatch.event1.source2","bestmatch.event2.source1","bestmatch.event2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=8,append=TRUE)
#	write(rep(NA,8),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=8,append=TRUE)
#	write("### 2-DATE FIT EVIDENCE, DATE ESTIMATES, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write(c("gen.2dates.date1","gen.2dates.date2","maxScore.2events","bestmatch.date1.source1","bestmatch.date1.source2","bestmatch.date2.source1","bestmatch.date2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=7,append=TRUE)
#	write(rep(NA,7),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=8,append=TRUE)
#	write("### 1-DATE FIT SOURCES, PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write.table(rbind(c("proportion",donor.pops),c(NA,ourmix)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write.table(rbind(c("fit.quality"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write("### 1-DATE FIT SOURCES, PC2:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write.table(rbind(c("proportion"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write.table(rbind(c("fit.quality"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write("### 2-DATE FIT SOURCES, DATE1-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write.table(rbind(c("proportion"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write.table(rbind(c("fit.quality"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write("### 2-DATE FIT SOURCES, DATE2-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
#	write.table(rbind(c("proportion"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	write.table(rbind(c("fit.quality"),c(NA)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
#	q(save='no')
#     }

			## get date/intercepts:
     admixtimes=as.double(read.table(paste(save.file.props,prop.save.post,sep=''),skip=3,nrows=1,header=TRUE)$gen.1date)
     admixtimes2=c(as.double(read.table(paste(save.file.props,prop.save.post,sep=''),skip=6,nrows=1,header=TRUE)$gen.2dates.date1),as.double(read.table(paste(save.file.props,prop.save.post,sep=''),skip=5,nrows=1,header=TRUE)$gen.2dates.date2))
     x=read.table(paste(save.file.props,"_curves.txt",sep=''),skip=1,header=TRUE,check.names=FALSE)
     ##dists=as.double(names(x)[7:length(names(x))])  ### !!!! error with fastGT (as it can make new dists per each run) !!!!
     times=dists
     intercept.est=as.double(x$intercept.fit)[seq(2,dim(x)[1],4)]
     intercept.est.fitted=as.double(x$intercept.fit)[seq(3,dim(x)[1],4)]
     intercept.est.2date.date1=as.double(x$intercept.fit)[seq(4,dim(x)[1],4)]
     intercept.est.2date.date2=as.double(x$intercept.fit.date2)[seq(4,dim(x)[1],4)]
     donors.infile=unique(as.character(x$surrogate1))
     match.check=match(donor.pops,donors.infile)
     if (sum(is.na(match.check))>0)
     {
	print(paste("THE FOLLOWING DONORS CANNOT BE FOUND IN USER-SUPPLIED FILE ",paste(save.file.props,"_curves.txt",sep=''),": ",setdiff(donor.pops,donors.infile),". PERHAPS YOU HAVE LOWERED YOUR THRESHOLD VALUE OF 'props.cutoff' IN ",param.infile," FROM BEFORE? Exiting....",sep=''))
	q(save='no')
     }
     donorpop.comb=NULL
     for (i in 1:length(donor.pops))
     	 donorpop.comb=c(donorpop.comb,paste(donor.pops[i]," ",donor.pops,sep=''))
     pop.comb.orig=paste(as.character(x$surrogate1)," ",as.character(x$surrogate2),sep='')[seq(1,dim(x)[1],4)]
     intercept.mat=matrix(intercept.est[match(donorpop.comb,pop.comb.orig)],ncol=length(donor.pops))
     intercept.mat.fitted=matrix(intercept.est.fitted[match(donorpop.comb,pop.comb.orig)],ncol=length(donor.pops))
     intercept.2date.mat1=matrix(intercept.est.2date.date1[match(donorpop.comb,pop.comb.orig)],ncol=length(donor.pops),byrow=TRUE)
     intercept.2date.mat2=matrix(intercept.est.2date.date2[match(donorpop.comb,pop.comb.orig)],ncol=length(donor.pops),byrow=TRUE)	
     colnames(intercept.mat)=donor.pops
     rownames(intercept.mat)=donor.pops
     colnames(intercept.mat.fitted)=donor.pops
     rownames(intercept.mat.fitted)=donor.pops
     colnames(intercept.2date.mat1)=donor.pops
     rownames(intercept.2date.mat1)=donor.pops
     colnames(intercept.2date.mat2)=donor.pops
     rownames(intercept.2date.mat2)=donor.pops

                         ## get R^2, etc:
     r.squared.est=as.double(x$rsquared.date.fit)[seq(2,dim(x)[1],4)][match(donorpop.comb,pop.comb.orig)]
     r.squared.est.2event=as.double(x$rsquared.date.fit)[seq(4,dim(x)[1],4)][match(donorpop.comb,pop.comb.orig)]
     r.squared.est.2=r.squared.est
     r.squared.est.2event.2=r.squared.est.2event
     donor.pops.bestrep=strsplit(pop.comb.orig[intercept.est<0][r.squared.est[intercept.est<0]==max(r.squared.est[intercept.est<0])][1],split=" ")[[1]]
     pop.comb=donorpop.comb

                        ## checks:
     if (is.na(prop.val.est) || is.na(admixtimes) || sum(is.na(admixtimes2))>0 || sum(is.na(mixture.props.1))>0 || sum(is.na(mixture.props.2)) || sum(is.na(ourmix))>0)
     {
	print(paste("SOMETHING WRONG WITH INPUT FILE ",paste(save.file.props,prop.save.post,sep='')," -- NA VALUES IN REQUIRED DONOR POPS. Exiting....",sep=''))
	q(save='no')
     }
     if (sum(is.na(c(intercept.mat)))>0 || sum(is.na(c(intercept.mat.fitted)))>0 || sum(is.na(c(intercept.2date.mat1)))>0 || sum(is.na(c(intercept.2date.mat2)))>0 || sum(is.na(r.squared.est))>0 || sum(is.na(r.squared.est.2event))>0)
     {
	print(paste("SOMETHING WRONG WITH INPUT FILE ",paste(save.file.props,"_curves.txt",sep='')," -- NA VALUES IN REQUIRED DONOR POPS. Exiting....",sep=''))
	q(save='no')
     }

			## make weights.mat:
     weights.mat=matrix(predmat[match(donor.pops,donor.pops.all2),],ncol=dim(predmat)[2])
     rownames(weights.mat)=donor.pops
     for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]*ourmix
     for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]/sum(weights.mat[,i])
     #weights.mat2=cbind(weights.mat,0)   ## for self-copying (only necessary for curves)
     #rownames(weights.mat2)=donor.pops
     #weights.mat2[is.nan(weights.mat2)]=0.0;
     names3=rownames(weights.mat)

     newdim=nrow(weights.mat)
     olddim=length(donor.pops.all)
     tempweights = matrix(mutiply.temp.weight(olddim,newdim,weights.mat),byrow=TRUE,ncol=length(donor.pops.all)*length(donor.pops.all))
  
      if (prop.continue.ind==1 && num.iterations==1)
      {
			## generate (weighted) coancestry curves:
	  dists.shift=abs(seq(0,grid.range[2],bin.width)-grid.range[1])
  	if (mode == 1 || mode== 2)
  	{
  	means=matrix(coancestry.curves(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
  	}
  	if (mode == 3)
  	{
  	means=matrix(coancestry.curves.mode3(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
  	}
          pop.comb=paste(rep(names3,each=length(names3)),rep(names3,length(names3)),sep=' ')
          rownames(means)=pop.comb
      	  if (null.ind==1 && null.calc==0)
      	  {
		null.calc=1
	  	tempchromindres.null=matrix(coancestry.curves.null(ploidy,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,simname,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff),byrow=TRUE,ncol=num.bins)
          }
      	  if (null.ind==1)
      	  {
		expindres=matrix(0,nrow=nrow(weights.mat)^2,ncol=length(dists))
	   	results=tempweights %*% tempchromindres.null
	   	for (k in 1:length(dists)) results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))

	               #####expectations
	        sums=rep(0,length(dists))
	   	sums2=matrix(0, nrow=nrow(weights.mat),ncol=length(dists))
	   	sums3=sums2
	   	for(k in 1:nrow(weights.mat))
	   	{
			for(l in 1:nrow(weights.mat))
	      		{
				sums=sums+results[(k-1)*nrow(weights.mat)+l,]
				sums2[k,]=sums2[k,]+results[(k-1)*nrow(weights.mat)+l,]
				sums3[l,]=sums3[l,]+results[(k-1)*nrow(weights.mat)+l,]		
	      		}
	   	}
	   	for(k in 1:nrow(weights.mat))
	   	{
			for(l in 1:nrow(weights.mat)) expindres[(k-1)*nrow(weights.mat)+l,]=sums2[k,]*sums3[l,]/sums
	   	}
	   	for (k in 1:length(dists))
  	   	{
			#results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
			expindres[,k]=c(t((matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
 	   	}
	   	results=results/expindres
	   	means=means/results
	  }
      }
}
if ((prop.ind==1 || prop.continue.ind==1) && num.iterations>0)
{
     pca.frac.val=1.0    ## shrink PC intercepts by this factor
     sd.fac=1/20         ## for proportion CI
     alpha.mat=NULL
     alpha.vec=seq(0.01,0.99,0.01)
     
     for (a in 1:num.iterations)
     { 
		 ptm4 = proc.time()
        print(paste("Iteration ",a," of ",num.iterations,"....",sep=''))
	get.fitted.results=getfitted(weights.mat,ourmix,intercept.mat,pcnum=1)
 	intercepts.pca=get.fitted.results[2:length(get.fitted.results)]
  	fit.quality.orig=get.fitted.results[1]
         ##############who admixed?
  	results=fitwhoadmixed(weights.mat,predmat,intercepts.pca,raw.copyprops.all)
	results.fit=results

                               ## EVALUATE 1-DATE MULTIPLE-EVENTS AND 2-DATE EVENTS IF WE'VE DONE ALL ITERATIONS:
       
		if (a==num.iterations)
        {
	        ptm6 = proc.time()
                            ## if necessary, get second date:
							
			ptm9 = proc.time()
            check2=getmax(means,times,nparam=2,popfracs=ourmix,returnall=T,ourpops=rownames(means))
            admixtimes2=check2$par
            pred2=matrix(nrow=length(times),ncol=2)
	    for(k in 1:2)
            {
                  vec2=exp(-times*admixtimes2[k]/100)
                  pred2[,k]=vec2
            }
                           ## get intercepts for each date:
            r.squared.est.2event=rep(NA,nrow(means))
            intercept.est.2date=matrix(NA,nrow=nrow(means),ncol=2)
            for(i in 1:nrow(means))
            {
                  temp2=lm(means[i,]~pred2)
                  r.squared.est.2event[i]=summary(temp2)$r.squared
                  intercept.est.2date[i,]=round(temp2$coeff[2:length(temp2$coeff)],10)
            }
            intercept.2date.mat1=matrix(intercept.est.2date[,1],ncol=length(donor.pops),byrow=TRUE)
            intercept.2date.mat2=matrix(intercept.est.2date[,2],ncol=length(donor.pops),byrow=TRUE)
                         ## assuming 1-date with multiple events, get mixing info on second event:
            fit.quality.1date.PC2=0
            if (length(ourmix)>2)
	    { 	
		get.fitted.results.1date.PC2=getfitted(weights.mat,ourmix,intercept.mat,pcnum=2)
#            fit.quality2=get.fitted.results.1date.PC2[1]
 	    	intercepts.pca.1date.PC2=get.fitted.results.1date.PC2[2:length(get.fitted.results.1date.PC2)]
  	    	fit.quality.1date.PC2=get.fitted.results.1date.PC2[1]
	    }
	    if (is.na(fit.quality.1date.PC2)) fit.quality.1date.PC2=0
	    if (fit.quality.1date.PC2==0 || length(ourmix)<=2)
	    {
		prop.val.est.1date.PC2=NA
     		prop.val.est.ciII.1date.PC2=NA
		mse.score.1date.PC2=NA
		mixture.props.1.1date.PC2=NA
		names(mixture.props.1.1date.PC2)=NA
		mixture.props.2.1date.PC2=NA	
		names(mixture.props.2.1date.PC2)=NA
     	    }
	    if (fit.quality.1date.PC2>0)
	    {
		results.1date.PC2=fitwhoadmixed(weights.mat,predmat,pca.frac.val*intercepts.pca.1date.PC2,raw.copyprops.all)
  	    	prop.val.est.1date.PC2=(results.1date.PC2[names(results.1date.PC2)=="alpha"][which(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"])== min(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]),na.rm=TRUE))])$alpha
  	    	mse.score.all.1date.PC2=as.double(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]))
  	    	mse.score.1date.PC2=(results.1date.PC2[names(results.1date.PC2)=="avrelmse"][which(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"])== min(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]),na.rm=TRUE))])$avrelmse
  	    	mu.fitted.1date.PC2=(results.1date.PC2[names(results.1date.PC2)=="mufitted"][which(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"])== min(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]),na.rm=TRUE))])$mufitted
  	    	mixture.props.1.1date.PC2=sort((results.1date.PC2[names(results.1date.PC2)=="pop1fit"][which(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"])== min(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]),na.rm=TRUE))])$pop1fit)
  	    	mixture.props.2.1date.PC2=sort((results.1date.PC2[names(results.1date.PC2)=="pop2fit"][which(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"])== min(unlist(results.1date.PC2[names(results.1date.PC2)=="avrelmse"]),na.rm=TRUE))])$pop2fit)
	    	results.fit.1date.PC2=results.1date.PC2
            	if (prop.val.est.1date.PC2>0.5)
  	    	{
			prop.val.est.1date.PC2=1-prop.val.est.1date.PC2
    			temp=mixture.props.1.1date.PC2
    			mixture.props.1.1date.PC2=mixture.props.2.1date.PC2
    			mixture.props.2.1date.PC2=temp
    			mse.score.all.1date.PC2=mse.score.all.1date.PC2[length(mse.score.all.1date.PC2):1]
  	    	}
	    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops.all2))
     	    	ourmix.pop1[match(names(mixture.props.1.1date.PC2),donor.pops.all2)]=mixture.props.1.1date.PC2
     	    	ourmix.pop2[match(names(mixture.props.2.1date.PC2),donor.pops.all2)]=mixture.props.2.1date.PC2
     	    	ourmix.1date.PC2=prop.val.est.1date.PC2*ourmix.pop1+(1-prop.val.est.1date.PC2)*ourmix.pop2
     	    	names(ourmix.1date.PC2)=donor.pops.all2
     	    	ourmix.1date.PC2=ourmix.1date.PC2[ourmix.1date.PC2>prop.cutoff]
     	    	ourmix.1date.PC2=ourmix.1date.PC2/sum(ourmix.1date.PC2)
     	    	donor.pops.1date.PC2=names(ourmix.1date.PC2)
            	prop.val.est.1date.PC2=round(prop.val.est.1date.PC2,2)
            	max.height.1date.PC2=min(mse.score.all.1date.PC2)[1]+sd.fac*sd(mse.score.all.1date.PC2)
           	proportion.est.ci.1date.PC2=c(max(alpha.vec[alpha.vec < alpha.vec[mse.score.all.1date.PC2==min(mse.score.all.1date.PC2)] & mse.score.all.1date.PC2>=max.height.1date.PC2]),min(alpha.vec[alpha.vec > alpha.vec[mse.score.all.1date.PC2==min(mse.score.all.1date.PC2)] & mse.score.all.1date.PC2>=max.height.1date.PC2]))
           	 if (proportion.est.ci.1date.PC2[1]<0) proportion.est.ci.1date.PC2[1]=alpha.vec[1]
            	 prop.val.est.ciII.1date.PC2=paste("(",proportion.est.ci.1date.PC2[1],"-",proportion.est.ci.1date.PC2[2],")",sep='')
	    }

                          ## assuming 2-dates, get mixing info on first event:
            get.fitted.results.2date.PC1=getfitted(weights.mat,ourmix,intercept.2date.mat1,pcnum=1)
 	    intercepts.pca.2date.PC1=get.fitted.results.2date.PC1[2:length(get.fitted.results.2date.PC1)]
  	    fit.quality.2date.PC1=get.fitted.results.2date.PC1[1]
	    if (is.na(fit.quality.2date.PC1)) fit.quality.2date.PC1=0
	    if (fit.quality.2date.PC1==0)
	    {
		prop.val.est.2date.PC1=NA
     		prop.val.est.ciII.2date.PC1=NA
		mse.score.2date.PC1=NA
		mixture.props.1.2date.PC1=NA
		names(mixture.props.1.2date.PC1)=NA
		mixture.props.2.2date.PC1=NA	
		names(mixture.props.2.2date.PC1)=NA
     	    }
	    if (fit.quality.2date.PC1>0)
	    {
		results.2date.PC1=fitwhoadmixed(weights.mat,predmat,pca.frac.val*intercepts.pca.2date.PC1,raw.copyprops.all)
  	    	prop.val.est.2date.PC1=(results.2date.PC1[names(results.2date.PC1)=="alpha"][which(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"])== min(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]),na.rm=TRUE))])$alpha
  	    	mse.score.all.2date.PC1=as.double(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]))
  	    	mse.score.2date.PC1=(results.2date.PC1[names(results.2date.PC1)=="avrelmse"][which(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"])== min(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]),na.rm=TRUE))])$avrelmse
  	    	mu.fitted.2date.PC1=(results.2date.PC1[names(results.2date.PC1)=="mufitted"][which(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"])== min(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]),na.rm=TRUE))])$mufitted
  	    	mixture.props.1.2date.PC1=sort((results.2date.PC1[names(results.2date.PC1)=="pop1fit"][which(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"])== min(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]),na.rm=TRUE))])$pop1fit)
  	    	mixture.props.2.2date.PC1=sort((results.2date.PC1[names(results.2date.PC1)=="pop2fit"][which(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"])== min(unlist(results.2date.PC1[names(results.2date.PC1)=="avrelmse"]),na.rm=TRUE))])$pop2fit)
	    	results.fit.2date.PC1=results.2date.PC1
		if (prop.val.est.2date.PC1>0.5)
  	    	{
			prop.val.est.2date.PC1=1-prop.val.est.2date.PC1
    			temp=mixture.props.1.2date.PC1
    			mixture.props.1.2date.PC1=mixture.props.2.2date.PC1
    			mixture.props.2.2date.PC1=temp
    			mse.score.all.2date.PC1=mse.score.all.2date.PC1[length(mse.score.all.2date.PC1):1]
  	    	}
	    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops.all2))
     	    	ourmix.pop1[match(names(mixture.props.1.2date.PC1),donor.pops.all2)]=mixture.props.1.2date.PC1
     	    	ourmix.pop2[match(names(mixture.props.2.2date.PC1),donor.pops.all2)]=mixture.props.2.2date.PC1
     	    	ourmix.2date.PC1=prop.val.est.2date.PC1*ourmix.pop1+(1-prop.val.est.2date.PC1)*ourmix.pop2
     	    	names(ourmix.2date.PC1)=donor.pops.all2
     	    	ourmix.2date.PC1=ourmix.2date.PC1[ourmix.2date.PC1>prop.cutoff]
     	    	ourmix.2date.PC1=ourmix.2date.PC1/sum(ourmix.2date.PC1)
     	    	donor.pops.2date.PC1=names(ourmix.2date.PC1)
            	prop.val.est.2date.PC1=round(prop.val.est.2date.PC1,2)
            	max.height.2date.PC1=min(mse.score.all.2date.PC1)[1]+sd.fac*sd(mse.score.all.2date.PC1)
            	proportion.est.ci.2date.PC1=c(max(alpha.vec[alpha.vec < alpha.vec[mse.score.all.2date.PC1==min(mse.score.all.2date.PC1)] & mse.score.all.2date.PC1>=max.height.2date.PC1]),min(alpha.vec[alpha.vec > alpha.vec[mse.score.all.2date.PC1==min(mse.score.all.2date.PC1)] & mse.score.all.2date.PC1>=max.height.2date.PC1]))
            	if (proportion.est.ci.2date.PC1[1]<0) proportion.est.ci.2date.PC1[1]=alpha.vec[1]
            	prop.val.est.ciII.2date.PC1=paste("(",proportion.est.ci.2date.PC1[1],"-",proportion.est.ci.2date.PC1[2],")",sep='')
	    }

                          ## assuming 2-dates, get mixing info on second event:
	    if (fit.quality.2date.PC1==1)
	    {
		fit.quality.2date.PC2=0
	    }
	    if (fit.quality.2date.PC1<1)
	    {
		get.fitted.results.2date.PC2=getfitted(weights.mat,ourmix,intercept.2date.mat2,pcnum=1)
 	    	intercepts.pca.2date.PC2=get.fitted.results.2date.PC2[2:length(get.fitted.results.2date.PC2)]
  	    	fit.quality.2date.PC2=get.fitted.results.2date.PC2[1]
	    }
	    if (fit.quality.2date.PC2==0)
	    {
		prop.val.est.2date.PC2=NA
     		prop.val.est.ciII.2date.PC2=NA
		mse.score.2date.PC2=NA
		mixture.props.1.2date.PC2=NA
		names(mixture.props.1.2date.PC2)=NA
		mixture.props.2.2date.PC2=NA	
		names(mixture.props.2.2date.PC2)=NA
     	    }
	    if (fit.quality.2date.PC2>0)
	    {
		results.2date.PC2=fitwhoadmixed(weights.mat,predmat,pca.frac.val*intercepts.pca.2date.PC2,raw.copyprops.all)
  	    	prop.val.est.2date.PC2=(results.2date.PC2[names(results.2date.PC2)=="alpha"][which(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"])== min(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]),na.rm=TRUE))])$alpha
  	    	mse.score.all.2date.PC2=as.double(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]))
  	    	mse.score.2date.PC2=(results.2date.PC2[names(results.2date.PC2)=="avrelmse"][which(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"])== min(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]),na.rm=TRUE))])$avrelmse
  	    	mu.fitted.2date.PC2=(results.2date.PC2[names(results.2date.PC2)=="mufitted"][which(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"])== min(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]),na.rm=TRUE))])$mufitted
  	    	mixture.props.1.2date.PC2=sort((results.2date.PC2[names(results.2date.PC2)=="pop1fit"][which(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"])== min(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]),na.rm=TRUE))])$pop1fit)
  	    	mixture.props.2.2date.PC2=sort((results.2date.PC2[names(results.2date.PC2)=="pop2fit"][which(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"])== min(unlist(results.2date.PC2[names(results.2date.PC2)=="avrelmse"]),na.rm=TRUE))])$pop2fit)
	    	results.fit.2date.PC2=results.2date.PC2
            	if (prop.val.est.2date.PC2>0.5)
  	    	{
			prop.val.est.2date.PC2=1-prop.val.est.2date.PC2
    			temp=mixture.props.1.2date.PC2
    			mixture.props.1.2date.PC2=mixture.props.2.2date.PC2
    			mixture.props.2.2date.PC2=temp
    			mse.score.all.2date.PC2=mse.score.all.2date.PC2[length(mse.score.all.2date.PC2):1]
  	    	}
	    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops.all2))
     	    	ourmix.pop1[match(names(mixture.props.1.2date.PC2),donor.pops.all2)]=mixture.props.1.2date.PC2
     	    	ourmix.pop2[match(names(mixture.props.2.2date.PC2),donor.pops.all2)]=mixture.props.2.2date.PC2
     	    	ourmix.2date.PC2=prop.val.est.2date.PC2*ourmix.pop1+(1-prop.val.est.2date.PC2)*ourmix.pop2
     	    	names(ourmix.2date.PC2)=donor.pops.all2
     	    	ourmix.2date.PC2=ourmix.2date.PC2[ourmix.2date.PC2>prop.cutoff]
     	    	ourmix.2date.PC2=ourmix.2date.PC2/sum(ourmix.2date.PC2)
     	    	donor.pops.2date.PC2=names(ourmix.2date.PC2)
            	prop.val.est.2date.PC2=round(prop.val.est.2date.PC2,2)
            	max.height.2date.PC2=min(mse.score.all.2date.PC2)[1]+sd.fac*sd(mse.score.all.2date.PC2)
            	proportion.est.ci.2date.PC2=c(max(alpha.vec[alpha.vec < alpha.vec[mse.score.all.2date.PC2==min(mse.score.all.2date.PC2)] & mse.score.all.2date.PC2>=max.height.2date.PC2]),min(alpha.vec[alpha.vec > alpha.vec[mse.score.all.2date.PC2==min(mse.score.all.2date.PC2)] & mse.score.all.2date.PC2>=max.height.2date.PC2]))
            	if (proportion.est.ci.2date.PC2[1]<0) proportion.est.ci.2date.PC2[1]=alpha.vec[1]
            	prop.val.est.ciII.2date.PC2=paste("(",proportion.est.ci.2date.PC2[1],"-",proportion.est.ci.2date.PC2[2],")",sep='')
	    }
		# print("Iter5:")
		
		#  print((proc.time() - ptm6)/60)
	}

                       ## (ii) minimise relative mean-squared-error between the fitted intercept vector and the target intercept vector, as well as minimizing proportion fit -- i.e. explaining as much of "ourmix" as possible maybe? (minimize average of the two?):  (!!!! WORKS BEST? !!!!)
  	lik.factor.ci=2         ## 2,3(ii), 5(iii), 10(iv)
  	prop.val.est=(results.fit[names(results.fit)=="alpha"][which(unlist(results.fit[names(results.fit)=="avrelmse"])== min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE))])$alpha
  	range.prop.val.est=as.double(unlist(results.fit[names(results.fit)=="alpha"][which(unlist(results.fit[names(results.fit)=="avrelmse"]) <= (min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE)*lik.factor.ci))]))
  	mse.score.all=as.double(unlist(results.fit[names(results.fit)=="avrelmse"]))
  	mse.score=(results.fit[names(results.fit)=="avrelmse"][which(unlist(results.fit[names(results.fit)=="avrelmse"])== min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE))])$avrelmse
  	mu.fitted=(results.fit[names(results.fit)=="mufitted"][which(unlist(results.fit[names(results.fit)=="avrelmse"])== min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE))])$mufitted
  	mixture.props.1=sort((results.fit[names(results.fit)=="pop1fit"][which(unlist(results.fit[names(results.fit)=="avrelmse"])== min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE))])$pop1fit)
  	mixture.props.2=sort((results.fit[names(results.fit)=="pop2fit"][which(unlist(results.fit[names(results.fit)=="avrelmse"])== min(unlist(results.fit[names(results.fit)=="avrelmse"]),na.rm=TRUE))])$pop2fit)
	if (is.na(prop.val.est)) prop.val.est=0
  	if (prop.val.est>0.5)
  	{
		prop.val.est=1-prop.val.est
    		range.prop.val.est=1-range.prop.val.est
    		temp=mixture.props.1
    		mixture.props.1=mixture.props.2
    		mixture.props.2=temp
    		mse.score.all=mse.score.all[length(mse.score.all):1]
  	}
  	alpha.mat=cbind(alpha.mat,mse.score.all)
      
                       ## make new weights.mat and mixture vec ("ourmix"):
		cutoff = prop.cutoff
		if (a >= num.iterations-1) {cutoff = 0.01}
		mixture.props.1 = mixture.props.1[mixture.props.1>cutoff]
		mixture.props.2 = mixture.props.2[mixture.props.2>cutoff]
		mixture.props.1 = mixture.props.1/sum(mixture.props.1)
		mixture.props.2 = mixture.props.2/sum(mixture.props.2)
		ourmix.pop1=ourmix.pop2 =rep(0,length(donor.pops.all2))
     	ourmix.pop1[match(names(mixture.props.1),donor.pops.all2)]=mixture.props.1
     	ourmix.pop2[match(names(mixture.props.2),donor.pops.all2)]=mixture.props.2
     	ourmix=prop.val.est*ourmix.pop1+(1-prop.val.est)*ourmix.pop2
     	names(ourmix)=donor.pops.all2
     	ourmix=ourmix[ourmix>prop.cutoff]
     	ourmix=ourmix/sum(ourmix)
     	donor.pops=names(ourmix)
	if (length(donor.pops)==1)
	{
		print(paste("A SINGLE DONOR POP FITS MIXTURE AFTER ",a," ITERATIONS OF SOURCE/DATE ESTIMATION! Exiting early....",sep=''))
		write("### INFERRED SOURCES AND DATES ('best-guess' conclusion: Unknown",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT EVIDENCE, DATE ESTIMATE, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.1date","proportion.source1","maxR2fit.1date","fit.quality.1event","fit.quality.2events","bestmatch.event1.source1","bestmatch.event1.source2","proportion.event2.source1","bestmatch.event2.source1","bestmatch.event2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
		write(rep(NA,10),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT EVIDENCE, DATE ESTIMATES, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.2dates.date1","gen.2dates.date2","maxScore.2events","proportion.date1.source1","bestmatch.date1.source1","bestmatch.date1.source2","proportion.date2.source1","bestmatch.date2.source1","bestmatch.date2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
		write(rep(NA,9),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(donor.pops,ourmix),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		write.table(rbind(NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC2:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE1-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE2-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
 		q(save='no')
	}

	weights.mat=matrix(predmat[match(donor.pops,donor.pops.all2),],ncol=dim(predmat)[2])
	rownames(weights.mat)=donor.pops
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]*ourmix
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]/sum(weights.mat[,i])

			## get fitted intercepts:
	test=weights.mat%*%t(predmat)
  	fittedpredmat=cbind(sqrt(prop.val.est)*sqrt(1-prop.val.est)*test,-sqrt(prop.val.est)*sqrt(1-prop.val.est)*test)
	fittedpca.thing=fittedpredmat%*%c(ourmix.pop1,ourmix.pop2)
	scaledintercept.mat.fitted=fittedpca.thing%*%t(fittedpca.thing)
	intercept.mat.fitted=scaledintercept.mat.fitted/(2*ourmix%*%t(ourmix))

			## make weights:
	weights.mat=matrix(predmat[match(donor.pops,donor.pops.all2),],ncol=dim(predmat)[2])
	rownames(weights.mat)=donor.pops
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]*ourmix
	for(i in 1:ncol(weights.mat)) weights.mat[,i]=weights.mat[,i]/sum(weights.mat[,i])
        #weights.mat2=cbind(weights.mat,0)   ## for self-copying (only necessary for curves)
        #rownames(weights.mat2)=donor.pops
        #weights.mat2[is.nan(weights.mat2)]=0.0;
        newdim=nrow(weights.mat)
        olddim=length(donor.pops.all)
  #  ptm3 = proc.time()
	#print('temweight:')
    tempweights = matrix(mutiply.temp.weight(olddim,newdim,weights.mat),byrow=TRUE,ncol=length(donor.pops.all)*length(donor.pops.all))
	
			## generate (weighted) coancestry curves:
	dists.shift=abs(seq(0,grid.range[2],bin.width)-grid.range[1])
	ptm2 = proc.time()
	# print('Coancestry:')
	if (mode == 1 || mode== 2)
	{
	means=matrix(coancestry.curves(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
	}
	if (mode == 3)
	{
	means=matrix(coancestry.curves.mode3(ploidy,ind.id.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
	}
	if (null.ind==1 && null.calc==0)
	{
		null.calc=1
		ptm22 = proc.time()
		# print('Null Curve:')
		tempchromindres.null=matrix(coancestry.curves.null(ploidy,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,simname,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff),byrow=TRUE,ncol=num.bins)
	}
	if (null.ind==1)
	{
		expindres=matrix(0,nrow=nrow(weights.mat)^2,ncol=length(dists))
	        results=tempweights %*% tempchromindres.null
		for (k in 1:length(dists)) results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))

	               #####expectations
		sums=rep(0,length(dists))
		sums2=matrix(0, nrow=nrow(weights.mat),ncol=length(dists))
		sums3=sums2
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat))
	      		{
			    sums=sums+results[(k-1)*nrow(weights.mat)+l,]
			    sums2[k,]=sums2[k,]+results[(k-1)*nrow(weights.mat)+l,]
			    sums3[l,]=sums3[l,]+results[(k-1)*nrow(weights.mat)+l,]		
	      		}
		}
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat)) expindres[(k-1)*nrow(weights.mat)+l,]=sums2[k,]*sums3[l,]/sums
		}
		for (k in 1:length(dists))
  		{
			#results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
			expindres[,k]=c(t((matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
 		}
		results=results/expindres
		means=means/results
	}

                 ## infer new date and intercept.mat:
      	#names3=rownames(weights.mat2)
      	names3=rownames(weights.mat)
		times=dists
		if (a== num.iterations){
		new.grid.range = findpeak(ourmix,means,bin.width)
		#print(new.grid.range)
        new.dists=seq(grid.range[1]+new.grid.range*bin.width,grid.range[2],bin.width)
		means = means[,(dim(means)[2]-length(new.dists)+1):dim(means)[2]]
		times=new.dists
		}
		
	   # }
      	
        pop.comb=paste(rep(names3,each=length(names3)),rep(names3,length(names3)),sep=' ')
        rownames(means)=pop.comb
        check=getmax(means,times,nparam=1,popfracs=ourmix,returnall=T,ourpops=rownames(means))
        admixtimes=check$par
        pred=matrix(nrow=length(times),ncol=1)
        vec=exp(-times*admixtimes[1]/100)
        pred[,1]=vec
        r.squared.est=intercept.est=rep(NA,nrow(means))
        for(i in 1:nrow(means))
 	{
            temp=lm(means[i,]~pred)
            r.squared.est[i]=summary(temp)$r.squared
            intercept.est[i]=round(temp$coeff[2:length(temp$coeff)],10)
        }
        intercept.mat=matrix(intercept.est,ncol=length(donor.pops),byrow=TRUE)
        donor.pops.bestrep=strsplit(pop.comb[intercept.est<0][r.squared.est[intercept.est<0]==max(r.squared.est[intercept.est<0])][1],split=" ")[[1]]

                               ## get final curve info:
        if (a==num.iterations)
        {
	  	    
			ptm7 = proc.time()
                        ## get fit quality scores:
	    fit.quality=getfitted(weights.mat,ourmix,intercept.mat,pcnum=1)[1]
	    fit.quality2=getfitted(weights.mat,ourmix,intercept.mat,pcnum=2)[1]
	    if (is.na(fit.quality2)) fit.quality2=0

                        ## get intercepts based on first pc only:
	    get.fitted.results.new=getfitted(weights.mat,ourmix,intercept.mat,pcnum=1)
		 
  	    intercepts.pca.new=get.fitted.results.new[2:length(get.fitted.results.new)]
	    scaledintercept.mat.firstpc=as.double(intercepts.pca.new)%*%t(as.double(intercepts.pca.new))
	    intercept.mat.firstpc=scaledintercept.mat.firstpc/(2*ourmix%*%t(ourmix))

                        ## get second date and second date intercepts, r^2:
			 ptm8 = proc.time()
			 check2=getmax(means,times,nparam=2,popfracs=ourmix,returnall=T,ourpops=rownames(means))
			admixtimes2=check2$par
            pred2=matrix(nrow=length(times),ncol=2)
	    for(k in 1:2)
            {
                  vec2=exp(-times*admixtimes2[k]/100)
                  pred2[,k]=vec2
            }
                       ## get intercepts for each date:
            r.squared.est.2event=rep(NA,nrow(means))
            intercept.est.2date=matrix(NA,nrow=nrow(means),ncol=2)
            for(i in 1:nrow(means))
            {
                  temp2=lm(means[i,]~pred2)
                  r.squared.est.2event[i]=summary(temp2)$r.squared
                  intercept.est.2date[i,]=round(temp2$coeff[2:length(temp2$coeff)],10)
            }
            intercept.2date.mat1=matrix(intercept.est.2date[,1],ncol=length(donor.pops),byrow=TRUE)
            intercept.2date.mat2=matrix(intercept.est.2date[,2],ncol=length(donor.pops),byrow=TRUE)
	}
	# end of a for loop
    #   print((proc.time() - ptm4)/60)
	  }
	 ptm5 = proc.time()
	 
     r.squared.est.2=r.squared.est
     r.squared.est.2event.2=r.squared.est.2event

			  ## GET "BEST-GUESS" ADMIXTURE CONCLUSION:
     fit.qualityBOTH=fit.quality+fit.quality2
     score.twoevent=max(((r.squared.est.2event-r.squared.est)/(1.0-r.squared.est))[r.squared.est<0.975])
     if (is.na(fit.qualityBOTH))
     {
		print(paste("ERROR IN TAKING EIGEN-DECOMPOSITION AFTER ",a," ITERATIONS OF SOURCE/DATE ESTIMATION! SUGGESTS 'NO ADMIXTURE' SIGNAL. Exiting early....",sep=''))
		write("### INFERRED SOURCES AND DATES ('best-guess' conclusion: Unknown",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT EVIDENCE, DATE ESTIMATE, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.1date","proportion.source1","maxR2fit.1date","fit.quality.1event","fit.quality.2events","bestmatch.event1.source1","bestmatch.event1.source2","proportion.event2.source1","bestmatch.event2.source1","bestmatch.event2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
		write(rep(NA,10),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT EVIDENCE, DATE ESTIMATES, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     		write(c("gen.2dates.date1","gen.2dates.date2","maxScore.2events","proportion.date1.source1","bestmatch.date1.source1","bestmatch.date1.source2","proportion.date2.source1","bestmatch.date2.source1","bestmatch.date2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
		write(rep(NA,9),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(donor.pops,ourmix),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
		write.table(rbind(NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 1-DATE FIT SOURCES, PC2:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE1-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     		write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write("### 2-DATE FIT SOURCES, DATE2-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
		write.table(rbind(NA,NA,NA,NA),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
 		q(save='no')
     }

     if (fit.qualityBOTH<0.985) admix.conclusion="uncertain"
     if (fit.qualityBOTH>=0.985 && score.twoevent>=0.35) admix.conclusion="multiple-dates"
     if (fit.qualityBOTH>=0.985 && score.twoevent<0.35 && fit.quality >= 0.975) admix.conclusion="one-date"
     if (fit.qualityBOTH>=0.985 && score.twoevent<0.35 && fit.quality < 0.975) admix.conclusion="one-date-multiway"
     if (max(r.squared.est.2event)<0.3) admix.conclusion="unclear signal -- check curves/bootstraps"

 			  ## FIND DONOR WITH BEST FITTING COPY-VECTOR TO EACH INFERRED SOURCE:
                                        ## 1-date, PC1:
     copyvec.source1=copyvec.source2=rep(0,dim(predmat)[2])
     for (j in 1:length(names(mixture.props.1))) copyvec.source1=copyvec.source1+mixture.props.1[j]*predmat[rownames(predmat)==names(mixture.props.1)[j],]
     for (j in 1:length(names(mixture.props.2))) copyvec.source2=copyvec.source2+mixture.props.2[j]*predmat[rownames(predmat)==names(mixture.props.2)[j],]
     corr.source1=apply(predmat,1,corr.find,vec.tocheck=copyvec.source1)
     corr.source2=apply(predmat,1,corr.find,vec.tocheck=copyvec.source2)
     best.fit.donors.date1.PC1=c(rownames(predmat)[corr.source1==max(corr.source1)][1],rownames(predmat)[corr.source2==max(corr.source2)][1])
                                        ## 1-date, PC2:
     if (fit.quality.1date.PC2==0) best.fit.donors.date1.PC2=c(NA,NA)
     if (fit.quality.1date.PC2>0)
     {
	copyvec.source1=copyvec.source2=rep(0,dim(predmat)[2])
     	for (j in 1:length(names(mixture.props.1.1date.PC2))) copyvec.source1=copyvec.source1+mixture.props.1.1date.PC2[j]*predmat[rownames(predmat)==names(mixture.props.1.1date.PC2)[j],]
     	for (j in 1:length(names(mixture.props.2.1date.PC2))) copyvec.source2=copyvec.source2+mixture.props.2.1date.PC2[j]*predmat[rownames(predmat)==names(mixture.props.2.1date.PC2)[j],]
     	corr.source1=apply(predmat,1,corr.find,vec.tocheck=copyvec.source1)
     	corr.source2=apply(predmat,1,corr.find,vec.tocheck=copyvec.source2)
     	best.fit.donors.date1.PC2=c(rownames(predmat)[corr.source1==max(corr.source1)][1],rownames(predmat)[corr.source2==max(corr.source2)][1])
     }
                                        ## 2-date, event 1-PC1:
     if (fit.quality.2date.PC1==0) best.fit.donors.date2.PC1=c(NA,NA)
     if (fit.quality.2date.PC1>0)
     {
        copyvec.source1=copyvec.source2=rep(0,dim(predmat)[2])
     	for (j in 1:length(names(mixture.props.1.2date.PC1))) copyvec.source1=copyvec.source1+mixture.props.1.2date.PC1[j]*predmat[rownames(predmat)==names(mixture.props.1.2date.PC1)[j],]
     	for (j in 1:length(names(mixture.props.2.2date.PC1))) copyvec.source2=copyvec.source2+mixture.props.2.2date.PC1[j]*predmat[rownames(predmat)==names(mixture.props.2.2date.PC1)[j],]
     	corr.source1=apply(predmat,1,corr.find,vec.tocheck=copyvec.source1)
     	corr.source2=apply(predmat,1,corr.find,vec.tocheck=copyvec.source2)
     	best.fit.donors.date2.PC1=c(rownames(predmat)[corr.source1==max(corr.source1)][1],rownames(predmat)[corr.source2==max(corr.source2)][1])
     }
                                        ## 2-date, event 2-PC1:
     if (fit.quality.2date.PC2==0) best.fit.donors.date2.PC2=c(NA,NA)
     if (fit.quality.2date.PC2>0)
     {
	copyvec.source1=copyvec.source2=rep(0,dim(predmat)[2])
     	for (j in 1:length(names(mixture.props.1.2date.PC2))) copyvec.source1=copyvec.source1+mixture.props.1.2date.PC2[j]*predmat[rownames(predmat)==names(mixture.props.1.2date.PC2)[j],]
     	for (j in 1:length(names(mixture.props.2.2date.PC2))) copyvec.source2=copyvec.source2+mixture.props.2.2date.PC2[j]*predmat[rownames(predmat)==names(mixture.props.2.2date.PC2)[j],]
     	corr.source1=apply(predmat,1,corr.find,vec.tocheck=copyvec.source1)
     	corr.source2=apply(predmat,1,corr.find,vec.tocheck=copyvec.source2)
     	best.fit.donors.date2.PC2=c(rownames(predmat)[corr.source1==max(corr.source1)][1],rownames(predmat)[corr.source2==max(corr.source2)][1])
     }
                          ## PROPORTION CI:
     prop.val.est=round(prop.val.est,2)
     #max.height=min(alpha.mat[,dim(alpha.mat)[2]])[1]+sd.fac*sd(alpha.mat[,dim(alpha.mat)[2]])
     #proportion.est.ci=c(max(alpha.vec[alpha.vec < alpha.vec[alpha.mat[,dim(alpha.mat)[2]]==min(alpha.mat[,dim(alpha.mat)[2]])] & alpha.mat[,dim(alpha.mat)[2]]>=max.height]),min(alpha.vec[alpha.vec > alpha.vec[alpha.mat[,dim(alpha.mat)[2]]==min(alpha.mat[,dim(alpha.mat)[2]])] & alpha.mat[,dim(alpha.mat)[2]]>=max.height]))
     #if (proportion.est.ci[1]<0) proportion.est.ci[1]=alpha.vec[1]
     #prop.val.est.ciII=paste("(",proportion.est.ci[1],"-",proportion.est.ci[2],")",sep='')

                          ## PRINT OUT PROPORTION, SOURCE, DATE ESTIMATES:
     #save(results.fit,fit.quality,fit.quality2,file=paste(save.file.props,prop.save.addon,prop.save.post,"ALLResults",sep=''))
     write(paste("### INFERRED SOURCES AND DATES ('best-guess' conclusion: ",admix.conclusion,")",sep=''),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 1-DATE FIT EVIDENCE, DATE ESTIMATE, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write(c("gen.1date","proportion.source1","maxR2fit.1date","fit.quality.1event","fit.quality.2events","bestmatch.event1.source1","bestmatch.event1.source2","proportion.event2.source1","bestmatch.event2.source1","bestmatch.event2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
     write(c(admixtimes,prop.val.est,max(r.squared.est),fit.quality,fit.quality+fit.quality2,best.fit.donors.date1.PC1,prop.val.est.1date.PC2,best.fit.donors.date1.PC2),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=10,append=TRUE)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 2-DATE FIT EVIDENCE, DATE ESTIMATES, SINGLE BEST-FITTING DONORS",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write(c("gen.2dates.date1","gen.2dates.date2","maxScore.2events","proportion.date1.source1","bestmatch.date1.source1","bestmatch.date1.source2","proportion.date2.source1","bestmatch.date2.source1","bestmatch.date2.source2"),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
     write(c(admixtimes2,score.twoevent,prop.val.est.2date.PC1,best.fit.donors.date2.PC1,prop.val.est.2date.PC2,best.fit.donors.date2.PC2),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=9,append=TRUE)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 1-DATE FIT SOURCES, PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.1)),c(prop.val.est,mixture.props.1)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.2)),c(1-prop.val.est,mixture.props.2)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 1-DATE FIT SOURCES, PC2:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.1.1date.PC2)),c(prop.val.est.1date.PC2,mixture.props.1.1date.PC2)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.2.1date.PC2)),c(1-prop.val.est.1date.PC2,mixture.props.2.1date.PC2)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 2-DATE FIT SOURCES, DATE1-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.1.2date.PC1)),c(prop.val.est.2date.PC1,mixture.props.1.2date.PC1)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.2.2date.PC1)),c(1-prop.val.est.2date.PC1,mixture.props.2.2date.PC1)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write("#######################################",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write("### 2-DATE FIT SOURCES, DATE2-PC1:",file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),ncolumns=1,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.1.2date.PC2)),c(prop.val.est.2date.PC2,mixture.props.1.2date.PC2)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)
     write.table(rbind(c("proportion",names(mixture.props.2.2date.PC2)),c(1-prop.val.est.2date.PC2,mixture.props.2.2date.PC2)),file=paste(save.file.props,prop.save.addon,prop.save.post,sep=''),row.names=FALSE,col.names=FALSE,quote=FALSE,append=TRUE)

                          ## PRINT OUT CURVE INFORMATION AND PLOT:
     pdf(file=paste(save.file.props,prop.save.addon,".pdf",sep=''),width=16, height=11,paper="a4r")
     par(mfrow=c(numpops.toplot,numpops.toplot))
     write(paste("### INFERRED COANCESTRY CURVES ('best-guess' conclusion: ",admix.conclusion,")",sep=''),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=1)
     write(c("surrogate1","surrogate2","curve.description","rsquared.date.fit","intercept.fit","intercept.fit.date2",times),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=length(times)+6,append=TRUE)
     blue.col=col2rgb("cornflowerblue")
     red.col=col2rgb(2)
     green.col=col2rgb(3)
     print("Writing out and plotting....")
     for (i in 1:length(rownames(means)))
     {
	pop.comb.i=as.character(strsplit(pop.comb[i],split=' ')[[1]])
                                      ## GENERATE FITTED CURVES:
	temp=lm(means[i,]~pred)
	temp2=lm(means[i,]~pred2)
	news=temp$coeff
	news2=temp2$coeff
	predline=news[1]+pred*news[2]
	predline.fitted=news[1]+pred*c(t(intercept.mat.fitted))[i]
	predline.2date=news2[1]+pred2[,1]*news2[2]+pred2[,2]*news2[3]
                                      ## PRINT OUT CURVES:
	write(c(pop.comb.i,"scaled.data",NA,NA,NA,means[i,]),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=length(times)+6,append=TRUE)
	write(c(pop.comb.i,"gen.fit.1date",r.squared.est[i],intercept.est[i],NA,predline),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=length(times)+6,append=TRUE)
	write(c(pop.comb.i,"source.fit.1date",r.squared.est[i],c(t(intercept.mat.fitted))[i],NA,predline.fitted),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=length(times)+6,append=TRUE)
	write(c(pop.comb.i,"gen.fit.2date",r.squared.est.2event[i],c(t(intercept.2date.mat1))[i],c(t(intercept.2date.mat2))[i],predline.2date),file=paste(save.file.props,prop.save.addon,"_curves.txt",sep=''),ncolumns=length(times)+6,append=TRUE)
                                      ## PLOT CURVES:
    	plot.subtitle=simname
	if ((length(xlim.plot)==1) && (xlim.plot[1]<=0)) xlim.plot=c(min(times),max(times))
	plot(times,means[i,],type="l",lwd=3,main=paste(pop.comb.i[1]," vs ",pop.comb.i[2],sep=''),sub=plot.subtitle,xlab="distance (cM)",ylab="",xlim=xlim.plot)
	lines(times,predline.fitted,col=rgb(blue.col[1],blue.col[2],blue.col[3],maxColorValue=255,alpha=155),lwd=4)
	lines(times,predline.2date,col=rgb(red.col[1],red.col[2],red.col[3],maxColorValue=255,alpha=155),lwd=4)
	lines(times,predline,col=rgb(green.col[1],green.col[2],green.col[3],maxColorValue=255,alpha=155),lwd=4)
	abline(h=1,col=1,lty="dotted")
	if (i==1) legend(max(xlim.plot),max(means[i,]),legend=c("data","1-date","1-date (source)","2-date"),lty=rep(1,4),lwd=rep(2,4),col=c(1,3,"cornflowerblue",2),xjust=1,yjust=1,bty='n',cex=1)
     }
     dev.off()
     #save.image("test71.RData")
     #print((proc.time() - ptm5)/60)
}
print((proc.time() - ptm)/60)
##############################################################
#############################################################
### (III) BOOTSTRAPPING DATE ESTIMATES
#############################################################
##############################################################

if (admix.curves.ind==1 && bootstrap.samp>0)
{
	num.bins = length(dists)
    y.intercept.val=matrix(intercept.mat,ncol=1)
    if (nparam.toplot>1)
    {
	admixtimes=admixtimes2
	r.squared.est.2event=r.squared.est.2event.2
	r.squared.est=r.squared.est.2
    }
    write(c("bootstrap.num",paste("date",1:nparam.toplot,".est.boot",sep=''),"maxR2fit.1date.boot","maxScore.2events.boot"),file=paste(save.file.bootstraps,prop.save.addon,prop.save.post,sep=''),ncolumns=nparam.toplot+3)
    boot.dates = matrix(NA,nrow=bootstrap.samp,ncol=nparam.toplot)
    boot.intercepts = array(NA,dim=c(bootstrap.samp,length(pop.comb),nparam.toplot))
    for (b in 1:bootstrap.samp)
    {
	   print(paste("BOOTSTRAP ",b," OF ",bootstrap.samp,"....",sep=''))
	   bootstrap.ind.vec=matrix(sample(1:ninds-1,num.samplefiles*ninds,replace=TRUE),byrow=FALSE,nrow=num.samplefiles)
	   for (i in 1:dim(bootstrap.ind.vec)[1]) bootstrap.ind.vec[i,]=sort(bootstrap.ind.vec[i,])

			## generate (weighted) coancestry curves:
	    dists.shift=abs(seq(0,grid.range[2],bin.width)-grid.range[1])
		
		if (mode == 1 || mode== 2)
		{
		means=matrix(coancestry.curves(ploidy,bootstrap.ind.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
		}
		if (mode == 3)
		{
		means=matrix(coancestry.curves.mode3(ploidy,bootstrap.ind.vec,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
		}
				
		if (null.ind==1 && null.calc==0)
	    {
		null.calc=1
		tempchromindres.null=matrix(coancestry.curves.null(ploidy,samples.filein,recomrates.filein,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,simname,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff),byrow=TRUE,ncol=num.bins)
            }
	    if (null.ind==1)
	    {
		expindres=matrix(0,nrow=nrow(weights.mat)^2,ncol=length(dists))
	        results=tempweights %*% tempchromindres.null
		for (k in 1:length(dists)) results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))

	               #####expectations
		sums=rep(0,length(dists))
		sums2=matrix(0, nrow=nrow(weights.mat),ncol=length(dists))
		sums3=sums2
		#save.image("test.RData")
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat))
	      		{
			    sums=sums+results[(k-1)*nrow(weights.mat)+l,]
			    sums2[k,]=sums2[k,]+results[(k-1)*nrow(weights.mat)+l,]
			    sums3[l,]=sums3[l,]+results[(k-1)*nrow(weights.mat)+l,]		
	      		}
		}
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat)) expindres[(k-1)*nrow(weights.mat)+l,]=sums2[k,]*sums3[l,]/sums
		}
		for (k in 1:length(dists))
  		{
			#results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
			expindres[,k]=c(t((matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
 		}
		results=results/expindres
		means=means/results
	    }

                 ## infer new date and intercept.mat:
            rownames(means)=pop.comb
			times=dists
			
			
			new.grid.range = findpeak(ourmix,means,bin.width)
	        new.dists=seq(grid.range[1]+new.grid.range*bin.width,grid.range[2],bin.width)
			#print(new.grid.range)
			#print(length(new.dists))
			#print(summary(new.dists))
			#print(num.bins)
			#print(dim(means))
			means = means[,(dim(means)[2]-length(new.dists)+1):dim(means)[2]]
			times=new.dists
			
			
	    #print(dim(means))
	    #print(length(times))
	    #print(summary(times))
            check.boot.1=getmax(means,times,nparam=1,popfracs=ourmix,returnall=TRUE,ourpops=rownames(means))
            check.boot.2=getmax(means,times,nparam=2,popfracs=ourmix,returnall=TRUE,ourpops=rownames(means))
            admixtimes.1=check.boot.1$par
            admixtimes.2=check.boot.2$par
            pred.1=matrix(nrow=length(times),ncol=1)
            vec=exp(-times*admixtimes.1[1]/100)
            pred.1[,1]=vec
            pred.2=matrix(nrow=length(times),ncol=2)
	    for(k in 1:2)
            {
                  vec2=exp(-times*admixtimes.2[k]/100)
                  pred.2[,k]=vec2
            }
            rsquared.1event.boot=rsquared.2event.boot=rep(NA,nrow(means))
            for(i in 1:nrow(means))
            {
                  temp1=lm(means[i,]~pred.1)
                  temp2=lm(means[i,]~pred.2)
                  rsquared.1event.boot[i]=summary(temp1)$r.squared
                  rsquared.2event.boot[i]=summary(temp2)$r.squared
            }

            if (nparam.toplot==1) check.boot=check.boot.1
            if (nparam.toplot==2) check.boot=check.boot.2
	    boot.dates[b,]=check.boot$par
	    boot.intercepts[b,,]=check.boot$y.intercept
	    score.twoevent.boot=max(((rsquared.2event.boot-rsquared.1event.boot)/(1.0-rsquared.1event.boot))[rsquared.1event.boot<0.975])
	    write(c(b,boot.dates[b,],max(rsquared.1event.boot),score.twoevent.boot),file=paste(save.file.bootstraps,prop.save.addon,prop.save.post,sep=''),ncolumns=length(boot.dates[b,])+3,append=TRUE)
            #write.table(cbind(curnames,"intercept: ",matrix(boot.intercepts[b,,],ncol=nparam.toplot),rsquared.1event.boot,rsquared.2event.boot),file=paste(save.file.bootstraps,prop.save.addon,prop.save.post,sep=''),append=TRUE,col.names=FALSE,row.names=FALSE,quote=FALSE)
     }
}

if ( admix.curves.ind == 2 )
{
	num.bins = length(dists)
    y.intercept.val=matrix(intercept.mat,ncol=1)
    if (nparam.toplot>1)
    {
	admixtimes=admixtimes2
	r.squared.est.2event=r.squared.est.2event.2
	r.squared.est=r.squared.est.2
    }
    write(c("jackknife.num",paste("date",1:nparam.toplot,".est.boot",sep=''),"maxR2fit.1date.boot","maxScore.2events.boot"),file=paste(save.file.bootstraps,prop.save.addon,prop.save.post,sep=''),ncolumns=nparam.toplot+3)
    boot.dates = matrix(NA,nrow=num.samplefiles,ncol=nparam.toplot)
    boot.intercepts = array(NA,dim=c(num.samplefiles,length(pop.comb),nparam.toplot))
    old_sam <- read.table(samples.filein)
	old_rec <- read.table(recomrates.filein)
	keep_sample_file = ""
	keep_rec_file = ""
	for (b in 1:num.samplefiles)
    {
	   print(paste("JACKKNIFING ",b," OF ",num.samplefiles,"....",sep=''))
	   temp.ind.vec = 1:(num.samplefiles)
	   temp.ind.vec = temp.ind.vec[-which(temp.ind.vec == b)]
	   jackknife.ind.vec=matrix(rep(1:ninds-1,each=num.samplefiles-1),byrow=FALSE,nrow=num.samplefiles-1)
	 
       jackknife_sample = paste("jackknife_sample_",toString(b),paste0(format(Sys.time(), "%Y%m%d_%H%M%S_")), ".txt", sep = "")
      	file_sample = jackknife_sample
       file_sample = sink(file_sample)
	   for (j in 1:length(temp.ind.vec))
       {
		cat(paste(old_sam[temp.ind.vec[j],], sep = ""),file_sample)
		cat('\n')
     	}
    	sink()
       jackknife_rec = paste("jackknife_recom_",toString(b),paste0(format(Sys.time(), "%Y%m%d_%H%M%S_")), ".txt", sep = "")
    	file_rec = jackknife_rec
       file_rec = sink(file_rec)
    	for (j in 1:length(temp.ind.vec))
       {
		cat(paste(old_rec[temp.ind.vec[j],], sep = ""),file_sample)
		cat('\n')
    	}
     	sink()
	 
			## generate (weighted) coancestry curves:
	    dists.shift=abs(seq(0,grid.range[2],bin.width)-grid.range[1])
		
		if (mode == 1 || mode== 2)
		{
		means=matrix(coancestry.curves(ploidy,jackknife.ind.vec,jackknife_sample,jackknife_rec,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
		}
		if (mode == 3)
		{
		means=matrix(coancestry.curves.mode3(ploidy,jackknife.ind.vec,jackknife_sample,jackknife_rec,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff,tempweights,nrow(weights.mat)),byrow=TRUE,ncol=num.bins)
		}
				
		if (null.ind==1 && null.calc==0)
	    {
		null.calc=1
		tempchromindres.null=matrix(coancestry.curves.null(ploidy,jackknife_sample,jackknife_rec,bin.width,num.bins,(1:length(dists.shift))[dists.shift==min(dists.shift)][1]-1,simname,donor.pops.all,donor.label.vec,recipient.label.vec,weightmin,gridweightMAXcutoff),byrow=TRUE,ncol=num.bins)
            }
	    if (null.ind==1)
	    {
		expindres=matrix(0,nrow=nrow(weights.mat)^2,ncol=length(dists))
	        results=tempweights %*% tempchromindres.null
		for (k in 1:length(dists)) results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))

	               #####expectations
		sums=rep(0,length(dists))
		sums2=matrix(0, nrow=nrow(weights.mat),ncol=length(dists))
		sums3=sums2
		#save.image("test.RData")
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat))
	      		{
			    sums=sums+results[(k-1)*nrow(weights.mat)+l,]
			    sums2[k,]=sums2[k,]+results[(k-1)*nrow(weights.mat)+l,]
			    sums3[l,]=sums3[l,]+results[(k-1)*nrow(weights.mat)+l,]		
	      		}
		}
		for(k in 1:nrow(weights.mat))
		{
			for(l in 1:nrow(weights.mat)) expindres[(k-1)*nrow(weights.mat)+l,]=sums2[k,]*sums3[l,]/sums
		}
		for (k in 1:length(dists))
  		{
			#results[,k]=c(t((matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(results[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
			expindres[,k]=c(t((matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))+t(matrix(expindres[,k],byrow=TRUE,ncol=nrow(weights.mat))))/2.0))
 		}
		results=results/expindres
		means=means/results
	    }

                 ## infer new date and intercept.mat:
            rownames(means)=pop.comb
			times=dists
			
			
			new.grid.range = findpeak(ourmix,means,bin.width)
	        new.dists=seq(grid.range[1]+new.grid.range*bin.width,grid.range[2],bin.width)
			#print(new.grid.range)
			means = means[,(dim(means)[2]-length(new.dists)+1):dim(means)[2]]
			times=new.dists
			
			
            check.boot.1=getmax(means,times,nparam=1,popfracs=ourmix,returnall=TRUE,ourpops=rownames(means))
            check.boot.2=getmax(means,times,nparam=2,popfracs=ourmix,returnall=TRUE,ourpops=rownames(means))
            admixtimes.1=check.boot.1$par
            admixtimes.2=check.boot.2$par
            pred.1=matrix(nrow=length(times),ncol=1)
            vec=exp(-times*admixtimes.1[1]/100)
            pred.1[,1]=vec
            pred.2=matrix(nrow=length(times),ncol=2)
	    for(k in 1:2)
            {
                  vec2=exp(-times*admixtimes.2[k]/100)
                  pred.2[,k]=vec2
            }
            rsquared.1event.boot=rsquared.2event.boot=rep(NA,nrow(means))
            for(i in 1:nrow(means))
            {
                  temp1=lm(means[i,]~pred.1)
                  temp2=lm(means[i,]~pred.2)
                  rsquared.1event.boot[i]=summary(temp1)$r.squared
                  rsquared.2event.boot[i]=summary(temp2)$r.squared
            }

            if (nparam.toplot==1) check.boot=check.boot.1
            if (nparam.toplot==2) check.boot=check.boot.2
	    boot.dates[b,]=check.boot$par
	    boot.intercepts[b,,]=check.boot$y.intercept
	    score.twoevent.boot=max(((rsquared.2event.boot-rsquared.1event.boot)/(1.0-rsquared.1event.boot))[rsquared.1event.boot<0.975])
	    write(c(b,boot.dates[b,],max(rsquared.1event.boot),score.twoevent.boot),file=paste(save.file.bootstraps,prop.save.addon,prop.save.post,sep=''),ncolumns=length(boot.dates[b,])+3,append=TRUE)
        file.remove(jackknife_rec,jackknife_sample)
		} 
		
		
}


warnings()
print("Finished!")

q(save = 'no')
