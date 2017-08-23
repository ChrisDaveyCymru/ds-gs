#R script for genomic selection.


#Areas to up-date for a given run are high-lighted by #-------- around comments.


#------------------------------
#Input file formats.
#There are 3 input files:
#(1) Phenotype data. File has column headings. Separate row for each miscanthus genotype. The miscanthus genotype column
#    must have name geno. The columns are broken into 2 sets of contiguous blocks. The first set
#    are information columns, e.g. geno and things like plate_id, order, longitude etc. The second
#    set are the trait phenotype data as BLUPs with the same trait measured in different years having
#    separate columns. The miscanthus genotypes must be in the same order as in the SNP genotype file.
#    The traits must be in the same order as in the H2 file. NA used to mark missing BLUP trait data.
#(2) SNP genotype file. No column headings. One row for each SNP. One column for each miscanthus 
#    genotype. SNP genotypes in -1/0/1 (NA) format. 
#(3) H2 file. Two columns with headings. The H2 column must be called H2. The other column contains the
#    trait names. 
#
#Output file formats.
#There are 4 output files.
#(1) Predictive abilities file. The row and column names give the trait names. There is a row for each trait.
#    There are two blocks of columns. The first block are the mean predictive abilities for each of the traits.
#    The second block are the SDs of the predictive abilities.
#(2) The accuracies file has exactly the same format as the predictive abilities file except the mean
#    and SD of the accuracies are recorded.
#(3) Regression coefficients file.
#(4) Workspace file for R.
#------------------------------

#------------------------------
#Enter details of input files.
#Enter phenotype file name.
pheno_input_name <- "si8_gs_si_blups.csv"

#Number of information columns at start of phenotype file.
blank_cols = 7 
#Number of trait phenotype columns in phenotype file entered as 1 to number of traits.
traits<-1:42 # Could also include latitude, longitude, PC1-10

#Enter SNP genotype file name
snp_input_name <- "flora_Sb_liberal_rrblup"

#Enter H2 file name.
h2_input_name <- "si8_gs_si_h2.csv"
#------------------------------

#------------------------------
#Configure the run.
#Number reps. Used so can get useful means and SDs
#of the predictive abilities and accuracies.
reps <-100

#For a given rep and trait the number of folds controls
#the number of miscanthus genotypes trait phenotypes (BLUPs) that are predicted from 
#the other miscanthus genotypes data in one go. When the fold is
#10 and there are 138 miscanthus genotypes then 10 blocks of
#phenotype predictions are made. Each block (apart from the last)
#predicts 13 (138/10) genotypes from the others genotypes phenotype (BLUP trait) data.
#Each miscanthus genotypes phenotype data is predicted only once for a given trait on 
#a given rep.
folds<-10
#------------------------------

#------------------------------
#Enter details of output files.
#Enter name of predictive abilities file.
predict_output_name <- "si8_gs_si_out_pred.csv"

#Enter name of accuracies file.
accuracy_output_name <- "si8_gs_si_out_acc.csv"

#Enter name of regression coefficients file.
reg_coeff_output_name <- "si8_gs_si_out_coeff.csv"

#Enter name of R workspace file.
work_sp_output_name <- "si8_gs_si_out_ws.RData"
#------------------------------

#------------------------------
#The above information entry is all that is needed to run this script.
#------------------------------

#On ecb will need to install this library for each users account.
library(rrBLUP)



#Load SNP genos
G<-as.matrix(read.table(snp_input_name ,header=F))
G<-t(G)
print(dim(G))

#Load phenotype and H2 data
pheno_data<-data.frame(read.csv(pheno_input_name, header=T))
H2 <- data.frame(read.csv(h2_input_name, header=T))


#Calc kinship matrix using rrBLUP. Missing data are imputed.
K <- A.mat(G)
rownames(K)<-pheno_data$geno

#Create storage 3D arrays for results across reps.
r_yy<- array(NA,c(length(traits),length(traits),reps)) # predictive abilities
r_gg<- array(NA,c(length(traits),length(traits),reps)) # accuracies
coeffs <- array(NA,c(length(traits),2,reps)) # slopes and intercepts of regression 




#Loop round the reps
for (repno in 1:reps) {

	print(paste("rep #",repno))
	
	#For this rep calc the number of miscanthus genotypes predicted in fold.
	pred_n<-floor(nrow(G)/folds)

	blup_pred <- as.matrix(rep(NA,nrow(G)))

	#For this rep the random order the miscanthus genotypes predicted in.
	perm_order<-sample(1:nrow(G))


	#Loop through each trait and predict.
	for (trait in traits) {
				# Pick a trait
				P<-data.frame(pheno_data[,blank_cols+trait])
				colnames(P)<-"pheno"
				
				traitname <- colnames(pheno_data)[blank_cols+trait]
				print(traitname)
				
						
						#Loop round the folds and do the actual predictions.
						for (i in 1:folds) {
						
							#Select the miscanthus genotypes to predict on this fold using
							#the other genotypes.
							if (i<folds){ 
								pred_idx<-perm_order[(pred_n*i-(pred_n-1)):(pred_n*i)] 
							}else{ 
								pred_idx<-perm_order[(pred_n*i-(pred_n-1)):(nrow(P))] 
							}#End of if
							
							
							#Set-up data to make predictions.
							train_idx <- setdiff((1:nrow(P)),pred_idx)
							P_train<-P$pheno
							P_train[pred_idx]<-NA
							d<-data.frame(pheno=P_train,id=pheno_data$geno)
							
							#rrBLUP function for genotypic prediction based on kinship
							out<-kin.blup(d,pheno="pheno",geno="id",K=K) 
							
							#Cumulative storage of predicted trait values for miscanthus genotypes.
							blup_pred[pred_idx] = out$g[pred_idx]
						}#End of fold loop.
						
				
				#Correlate the predicted values with the actual values and store 
				#resulting predictive abilities in 3D array at this rep.
				r_yy[trait,,repno] <-  cor(blup_pred,pheno_data[,(blank_cols+1):(blank_cols+length(traits))],use="pairwise.complete.obs")
				
				#Calculate accuracies and store in 3D array at this rep.
				r_gg[trait,,repno] <-  r_yy[trait,,repno]/sqrt(H2$H2[trait])
				
				#Store slopes and intercepts of regression in 3D array at this rep
				coeffs[trait,1,repno] <-  coef(lm(P$pheno~blup_pred))[1] # Regression intercepts
				coeffs[trait,2,repno] <-  coef(lm(P$pheno~blup_pred))[2] # Regression slopes
	}#End of trait loop.
}#End of rep loop.



#For each trait-by-trait location find the mean predictive ability and the SD of
#the predictive abilities across all the reps and store ready for output.
r_yy_bar<-apply(r_yy,c(1,2),mean) # Ave r_yy
r_yy_sd<-apply(r_yy,c(1,2),sd)	# SD r_yy
r_yy_summary<-(cbind(r_yy_bar, r_yy_sd))
rownames(r_yy_summary)<-colnames(pheno_data)[(blank_cols+1):(blank_cols+length(traits))]
colnames(r_yy_summary)<-rep(rownames(r_yy_summary),2)


#For each trait-by-trait location find the mean accuracy and the SD of
#the accuracies across all the reps and store ready for output. 
r_gg_bar<-apply(r_gg,c(1,2),mean) # Ave r_gg
r_gg_sd<-apply(r_gg,c(1,2),sd)	# SD r_gg
r_gg_summary<-(cbind(r_gg_bar, r_gg_sd))
rownames(r_gg_summary)<-colnames(pheno_data)[(blank_cols+1):(blank_cols+length(traits))]
colnames(r_gg_summary)<-rep(rownames(r_gg_summary),2)


#Means and SDs of the regression coefficients ready for output.
coeffs_bar<-apply(coeffs,c(1,2),mean) # Ave coeffs
coeffs_sd<-apply(coeffs,c(1,2),sd)	# SD coeffs
coeffs_summary<-(cbind(coeffs_bar, coeffs_sd))
rownames(coeffs_summary)<-colnames(pheno_data)[(blank_cols+1):(blank_cols+length(traits))]
colnames(coeffs_summary)<-c('bo','b1','sd_bo','sd_b1')


# Print/save files
write.csv(r_yy_summary, predict_output_name)
write.csv(r_gg_summary, accuracy_output_name)
write.csv(coeffs_summary, reg_coeff_output_name)
rm(G)

print("***********Program run complete***********")
save.image(work_sp_output_name)

