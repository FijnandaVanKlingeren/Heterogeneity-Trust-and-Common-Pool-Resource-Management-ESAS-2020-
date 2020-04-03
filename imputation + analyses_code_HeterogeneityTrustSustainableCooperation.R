rm(list=ls())
options(scipen=999)
library(foreign)
library(readstata13)
library(ranger)
library(foreach)
library(missRanger)
library(ltm)
library(foreach)
library(mice)
library(MASS)
library(arm)
library(lmtest)
library(DescTools)
library(pastecs)
library(Hmisc)
library(MASS)
n.sims = 100

#-----------------------------------------------------------------------------------------------------------------------------
# Prepare data for Multiple Imputation procedure
#-----------------------------------------------------------------------------------------------------------------------------

dt = read.table(file = 'C://Users/Fijnanda/Dropbox/Nuffield DPhil/CPR Database/Dataset/CPR_12012018.txt',sep = '\t',header=TRUE)
head(dt)
# remove CL variables - these indicate confidence level and ar enot directly related to the variables of interest
dt_temp = dt
# remove unstandardized survey answers, and other variables which ar convolutedly coded - keep lowercase so you can look in stata for variable
dt_temp = dt_temp[,-grep('oplevelrealyes',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('subgroupsbgpdes',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('oplevelsg',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('subgroupsethid2',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('subgroupsclanid2',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('subgroupscaste2',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('subgroupsgender2',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('screenerslocsize',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('screenerscreenid',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('screenernumber',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('screenertype',gsub("\\.","",tolower(names(dt_temp))))]
dt_temp = dt_temp[,-grep('oplevelconditon',gsub("\\.","",tolower(names(dt_temp))))]

names(dt_temp)=gsub("\\.","",tolower(names(dt_temp)))

# getting rid of variables which do not inform the predictive model of any of the others - useless admin vars
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelbegdate')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelenddate')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelestndate')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelsubgp1')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelsubgp2')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelsubgp3')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelsubgp4')] 
dt_temp=dt_temp[,-which(names(dt_temp)=='subgroupsubgrpid')]
dt_temp=dt_temp[,-which(names(dt_temp)=="resourceresourid")]
dt_temp=dt_temp[,-which(names(dt_temp)=="resourcefk_screenid")]
#dt_temp=dt_temp[,-which(names(dt_temp)=="opleveloplevid")]
dt_temp=dt_temp[,-which(names(dt_temp)=="oplevelfk_resourid")]
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelduration')]
dt_temp=dt_temp[,-which(names(dt_temp)=='oplevelinternal')]

head(dt_temp)
# correct typo 
dt_temp$screenerregion[grep("5rient",dt_temp$screenerregion)]="Orient"
# remove .number - it means confidence level
# if in last two strings is a dot, remove the last two strings
#oplevelid actually has meaningful post-dot numbers
temp = dt_temp$opleveloplevid
dt_temp_clean = foreach(i = 1:dim(dt_temp)[2],.combine = 'cbind') %do% {
if(sum(grepl("\\.",substr(as.character(unlist(dt_temp[,i])),nchar(as.character(unlist(dt_temp[,i])))-1,nchar(as.character(unlist(dt_temp[,i]))))))>0){
 substr(as.character(unlist(dt_temp[,i])),1,nchar(as.character(unlist(dt_temp[,i])))-2)
}else{as.character(unlist(dt_temp[,i]))} }
dt_temp_clean = as.data.frame(dt_temp_clean)
names(dt_temp_clean)=names(dt_temp)
dt_temp_clean$opleveloplevid = temp

# and now make zeroes, -1, -2 into NA
dt_temp_clean = apply(dt_temp_clean,2,function(x){as.character(unlist(ifelse(as.character(unlist(x))==""|
                                                               as.character(unlist(x))=="\\."|
                                                               as.character(unlist(x))=="-"|
                                                               as.character(unlist(x))=="-1"|
                                                               as.character(unlist(x))=="-2"|
                                                               as.character(unlist(x))=="0",NA,as.character(unlist(x))
                                                             )))
                                                               })
dt_temp_clean = as.data.frame(dt_temp_clean)

dt_temp_clean[,which(names(dt_temp_clean)=='resourcesurfarea')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='resourcesurfarea')] )))
dt_temp_clean[,which(names(dt_temp_clean)=='resourcelength')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='resourcelength')])))
dt_temp_clean[,which(names(dt_temp_clean)=='resourcestorevol')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='resourcestorevol')])))
dt_temp_clean[,which(names(dt_temp_clean)=='resourceflowvol')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='resourceflowvol')]))) 
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate1')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate1')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate1')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate1')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate2')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate2')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate2')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate2')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate3')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelbegrate3')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate3')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelendrate3')])))
dt_temp_clean[,which(names(dt_temp_clean)=='oplevelnumsubgp')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='oplevelnumsubgp')])))
dt_temp_clean[,which(names(dt_temp_clean)=='subgroupbnumusr1')]  = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='subgroupbnumusr1')] )))
dt_temp_clean[,which(names(dt_temp_clean)=='subgroupenumusr1')] = as.numeric(as.character(unlist(dt_temp_clean[,which(names(dt_temp_clean)=='subgroupenumusr1')])))
# drop names because too wildly coded and too any categories
names_temp = dt_temp_clean[,grep("name",names(dt_temp_clean))]
dt_temp_clean = dt_temp_clean[,-grep("name",names(dt_temp_clean))]

# plot histogram of % missing for each variable in raw data
hist(sapply(1:dim(dt_temp_clean)[2],function(x){sum(is.na(dt_temp_clean[,x]))})/dim(dt_temp_clean)[1],
     main = "% missing raw data", xlab = "% NA",breaks = 100)

# if % missing < 75, keep, otherwise remove
dt_temp_clean_MI = dt_temp_clean[,which( (sapply(1:dim(dt_temp_clean)[2],function(x){ sum(is.na(dt_temp_clean[,x]))})/dim(dt_temp_clean)[1])<0.75)]

# calculate % missing per remaining variables
pct_miss = apply(dt_temp_clean_MI,2,function(x){sum(is.na(x))/length(x)})



table(dt_temp_clean_MI$oplevelendpoll)
table(dt_temp_clean_MI$subgroupkpresure)
table(dt_temp_clean_MI$subgroupfamincde)
table(dt_temp_clean_MI$resourcevarspace)
table(dt_temp_clean_MI$oplevelworstoff)
table(dt_temp_clean_MI$oplevelsocsanct)
table(dt_temp_clean_MI$subgroupsubnot)


# MI procedure, setting k = 5 as this is the default on mice. We'll perform sensitivity analysis for appendix
n.sims = 100
dt_MI = foreach(i = 1:n.sims) %do% missRanger(dt_temp_clean_MI,maxiter = 30,verbose = 2,returnOOB = TRUE,pmm.k = 5)



#-----------------------------------------------------------------------------------------------------------------------------
# Save data
#-----------------------------------------------------------------------------------------------------------------------------

#Save available case data
setwd("G:/Nuffield DPhil/Paper 1/Data")
save(dt_temp_clean_MI, file = "dt_temp_clean_MI.RData", compress = TRUE)

#Save imputed data
setwd("G:/Nuffield DPhil/Paper 1/Data")
save(dt_MI,file = "dt_MI.RData",compress=TRUE)

#Load imputed and available data to create variables 
rm(list=ls())
setwd("D:/Nuffield DPhil/Paper 1/Data")
load(file = "dt_MI.RData")
load(file = "dt_temp_clean_MI.RData")
n.sims = 100


#-----------------------------------------------------------------------------------------------------------------------------
# Creating variables
#-----------------------------------------------------------------------------------------------------------------------------

# change name of oplevelid to slice ID 
for(i in 1:n.sims){
names(dt_MI[[i]])[which(names(dt_MI[[i]])=="opleveloplevid")]="sliceid"
}
#creating balance
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelendblnc")]="balance_factor"
}

for(i in 1:n.sims){
  levels(dt_MI[[i]]$balance_factor)=c("0", "1", "2", "3")
}

table(dt_MI[[i]]$balance_factor, dt_MI[[i]]$sliceid) #check that balance does not have multiple different values over different timeslices


#creating unit quality 
for(i in 1:n.sims){
  dt_MI[[i]]$unit_quality_cont=5-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelendqual")])))
}

table(dt_MI[[i]]$unit_quality_cont, dt_MI[[i]]$sliceid)

#create trust 
for(i in 1:n.sims){
  dt_MI[[i]]$trust=3-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelendtrust")])))
}

table(dt_MI[[i]]$trust, dt_MI[[i]]$sliceid)


# create monitoring_bribery
for(i in 1:n.sims){
  dt_MI[[i]]$monitoring_bribery=6-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelbribery")])))
}
table(dt_MI[[i]]$monitoring_bribery, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,monitoring_bribery:=mean(as.numeric(as.character(unlist(monitoring_bribery)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}


table(dt_temp_clean_MI$oplevelbribery, dt_temp_clean_MI$opleveloplevid)


# create formal sanctions
for(i in 1:n.sims){
  dt_MI[[i]]$formal_sanctions=6-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelmonsanct")])))
}
table(dt_MI[[i]]$formal_sanctions, dt_MI[[i]]$sliceid)


# create pressure
for(i in 1:n.sims){
  dt_MI[[i]]$pressure=2-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="subgroupkpresure")])))
}
table(dt_MI[[i]]$pressure, dt_MI[[i]]$sliceid)

table(dt_temp_clean_MI$subgroupkpresure, dt_temp_clean_MI$opleveloplevid)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,pressure:=mean(as.numeric(as.character(unlist(pressure)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}


# create income dependent
for(i in 1:n.sims){
  dt_MI[[i]]$incdependent=4-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="subgroupfamincde")])))
}
table(dt_MI[[i]]$incdependent, dt_MI[[i]]$sliceid)

for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,incdependent:=mean(as.numeric(as.character(unlist(incdependent)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}

#create rulefollowing
for(i in 1:n.sims){
  dt_MI[[i]]$rulefollowing=6-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="subgrouprulefoll")])))
}
table(dt_MI[[i]]$rulefollowing, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,rulefollowing:=mean(as.numeric(as.character(unlist(rulefollowing)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}

# create social sanctions
for(i in 1:n.sims){
  dt_MI[[i]]$socialsanct=6-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelsocsanct")])))
}
table(dt_MI[[i]]$socialsanct, dt_MI[[i]]$sliceid)

table(dt_temp_clean_MI$oplevelsocsanct, dt_temp_clean_MI$opleveloplevid)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,socialsanct:=mean(as.numeric(as.character(unlist(socialsanct)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}

# create physical sanctions
for(i in 1:n.sims){
  dt_MI[[i]]$physicsanct=6-as.numeric(as.character(unlist(dt_MI[[i]][,which(names(dt_MI[[i]])=="oplevelphysanct")])))
}
table(dt_MI[[i]]$physicsanct, dt_MI[[i]]$sliceid)

table(dt_temp_clean_MI$oplevelphysanct, dt_temp_clean_MI$opleveloplevid)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,physicsanct:=mean(as.numeric(as.character(unlist(physicsanct)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}

#creating variance in income
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="subgroupsubvar")]="varinc"
}
table(dt_MI[[i]]$varinc, dt_MI[[i]]$sliceid) #varinc will later be used to make ecohet, where the max and mean will be taken 

#creating proportion exit options
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="subgroupsubnot")]="exit_prop"
}

table(dt_MI[[i]]$exit_prop, dt_MI[[i]]$sliceid)

for(i in 1:n.sims){
  dt_MI[[i]]$exit_prop = as.numeric(as.character(unlist(dt_MI[[i]]$exit_prop)))
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,exit_prop:=mean(exit_prop),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}


#creating number of users
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="subgroupenumusr1")]="numusers"
}
#creating estimated number of users
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="subgroupenumusr2")]="est_numusers"
}
#create closed access 
for(i in 1:n.sims){
levels(dt_MI[[i]]$opleveleapclose)=c("1", "7", "5", "4", "3" ,"2")
names(dt_MI[[i]])[which(names(dt_MI[[i]])=="opleveleapclose")]="closedaccess"
}
table(dt_MI[[i]]$closedaccess, dt_MI[[i]]$sliceid)


#create variation in flow over space 
for(i in 1:n.sims){
  levels(dt_MI[[i]]$resourcevarspace) = c("0", "1")
}

table(dt_MI[[i]]$resourcevarspace, dt_MI[[i]]$sliceid)

#create pollution
for(i in 1:n.sims){
  levels(dt_MI[[i]]$oplevelendpoll)=c("0", "1")
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelendpoll")]="pollution"
}
table(dt_MI[[i]]$pollution, dt_MI[[i]]$sliceid)



#create race id, religious id, ethnicity id, clan id, common language and sex id 
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelraceid")]="raceid"
}
table(dt_MI[[i]]$raceid, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelrelanims")]="relanims"
}
table(dt_MI[[i]]$relanims, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelethncid")]="ethncid"
}
table(dt_MI[[i]]$ethncid, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelclanid")]="clanid"
}
table(dt_MI[[i]]$clanid, dt_MI[[i]]$sliceid)
table(dt_temp_clean_MI$oplevelclanid, dt_temp_clean_MI$opleveloplevid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelcommlang")]="commlang"
}
table(dt_MI[[i]]$commlang, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelsex")]="sex"
}
table(dt_MI[[i]]$sex, dt_MI[[i]]$sliceid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelsocstrat")]="casteid"
}
table(dt_MI[[i]]$casteid, dt_MI[[i]]$sliceid)
table(dt_temp_clean_MI$oplevelsocstrat, dt_temp_clean_MI$opleveloplevid)
for(i in 1:n.sims){
  names(dt_MI[[i]])[which(names(dt_MI[[i]])=="oplevelcultvwr")]="cultvwr"
}
table(dt_MI[[i]]$cultvwr, dt_MI[[i]]$sliceid)


#Create sochet as latent variable with factor analysis
library(psych)
library(plyr)
sochet_alpha_drop=data.frame()
sochet_alpha = data.frame()
for(i in 1:n.sims){
dt_MI[[i]][c("raceid", "relanims", "ethncid", "clanid", "commlang", "sex")] = apply(dt_MI[[i]][c("raceid", "relanims", "ethncid", "clanid", "commlang", "sex")],2,function(x){as.numeric(as.character(unlist(x)))})
sochet_alpha = rbind(sochet_alpha, alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","clanid","commlang","sex")])$total)

drop_dt = as.data.frame(t(data.frame(drop_alpha = alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","clanid","commlang","sex")])$alpha.drop$raw_alpha)))
names(drop_dt) = row.names(alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","clanid","commlang","sex")])$alpha.drop)
sochet_alpha_drop = plyr::rbind.fill(sochet_alpha_drop,drop_dt)

dt_MI[[i]]$sochet_latent=rowMeans(dt_MI[[i]][c("raceid", "relanims", "ethncid", "clanid", "commlang", "sex")] )
}
# check which, after averaging over all simulations, is worth removing to make the scale better 
colMeans(sochet_alpha_drop)
# remove clanid
sochet_alpha_drop=data.frame()
sochet_alpha = data.frame()
for(i in 1:n.sims){
  dt_MI[[i]][c("raceid", "relanims", "ethncid", "commlang", "sex")] = apply(dt_MI[[i]][c("raceid", "relanims", "ethncid", "commlang", "sex")],2,function(x){as.numeric(as.character(unlist(x)))})
  sochet_alpha = rbind(sochet_alpha, alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","commlang","sex")])$total)
  
  drop_dt = as.data.frame(t(data.frame(drop_alpha = alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","commlang","sex")])$alpha.drop$raw_alpha)))
  names(drop_dt) = row.names(alpha(dt_MI[[i]][,c("raceid","relanims","ethncid","commlang","sex")])$alpha.drop)
  sochet_alpha_drop = plyr::rbind.fill(sochet_alpha_drop,drop_dt)
  
  dt_MI[[i]]$sochet_latent=rowMeans(dt_MI[[i]][c("raceid", "relanims", "ethncid", "commlang", "sex")] )
}


#Create sochet as maximum found in any of the seven categories
library(data.table)
for(i in 1:n.sims){
  dt_MI[[i]]$raceid_cont = as.numeric(as.character(unlist(dt_MI[[i]]$raceid)))
  dt_MI[[i]]$relanims_cont = as.numeric(as.character(unlist(dt_MI[[i]]$relanims)))
  dt_MI[[i]]$ethncid_cont = as.numeric(as.character(unlist(dt_MI[[i]]$ethncid)))
  dt_MI[[i]]$commlang_cont = as.numeric(as.character(unlist(dt_MI[[i]]$commlang)))
  dt_MI[[i]]$sex_cont = as.numeric(as.character(unlist(dt_MI[[i]]$sex)))
  dt_MI[[i]]$clanid_cont = as.numeric(as.character(unlist(dt_MI[[i]]$clanid)))
  dt_MI[[i]]$casteid_cont = as.numeric(as.character(unlist(dt_MI[[i]]$casteid)))
  dt_MI[[i]]$cultvwr_cont = as.numeric(as.character(unlist(dt_MI[[i]]$cultvwr)))
}

for(i in 1:n.sims){
  dt_MI[[i]]$sochet_max_temp = apply(dt_MI[[i]][,c("raceid_cont","relanims_cont","ethncid_cont","commlang_cont","sex_cont","clanid_cont", "casteid_cont")],1,max)
}

# Maximum value sochet BY sliceid
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,sochet_max:=max(as.numeric(as.character(unlist(sochet_max_temp)))),by=c("sliceid")]
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}

table(dt_MI[[i]]$sochet_max_temp, dt_MI[[i]]$sliceid)
table(dt_MI[[i]]$sochet_max, dt_MI[[i]]$sliceid)

#center to min = 0 
for(i in 1:n.sims){
  dt_MI[[i]]$sochet_max=(dt_MI[[i]][,which(names(dt_MI[[i]])=="sochet_max")])-1
}

dt_MI[[i]]$sochet_max


#create sochet as mean of all seven categories 
for(i in 1:n.sims){
  dt_MI[[i]]$sochet_mean_temp = apply(dt_MI[[i]][,c("raceid_cont","relanims_cont","ethncid_cont","commlang_cont","sex_cont","clanid_cont", "casteid_cont")],1,mean)
}
# sochet mean by sliceid
library(data.table)
for(i in 1:n.sims){
dt_MI[[i]] = as.data.table(dt_MI[[i]])[,sochet_mean:=mean(as.numeric(as.character(unlist(sochet_mean_temp)))),by=c("sliceid")]
dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}


table(dt_MI[[i]]$sochet_mean, dt_MI[[i]]$sliceid)

# max variance of average annual family income across families by sliceid = ecohet_max
library(data.table)
for(i in 1:n.sims){
dt_MI[[i]] = as.data.table(dt_MI[[i]])[,ecohet_max:=max(as.numeric(as.character(unlist(varinc)))),by=c("sliceid")]
}
for(i in 1:n.sims){
dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}


# check distribution of ecohet max
par(mfrow=c(1,3)) 
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$ecohet_max)}))[,1],xlab = "",breaks=10,main = "ecohet_max = 1")
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$ecohet_max)}))[,2],xlab = "",breaks=10,main = "ecohet_max = 2")
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$ecohet_max)}))[,3],xlab = "",breaks=10,main = "ecohet_max = 3")

#center to min = 0 

for(i in 1:n.sims){
  dt_MI[[i]]$ecohet_max=(dt_MI[[i]][,which(names(dt_MI[[i]])=="ecohet_max")])-1
}

# mean variance of average annual family income across families by sliceid = ecohet_mean
library(data.table)
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.table(dt_MI[[i]])[,ecohet_mean:=mean(as.numeric(as.character(unlist(varinc)))),by=c("sliceid")]
}
for(i in 1:n.sims){
  dt_MI[[i]] = as.data.frame(dt_MI[[i]])
}



# Have separate samples for fishery and water 
for(i in 1:n.sims){
dt_MI[[i]]$fishery = ifelse(dt_MI[[i]]$resourcesector1==1,1,0) 
dt_MI[[i]]$water = ifelse(dt_MI[[i]]$resourcesector1==5,1,0)

}


table(dt_MI[[i]]$water)
table(dt_MI[[i]]$fishery)


#remove duplicates of sliceid 
for(i in 1:n.sims){
  dt_MI[[i]] = dt_MI[[i]][-which(duplicated(dt_MI[[i]]$sliceid)),]
}
#remove if not water not fishery
for(i in 1:n.sims){
dt_MI[[i]] = dt_MI[[i]][-which(dt_MI[[i]]$fishery==0 & dt_MI[[i]]$water==0),]
}


# code for histogram of imputed values of subgroupsubvar across simulations
par(mfrow=c(1,3))
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$subgroupsubvar)}))[,1],xlab = "",breaks=10,main = "subgroup-subvar = 1")
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$subgroupsubvar)}))[,2],xlab = "",breaks=10,main = "subgroup-subvar = 2")
hist(do.call('rbind',lapply(1:n.sims,function(x){table(dt_MI[[x]]$subgroupsubvar)}))[,3],xlab = "",breaks=10,main = "subgroup-subvar = 3")
dev.off()
# code for out of bag error mean and sd for each variable
# turn list into dataframe
oob_dt = foreach(i = 1:n.sims,.combine = 'rbind') %do% attributes(dt_MI[[i]])$oob


oob_mean = apply(oob_dt,2,function(x){mean(x)})
oob_sd = apply(oob_dt,2,function(x){sd(x)})
names(oob_mean)


# plot oob error for a given variable
hist(do.call('c',lapply(1:n.sims,function(x){attributes(dt_MI[[x]])$oob["subgroupsubvar"]})),xlab = "",breaks=10,main = 'oob error over simulations - subgroup-subvar')


# calculate % missing per remaining variables
pct_miss = apply(dt_temp_clean_MI,2,function(x){sum(is.na(x))/length(x)})


# correlation between pct missing and oob error
plot(pct_miss[-which(pct_miss==0)],oob_mean,main = "% NA v. OOB error",bty = 'n',xlab = "% NA")
temp= loess(oob_mean~pct_miss[-which(pct_miss==0)])
lines(pct_miss[-which(pct_miss==0)][order(pct_miss[-which(pct_miss==0)])],temp$fitted[order(pct_miss[-which(pct_miss==0)])],
      col='darkgreen')

# mean and sd oob for key variables 


# have a continuous & a factor version of unit quality and balance and trust

for(i in 1:n.sims){
  dt_MI[[i]]$unit_quality_factor = as.factor(as.character(unlist(dt_MI[[i]]$unit_quality_cont))) 
  dt_MI[[i]]$balance_cont = as.numeric(as.character(unlist(dt_MI[[i]]$balance_factor)))
}

# time to run analyses
library("lmtest")

# first analysis: test whether we can use econ het and trust as continuous
# 1. Unit quality & Trust
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
model_cont = lm(formula = unit_quality_cont ~ trust,data = dt_MI[[i]] )
model_fact = lm(formula = unit_quality_cont ~ as.factor(trust),data = dt_MI[[i]])
optimal_cov_format_trust = rbind(optimal_cov_format_trust,
data.frame(
switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50)
abline(v = 0.05,col = "red",lty = 2,lwd=2)



#2. Balance & Trust 
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
model_cont = lm(formula = balance_cont ~ trust,data = dt_MI[[i]] )
model_fact = lm(formula = balance_cont ~ as.factor(trust),data = dt_MI[[i]])
optimal_cov_format_trust = rbind(optimal_cov_format_trust,
data.frame(
switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50,xlim = c(0,1))
abline(v = 0.05,col = "red",lty = 2,lwd=2)

#3. Unit quality & economic heterogeneity
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
  model_cont = lm(formula = unit_quality_cont ~ ecohet_max,data = dt_MI[[i]] )
  model_fact = lm(formula = unit_quality_cont ~ as.factor(ecohet_max),data = dt_MI[[i]])
  optimal_cov_format_trust = rbind(optimal_cov_format_trust,
                                   data.frame(
                                     switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
                                     p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50,xlim = c(0,1))
abline(v = 0.05,col = "red",lty = 2,lwd=2)

#4. Balance & economic heterogeneity
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
  model_cont = lm(formula = balance_cont ~ ecohet_max,data = dt_MI[[i]] )
  model_fact = lm(formula = balance_cont ~ as.factor(ecohet_max),data = dt_MI[[i]])
  optimal_cov_format_trust = rbind(optimal_cov_format_trust,
                                   data.frame(
                                     switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
                                     p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50,xlim = c(0,1))
abline(v = 0.05,col = "red",lty = 2,lwd=2)

#Unit quality and sociocultural heterogeneity
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
  model_cont = lm(formula = unit_quality_cont ~ sochet_max,data = dt_MI[[i]] )
  model_fact = lm(formula = unit_quality_cont ~ as.factor(sochet_max),data = dt_MI[[i]])
  optimal_cov_format_trust = rbind(optimal_cov_format_trust,
                                   data.frame(
                                     switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
                                     p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50,xlim = c(0,1))
abline(v = 0.05,col = "red",lty = 2,lwd=2)

#Balance and sociocultural heterogeneity
optimal_cov_format_trust = data.frame()
for(i in 1:n.sims){
  model_cont = lm(formula = balance_cont ~ sochet_max,data = dt_MI[[i]] )
  model_fact = lm(formula = balance_cont ~ as.factor(sochet_max),data = dt_MI[[i]])
  optimal_cov_format_trust = rbind(optimal_cov_format_trust,
                                   data.frame(
                                     switch_to = c("continuous","factor")[which(lrtest(model_cont, model_fact)$LogLik==max(lrtest(model_cont, model_fact)$LogLik))],
                                     p_value = lrtest(model_cont, model_fact)$`Pr(>Chisq)`[2]))
}
#Is there much difference? check p values
hist(optimal_cov_format_trust$p_value,main = "p values for likelihood ratio test\ndistribution over simulations",breaks=50,xlim = c(0,1))
abline(v = 0.05,col = "red",lty = 2,lwd=2)


#--------------------------------------------------------------------------------------------------------------------------
# Creating subsample datasets
#--------------------------------------------------------------------------------------------------------------------------


# Create only fishing ground sample
dt_MI_fishery =
  foreach(i = 1:n.sims) %do%
  dt_MI[[i]][which(dt_MI[[i]]$fishery==1),]
# how big is each sim? 
dimensions_fishery = foreach(i = 1:n.sims,.combine = 'rbind') %do% dim(dt_MI_fishery[[i]])

# Create only irrigation system sample
dt_MI_water =
  foreach(i = 1:n.sims) %do%
  dt_MI[[i]][which(dt_MI[[i]]$water==1),]
# how big is each sim? 
dimensions_water = foreach(i = 1:n.sims,.combine = 'rbind') %do% dim(dt_MI_water[[i]])



#--------------------------------------------------------------------------------------------------------------------------
# Descriptive stats
#--------------------------------------------------------------------------------------------------------------------------

stat.desc(dt_MI[[i]]$unit_quality_cont)
stat.desc(dt_MI[[i]]$balance_cont)
stat.desc(dt_MI[[i]]$sochet_latent)
hist(dt_MI[[i]]$sochet_latent)
stat.desc(dt_MI[[i]]$ecohet_max)
stat.desc(dt_MI[[i]]$trust)
stat.desc(dt_MI[[i]]$sochet_max)

table(is.na(dt_temp_clean$oplevelsocstrat))
table(is.na(dt_temp_clean_MI$oplevelrelanims))
table(is.na(dt_temp_clean_MI$subgroupsubvar))
table(is.na(dt_temp_clean_MI$oplevelclanid))

#--------------------------------------------------------------------------------------------------------------------------
# Correlations
#--------------------------------------------------------------------------------------------------------------------------

# combined sample
rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$trust,  type=c("spearman"))
rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$trust, type=c("spearman"))
rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$balance_cont, type=c("spearman"))


# Irrigation system sample
rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$trust, type=c("spearman"))
rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$trust, type=c("spearman"))
rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$balance_cont, type=c("spearman"))


# Fishing ground sample
rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$trust, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$trust, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$balance_cont, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$unit_quality_cont, type=c("spearman"))
rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$balance_cont, type=c("spearman"))


#--------------------------------------------------------------------------------------------------------------------------
# Fitting models
#--------------------------------------------------------------------------------------------------------------------------


# Fitting model to each sim and pooling results


#------------ Full model OLS results --------------
# pool scalar function - home made for it to work on 
# univariate parameters from non-mice multiple imputation lists.
# See https://stefvanbuuren.name/fimd/sec-whyandwhen.html#eq:within, 
# section 2, and specifically 2.3.5 and 2.3.6
pool.scalar.homemade = function(Q, U, n, k){
  
  # m = number of imputed datasets
  m <- length(Q)
  
  # qbar = average of the estimated parameter of interest across simulations
  qbar <- mean(Q)
  
  # ubar = average of estimated variance of estimated parameter of interest across simulations
  # estimate of the variance caused by the fact that we are taking a 
  # sample rather than observing the entire population
  ubar <- mean(U)
  
  # variance of estimated parameter across smulations 
  # (the extra variance caused by the fact that there are missing values in the sample)
  b <- var(Q)
  
  # total variance according to rubin, 
  # accounting for fact that qbar is itself estimated (hence the (m+1)/m correction on b)
  t <- ubar + (m + 1) * b/m
  
  # proportion of the variation attributable to the missing data
  lambda = (b + (b/m))/t;
  
  # relative increase in variance due to nonresponse
  r = lambda/(1-lambda);
  
  # calculation for adjusted degrees of freedom, accounting for missing data
  df_old = (m-1)*(1/(lambda^2));
  df_com = n - k;
  df_obs = ((df_com + 1)/(df_com + 3))*df_com*(1-lambda)
  df = (df_old*df_obs)/(df_old + df_obs)
  
  # fraction of information about qbar missing due to nonresponse
  gamma = (r+2/(df + 3))/(1+r); 
  
  return(list(
    m = m,
    qhat = Q,
    u = U,
    qbar = qbar,
    ubar = ubar,
    b=b,
    t=t,
    df = df,
    r = r,
    lambda = lambda,
    fmi = gamma))
  
}

# identical to pool.r.squared(..., adjusted = TRUE), just it takes as value a list of fits, not a .mids object
pooled.r.squared.homemade.adj <- function(fit_sims){
  
    n <- length(fit_sims[[1]]$fitted.values) # number of observations in dataset
    k <- 1 + n - fit_sims[[1]]$df.residual # number of parameters estimated
    m <- length(fit_sims)
   r2 <- matrix(NA, nrow = m, ncol = 3, dimnames = list(seq_len(m), c("R^2", "Fisher trans F^2", "se()")))
   
# Fill arrays
for (i in seq_len(m)) {
       fit <- fit_sims[[i]]
  r2[i, 1] <- sqrt(summary(fit)$adj.r.squared)
  r2[i, 2] <- 0.5 * log((r2[i, 1] + 1)/(1 - r2[i, 1]))
  r2[i, 3] <- 1/(length(summary(fit)$residuals) - 3)
}

       fit <- pool.scalar.homemade(r2[, 2], r2[, 3],n= n, k=k)

      qbar <- fit$qbar
     table <- array(((exp(2 * qbar) - 1)/(1 + exp(2 * qbar)))^2,dim = c(1, 4))

dimnames(table) <- list("adj R^2", c("est", "lo 95", "hi 95", "fmi"))

table[, 2] <- ((exp(2 * (qbar - 1.96 * sqrt(fit$t))) - 1)/(1 + exp(2 * (qbar - 1.96 * sqrt(fit$t)))))^2
table[, 3] <- ((exp(2 * (qbar + 1.96 * sqrt(fit$t))) - 1)/(1 + exp(2 * (qbar + 1.96 * sqrt(fit$t)))))^2
table[, 4] <- fit$fmi

return(table)
}

table


#--------------------------------------------------------------------------------------------------------------------------
# Fitting models: OLS regressions
#--------------------------------------------------------------------------------------------------------------------------


#Model 1: Unit quality : irrigation
fit_sims1 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ water,data = as.data.frame(dt_MI[[i]] ))
# look into
pool(fit_sims1)
summary(pool(fit_sims1))
# find pooled R2
pooled.r.squared.homemade.adj(fit_sims1)

# get AIC
AIC_sims1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ water,data = as.data.frame(dt_MI[[i]] )))
mean(unlist((AIC_sims1[[i]][2])))



#Model 2: Balance : irrigation 
fit_sims2 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water,data = dt_MI[[i]] )
# look into
pool(fit_sims2)
summary(pool(fit_sims2))

pooled.r.squared.homemade.adj(fit_sims2)

# get AIC
AIC_sims2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water,data = dt_MI[[i]] ))
mean(unlist((AIC_sims2[[i]][2])))



#Model 3: Unit quality: irrigation, sochet, ecohet
fit_sims3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ water + sochet_max + ecohet_max,data = dt_MI[[i]] )
# look into
pool(fit_sims3)
summary(pool(fit_sims3))
pooled.r.squared.homemade.adj(fit_sims3)

#AIC
AIC_sims3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ water + sochet_max + ecohet_max,data = dt_MI[[i]] ))
mean(unlist((AIC_sims3[[i]][2])))


#model 4: Balance: irrigation, sochet, ecohet 
fit_sims4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + sochet_max + ecohet_max,data = dt_MI[[i]] )
# look into
pool(fit_sims4)
summary(pool(fit_sims4))
pooled.r.squared.homemade.adj(fit_sims4)

#AIC
AIC_sims4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + sochet_max + ecohet_max,data = dt_MI[[i]] ))
mean(unlist((AIC_sims4[[i]][2])))


#Model 5: Unit quality: irrigation, sochet, ecohet, trust
fit_sims5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ water + trust + sochet_max + ecohet_max,data = dt_MI[[i]] )
# look into
pool(fit_sims5)
summary(pool(fit_sims5))
pooled.r.squared.homemade.adj(fit_sims5)

#AIC
AIC_sims5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ water + trust + sochet_max + ecohet_max,data = dt_MI[[i]] ))
mean(unlist((AIC_sims5[[i]][2])))


#Model 6: Balance: irrigation, sochet, ecohet, trust
fit_sims6 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + trust+ sochet_max + ecohet_max + sochet_max*ecohet_max,data = dt_MI[[i]] )
# look into
pool(fit_sims6)
summary(pool(fit_sims6))
pooled.r.squared.homemade.adj(fit_sims6)

#AIC
AIC_sims6 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + trust+ sochet_max + ecohet_max,data = dt_MI[[i]] ))
mean(unlist((AIC_sims6[[i]][2])))


#Model 7: Unit quality: irrigation x trust, sochet, ecohet 
fit_sims7 = 
foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + water + ecohet_max + sochet_max + water*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims7)
summary(pool(fit_sims7))
pooled.r.squared.homemade.adj(fit_sims7)

#for fishery
fit_sims7 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + fishery + ecohet_max + sochet_max + fishery*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims7)
summary(pool(fit_sims7))
pooled.r.squared.homemade.adj(fit_sims7)


#AIC
AIC_sims7 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + water + ecohet_max + sochet_max + water*trust,data = dt_MI[[i]] ))
mean(unlist((AIC_sims7[[i]][2])))



#Model 8: Balance: irrigation x trust, sochet, ecohet
fit_sims8 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + trust+ sochet_max+ ecohet_max + water*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims8)
summary(pool(fit_sims8))
pooled.r.squared.homemade.adj(fit_sims8)

#for fisheries
fit_sims8 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ fishery + trust+ sochet_max+ ecohet_max + fishery*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims8)
summary(pool(fit_sims8))
pooled.r.squared.homemade.adj(fit_sims8)

#AIC
AIC_sims8 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + trust+ sochet_max+ ecohet_max + water*trust,data = dt_MI[[i]] ))
mean(unlist((AIC_sims8[[i]][2])))

#Model 9: Trust: ecohet, sochet 
fit_sims9 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~ sochet_max + ecohet_max,data = dt_MI[[i]] )
# look into
pool(fit_sims9)
summary(pool(fit_sims9))
pooled.r.squared.homemade.adj(fit_sims9)

#AIC
AIC_sims9 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~ sochet_max + ecohet_max,data = dt_MI[[i]] ))
mean(unlist((AIC_sims9[[i]][2])))


#Model 10: Trust: ecohet, sochet , water 
fit_sims10 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~ sochet_max + ecohet_max + water ,data = dt_MI[[i]] )
# look into
pool(fit_sims10)
summary(pool(fit_sims10))
pooled.r.squared.homemade.adj(fit_sims9)

#AIC
AIC_sims10 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( lm(formula = trust ~ sochet_max + ecohet_max + water ,data = dt_MI[[i]] ))
mean(unlist((AIC_sims10[[i]][2])))

#----------------- Separate sample results: Fishery ---------------------------------------------------------------------


#Model 1: unit quality: sochet, ecohet 
fit_sims_fish1 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ sochet_max + ecohet_max  ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish1)
summary(pool(fit_sims_fish1))
pooled.r.squared.homemade.adj(fit_sims_fish1)

#AIC
AIC_sims_fish1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish1[[i]][2])))



#Model 2: balance: sochet, ecohet
fit_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ sochet_max + ecohet_max  ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish2)
summary(pool(fit_sims_fish2))
pooled.r.squared.homemade.adj(fit_sims_fish2)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))


#Model 3: unit quality: trust, sochet, ecohet 
fit_sims_fish3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish3)
summary(pool(fit_sims_fish3))
pooled.r.squared.homemade.adj(fit_sims_fish3)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))



#Model 4: balance: trust, sochet, ecohet 
fit_sims_fish4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish4)
summary(pool(fit_sims_fish4))
pooled.r.squared.homemade.adj(fit_sims_fish4)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))


#Model 5: Trust: sochet, ecohet 
fit_sims_fish5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~  sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish5)
summary(pool(fit_sims_fish5))
pooled.r.squared.homemade.adj(fit_sims_fish5)

#AIC
AIC_sims_fish5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~  sochet_max + ecohet_max ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish5[[i]][2])))


#----------------- Separate sample results: Irrigation ---------------------------------------------------------------------


#Model 1: unit quality: sochet, ecohet 
fit_sims_water1 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ sochet_max + ecohet_max ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water1)
summary(pool(fit_sims_water1))
pooled.r.squared.homemade.adj(fit_sims_water1)

#AIC
AIC_sims_water1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ sochet_max + ecohet_max ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water1[[i]][2])))



#Model 2: balance: sochet, ecohet
fit_sims_water2 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ sochet_max + ecohet_max  ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water2)
summary(pool(fit_sims_water2))
pooled.r.squared.homemade.adj(fit_sims_water2)

#AIC
AIC_sims_water2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ sochet_max + ecohet_max ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water2[[i]][2])))


#Model 3: unit quality: trust, sochet, ecohet 
fit_sims_water3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + sochet_max + ecohet_max  ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water3)
summary(pool(fit_sims_water3))
pooled.r.squared.homemade.adj(fit_sims_water3)

#AIC
AIC_sims_water3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water3[[i]][2])))


#Model 4: balance: trust, sochet, ecohet 
fit_sims_water4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ trust + sochet_max + ecohet_max  ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water4)
summary(pool(fit_sims_water4))
pooled.r.squared.homemade.adj(fit_sims_water4)

#AIC
AIC_sims_water4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ trust + sochet_max + ecohet_max ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water4[[i]][2])))


#Model 5: Trust: sochet, ecohet 
fit_sims_water5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~  sochet_max + ecohet_max ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water5)
summary(pool(fit_sims_water5))
pooled.r.squared.homemade.adj(fit_sims_water5)

#AIC
AIC_sims_water5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~  sochet_max + ecohet_max ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water5[[i]][2])))

#------------------------------------------------------------------------------------------------------------------------------
# Indirect effect tests
#------------------------------------------------------------------------------------------------------------------------------

# WHOLE SAMPLE MAIN TRUST = FISHERY

# Ecohet --> trust and trust --> unit quality for whole sample (main trust for fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.318,0.104) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.530, 0.131) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.019
p_value_Z #0.017


# Sochet --> trust and trust --> unit quality for whole sample (main trust fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.186,0.142) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,0.530, 0.131) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.229
p_value_Z #0.226

# Ecohet --> trust and trust --> balance for whole sample (main trust fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.318,0.104) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.433, 0.220) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.114
p_value_Z #0.111


# Sochet --> trust and trust --> balance for whole sample (main trust fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.186,0.142) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,0.433, 0.220) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.316
p_value_Z #0.314

# WHOLE SAMPLE MAIN TRUST = IRRIGATION

# Ecohet --> trust and trust --> unit quality for whole sample (main trust for irrigation)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.318,0.104) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,-0.045, 0.106) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.689
p_value_Z #0.688

# Sochet --> trust and trust --> unit quality for whole sample (main trust irrigation)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.186,0.142) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,-0.045, 0.106) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.747
p_value_Z #0.746

# Ecohet --> trust and trust --> balance for whole sample (main trust irrigation)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.318,0.104) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.619, 0.198) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.036
p_value_Z #0.033


# Sochet --> trust and trust --> balance for whole sample (main trust irrigation)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.186,0.142) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,0.619, 0.198) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.249
p_value_Z #0.246

#-----------------SEPARATE SAMPLE: FISHERIES----------------------------------------------

# Ecohet --> trust and trust --> unit quality 
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.280,0.148) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.540, 0.175) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.128
p_value_Z #0.120

# Sochet --> trust and trust --> unit quality
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.055,0.174) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,0.540, 0.175) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.769
p_value_Z #0.767

# Ecohet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.280,0.148) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.498, 0.256) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.213
p_value_Z #0.205

# Sochet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.055,0.174) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,0.498, 0.256) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.786
p_value_Z #0.785


#--------------SEPARATE SAMPLE: IRRIGATION------------------------------------------------------------------------------------

# Ecohet --> trust and trust --> unit quality 
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.286,0.140) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,-0.097, 0.073) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.309
p_value_Z #0.304

# Sochet --> trust and trust --> unit quality
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.368,0.156) #sochet --> trust 
B_sims = rnorm(n.sims.indirect,-0.097, 0.073) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.284
p_value_Z #0.279

# Ecohet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.286,0.140) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.591, 0.206) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.115
p_value_Z #0.109

# Sochet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.368,0.156) #sochet --> trust 
B_sims = rnorm(n.sims.indirect, 0.591, 0.206) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.085
p_value_Z #0.079


#----------------------------------------------------------------------------------------------------------------------
# OLS including control variables
#----------------------------------------------------------------------------------------------------------------------

#make extra controls numeric
for (i in 1:n.sims) {
  dt_MI[[i]]$worstoff_cont = as.numeric(as.character(unlist(dt_MI[[i]]$oplevelworstoff)))
  dt_MI[[i]]$varspace_cont = as.numeric(as.character(unlist(dt_MI[[i]]$resourcevarspace)))
  dt_MI[[i]]$closedaccess_cont = as.numeric(as.character(unlist(dt_MI[[i]]$closedaccess)))
  dt_MI[[i]]$exitprop_cont = as.numeric(as.character(unlist(dt_MI[[i]]$exit_prop)))
}


#Model 12: Unit quality: irrigation x trust, sochet, ecohet + various controls
fit_sims12 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + water + ecohet_max + sochet_max + 
       water*trust + exitprop_cont + formal_sanctions + socialsanct + 
       physicsanct + numusers + pollution + pressure + incdependent + 
       worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]] )
# look into
pool(fit_sims12)
summary(pool(fit_sims12))
pooled.r.squared.homemade.adj(fit_sims12)

#AIC
AIC_sims12 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + water + ecohet_max + sochet_max + 
                  water*trust + exitprop_cont + formal_sanctions + socialsanct + 
                  physicsanct + numusers + pollution + pressure + incdependent + 
                  worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]] ))
mean(unlist((AIC_sims12[[i]][2])))

#Ordinal logistic
fit_sims_controlsOR = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ trust_max + water + eco_maxval + sochet_max + 
       water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
       physicsanct + numusers + pollution + pressure + incdependent + 
       worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]], 
       Hess = TRUE, start = c(4.55, 1.83, 0.93, 0.02, -5.71, 0.03, 0.9, 0.21, -0.49, 
                              0.00, -19.94, 0.18, 0.08, 0.55, -0.24, 0.38, -1.65, -3.31, 2.35))

#one part
fit_sims_controlsOR = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ trust_max + water + eco_maxval + sochet_max + 
         water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
         physicsanct 
       ,data = dt_MI[[i]], 
       Hess = TRUE, start = c(1.11, -2.65, -0.23, -0.23, -5.01, 0.31, -0.14, 0.34, -0.57, 
                                                             -3.76, 0.47))
summary(pool(fit_sims_controlsOR))


#another  part --> taking water out fixes the model which can be run now 
fit_sims_controlsOR_new = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ trust_max + eco_maxval + sochet_max + 
         water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
         physicsanct + numusers + pollution + pressure + incdependent + 
         worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont ,data = dt_MI[[i]], 
       Hess = TRUE, start = c(4.11,  0.44, -0.71, -5.71, 0.03, 0.9, 0.21, -0.49, 
                              0.00, -19.94, 0.18, 0.08, 0.55, -0.24, 0.38, -1.65, -1.65, -3.31, 2.35))
summary(pool(fit_sims_controlsOR_new))

#AIC
AIC_controls1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(  polr(formula = unit_quality_factor ~ trust_max + eco_maxval + sochet_max + 
                      water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
                      physicsanct + numusers + pollution + pressure + incdependent + 
                      worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont ,data = dt_MI[[i]], 
                    Hess = TRUE, start = c(4.11,  0.44, -0.71, -5.71, 0.03, 0.9, 0.21, -0.49, 
                                           0.00, -19.94, 0.18, 0.08, 0.55, -0.24, 0.38, -1.65, -1.65, -3.31, 2.35)))
mean(unlist((AIC_controls1[[i]][2])))

#find startvalues
fit_sims_controlsOR = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ sochet_max, data = dt_MI[[i]])
summary(pool(fit_sims_controlsOR))

#startvalue trustmax = 1.11
#startvalue water = -2.65
#startvalue ecomaxval =-0.23
#startvalue sochet_max = -0.23
#startvalue water*trust = -5.01
#startvalue exitpropcont = 0.31
#startvalue formalsanctions = -0.14
#startvalue socialsanct = 0.34
#startvalye physicsanct = -0.57
#startvelue numusers = 0.00
#startvalue pollution = -19.94
#startvalue pressure = 0.18
#startvalue incdependent = 0.08
#startvalue worstoff = 0.55
#startvalue cultvwr = -0.24
#startvalue closedaccess = 0.38
#startvalue varspace = -1.65


#Model 13: Balance: irrigation x trust, sochet, ecohet + various controls
fit_sims13 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + trust+ sochet_max + ecohet_max + 
       water*trust + + exitprop_cont + formal_sanctions + socialsanct + 
       physicsanct + numusers + pollution + pressure + incdependent + 
       worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]] )
# look into
pool(fit_sims13)
summary(pool(fit_sims13))
pooled.r.squared.homemade.adj(fit_sims13)


#AIC
AIC_controls2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( lm(formula = balance_cont ~ water + trust+ sochet_max + ecohet_max + 
                   water*trust + + exitprop_cont + formal_sanctions + socialsanct + 
                   physicsanct + numusers + pollution + pressure + incdependent + 
                   worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]] ))
mean(unlist((AIC_controls2[[i]][2])))


#ordinal logit
fit_sims_controlsOR2 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ trust_max + water + eco_maxval + sochet_max + 
         water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
         physicsanct + numusers + pollution + pressure + incdependent + 
         worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]]) 
summary(pool(fit_sims_controlsOR2))


#AIC
AIC_controls1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(    polr(formula = balance_factor ~ trust_max + water + eco_maxval + sochet_max + 
                        water*trust_max + exitprop_cont + formal_sanctions + socialsanct + 
                        physicsanct + numusers + pollution + pressure + incdependent + 
                        worstoff_cont + cultvwr_cont + closedaccess_cont + varspace_cont,data = dt_MI[[i]]) )
mean(unlist((AIC_controls1[[i]][2])))


#--------------------------------------------------------------------------------------------------------------------------
# Fitting models: Ordinal Logistic regressions
#--------------------------------------------------------------------------------------------------------------------------

#-------------- Ordinal logistic regressions: whole sample -------------------

#Model 1: Unit quality : irrigation
fit_simsOR1 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR1)
summary(pool(fit_simsOR1))

# get AIC
AIC_simsOR1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = unit_quality_factor ~ water,data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR1[[i]][2])))

#Model 2: Balance : irrigation 
fit_simsOR2 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR2)
summary(pool(fit_simsOR2))

# get AIC
AIC_simsOR2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~ water,data = dt_MI[[i]]  ))
mean(unlist((AIC_simsOR2[[i]][2])))

#Model 3: Unit quality: irrigation, sochet, ecohet
fit_simsOR3 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ water + sochet_max + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR3)
summary(pool(fit_simsOR3))

#sochet as factor
fit_simsOR3 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~ water + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR3)
summary(pool(fit_simsOR3))


# get AIC
AIC_simsOR3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~ water + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR3[[i]][2])))


#model 4: Balance: irrigation, sochet, ecohet 
fit_simsOR4 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water + sochet_max + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR4)
summary(pool(fit_simsOR4))

#sochet as factor
fit_simsOR4 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR4)
summary(pool(fit_simsOR4))

# get AIC
AIC_simsOR4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~ water + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR4[[i]][2])))

#Model 5: Unit quality: irrigation, sochet, ecohet, trust 
fit_simsOR5 = 
  foreach(i = 1:n.sims) %do%
polr(formula = unit_quality_factor ~ water + as.factor(trust) + sochet_max + as.factor(ecohet_max),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26,-0.23,-0.48,--0.50 , -3.46, 0.47 ))
# look into
pool(fit_simsOR5)
summary(pool(fit_simsOR5))

# with sochet max as factor 
fit_simsOR5 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~ water + as.factor(trust) + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26, 0.80, 0.27, 9.12, -0.64,-0.48,-0.50 , -3.46, 0.47 ))
# look into
pool(fit_simsOR5)
summary(pool(fit_simsOR5))


# get AIC
AIC_simsOR5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~ water + as.factor(trust) + as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26, 0.80, 0.27, 9.12, -0.64,-0.48,-0.50 , -3.46, 0.47 )))
mean(unlist((AIC_simsOR5[[i]][2])))

# get start values for unit quality model
fit_sims_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ as.factor(sochet_max),data = dt_MI[[i]], Hess = TRUE)
# look into
pool(fit_sims_start)
summary(pool(fit_sims_start))
#startvalue is -2.64 for water 
# startvalues for trust are 3.01 and 3.26
# startvalue for sochet is -0.23
# startvalues for ecohet are -0.48 and -0.50 
#startvalues for sochet factor are 0.80, 0.27, 9.12, -0.64


#Model 6: Balance: irrigation, sochet, ecohet, trust
fit_simsOR6 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water + as.factor(trust)+ sochet_max + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR6)
summary(pool(fit_simsOR6))

#with sochet as factor
fit_simsOR6 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water + as.factor(trust)+ as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR6)
summary(pool(fit_simsOR6))

# get AIC
AIC_simsOR6 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~ water + as.factor(trust)+ as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR6[[i]][2])))


#Model 7: Unit quality: irrigation x trust, sochet, ecohet
fit_simsOR7 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~  water + as.factor(trust) + sochet_max + as.factor(ecohet_max)  + water*as.factor(trust),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26,-0.23,-0.48,-0.50 , -0.79, -5.21, -3.46, 0.47 ) )
# look into
pool(fit_simsOR7)
summary(pool(fit_simsOR7))

#sochet as factor
fit_simsOR7 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~  water + as.factor(trust) + as.factor(sochet_max) + as.factor(ecohet_max)  + water*as.factor(trust),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26,0.80, 0.27, 9.12, -0.64,-0.48,-0.50 , -0.79, -5.21, -3.46, 0.47 ) )
# look into
pool(fit_simsOR7)
summary(pool(fit_simsOR7))


# get AIC
AIC_simsOR7 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( bayespolr(formula = unit_quality_factor ~  water + as.factor(trust) + as.factor(sochet_max) + as.factor(ecohet_max)  + water*as.factor(trust),data = dt_MI[[i]], Hess = TRUE,start = c(-2.64, 3.01, 3.26,0.80, 0.27, 9.12, -0.64,-0.48,-0.50 , -0.79, -5.21, -3.46, 0.47 ) ))
mean(unlist((AIC_simsOR7[[i]][2])))

# get start values for unit quality model
fit_sims_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ as.factor(trust)*water,data = dt_MI[[i]], Hess = TRUE)
# look into
pool(fit_sims_start)
summary(pool(fit_sims_start))
#startvalue is -2.64 for water 
# startvalues for trust are 3.01 and 3.26
# startvalue for sochet is -0.23
# startvalues for ecohet are -0.48 and -0.50 
# startvalyes for trustxwater are -0.79 and -5.21 


#Model 8: Balance: irrigation x trust, sochet, ecohet
fit_simsOR8 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ water + as.factor(trust)+ sochet_max + as.factor(as.character(unlist(ecohet_max))) + water*as.factor(trust),data = dt_MI[[i]] )
# look into
pool(fit_simsOR8)
summary(pool(fit_simsOR8))

#sochet as factor
fit_simsOR8 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = balance_factor ~ water + as.factor(trust)+ as.factor(sochet_max) + as.factor(as.character(unlist(ecohet_max))) + water*as.factor(trust),data = dt_MI[[i]] )
# look into
pool(fit_simsOR8)
summary(pool(fit_simsOR8))

# get AIC
AIC_simsOR8 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = balance_factor ~ water + as.factor(trust)+ as.factor(sochet_max) + as.factor(as.character(unlist(ecohet_max))) + water*as.factor(trust),data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR8[[i]][2])))

#make factor variable from continuous variable trust 

for (i in 1:n.sims) {
  dt_MI[[i]]$trust_factor = as.factor(as.character(unlist(dt_MI[[i]]$trust)))
}

#Model 9: Trust: ecohet, sochet 
fit_simsOR9 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = trust_factor ~ sochet_max + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR9)
summary(pool(fit_simsOR9))

#sochet as factor
fit_simsOR9 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] )
# look into
pool(fit_simsOR9)
summary(pool(fit_simsOR9))

# get AIC
AIC_simsOR9 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max),data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR9[[i]][2])))

# Model 10: Trust: ecohet sochet water
fit_simsOR10 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = trust_factor ~ sochet_max + as.factor(ecohet_max) + water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR10)
summary(pool(fit_simsOR10))

# with sochet as factor
fit_simsOR10 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max) + water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR10)
summary(pool(fit_simsOR10))

# get AIC
AIC_simsOR10 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max) + water,data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR10[[i]][2])))

#Create factor variables from continuous variables
for (i in 1:n.sims) {
  dt_MI[[i]]$sochet_factor = as.factor(as.character(unlist(dt_MI[[i]]$sochet_max)))
  dt_MI[[i]]$ecohet_factor = as.factor(as.character(unlist(dt_MI[[i]]$ecohet_max)))
}

dt_MI[[i]]$sochet_max

#Model 10: Sochet: irrigation
fit_simsOR10 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = sochet_factor ~ water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR10)
summary(pool(fit_simsOR10))

#Model 10b OLS: sochet : irrigation
fit_simsOLS10 =
  foreach(i = 1:n.sims) %do%
  lm(formula = sochet_max ~water, data = dt_MI[[i]])
pool(fit_simsOLS10)
summary(pool(fit_simsOLS10))

#Model 11: Ecohet: irrigation
fit_simsOR11 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = ecohet_factor ~ water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR11)
summary(pool(fit_simsOR11))

#------------------- Ordinal logistic regressions: FISHERIES ----------------------------------------------------------------------------

# Make factor variable from factor variable trust 
for (i in 1:n.sims) {
  dt_MI_fishery[[i]]$trust_factor = as.factor(as.character(unlist(dt_MI_fishery[[i]]$trust)))
}

#Model 1: unit quality: sochet ecohet
fit_sims_fishOR1 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]], Hess = TRUE, start = c(0.17, -1.21, -0.89, -0.21, 2.37 ))
#look into
pool(fit_sims_fishOR1)
summary(pool(fit_sims_fishOR1))
# startvalue of sochet is -0.17
# startvalues of ecohet are -1.21 and -0.89

#sochet factor
fit_sims_fishOR1 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~ as.factor(sochet_max) + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR1)
summary(pool(fit_sims_fishOR1))


# get AIC
AIC_simsfishOR1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]], Hess = TRUE, start = c(0.17, -1.21, -0.89, -0.21, 2.37 )))
mean(unlist((AIC_simsfishOR1[[i]][2])))

#Model 2: balance: sochet ecohet
fit_sims_fishOR2 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR2)
summary(pool(fit_sims_fishOR2))

# get AIC
AIC_simsfishOR2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR2[[i]][2])))

#Model 3: unit quality: sochet ecohet trust 
fit_sims_fishOR3 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR3)
summary(pool(fit_sims_fishOR3))

#sochey factor
fit_sims_fishOR3 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~as.factor(sochet_max) + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR3)
summary(pool(fit_sims_fishOR3))

# get AIC
AIC_simsfishOR3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR3[[i]][2])))


#Model 4: unit quality:  trust 
fit_sims_fishOR4 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~as.factor(trust), data = dt_MI_fishery[[i]], Hess = TRUE)
#look into
pool(fit_sims_fishOR4)
summary(pool(fit_sims_fishOR4))

# get AIC
AIC_simsfishOR4= 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = unit_quality_factor ~as.factor(trust), data = dt_MI_fishery[[i]], Hess = TRUE))
mean(unlist((AIC_simsfishOR4[[i]][2])))

# finding startvalues for unit quality models
fit_sims_fish_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ as.factor(trust) , data = dt_MI_fishery[[i]], Hess = TRUE)
#look into
pool(fit_sims_fish_start)
summary(pool(fit_sims_fish_start))
# startvalues of trust are 1.95 and 3.45
# startvalue of sochet is -0.17
# startvalues of ecohet are -1.21 and -0.89


#Model 5: balance:  trust 
fit_sims_fishOR5 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~as.factor(trust), data = dt_MI_fishery[[i]], Hess = TRUE)
#look into
pool(fit_sims_fishOR5)
summary(pool(fit_sims_fishOR5))

# get AIC
AIC_simsfishOR5= 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = balance_factor ~as.factor(trust), data = dt_MI_fishery[[i]], Hess = TRUE))
mean(unlist((AIC_simsfishOR5[[i]][2])))

#Model 6: balance: sochet ecohet trust
fit_sims_fishOR__6 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR__6)
summary(pool(fit_sims_fishOR__6))

# get AIC
AIC_simsfishOR6= 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR6[[i]][2])))


#Model 7: trust: sochet ecohet (with Bayespolr otherwise standard errors are gigantic)
fit_sims_fishOR7 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = trust_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR7)
summary(pool(fit_sims_fishOR7))

# get AIC
AIC_simsfishOR7= 
  foreach(i = 1:n.sims) %do%
  extractAIC( bayespolr(formula = trust_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR7[[i]][2])))



# finding startvalyes for trust model
fit_sims_fish_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = trust_factor ~ as.factor(ecohet_max) , data = dt_MI_fishery[[i]], Hess = TRUE)
#look into
pool(fit_sims_fish_start)
summary(pool(fit_sims_fish_start))

# startvalue for sochet is -0.02
# startvalues for sochet are -14.95 and -15.97



#------------------------- Ordinal logistic regressions: IRRIGATION ----------------------------------------------------------------------------
# Make factor variables from continuous variable trust 
for (i in 1:n.sims) {
  dt_MI_water[[i]]$trust_factor = as.factor(as.character(unlist(dt_MI_water[[i]]$trust)))
}


#Model 1: unit quality: ecohet 
fit_sims_irrOR1 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ as.factor(ecohet_max), data = dt_MI_water[[i]], Hess = TRUE, start = c( 2.78, -4.48, -7.93, 5.23))
#look into
pool(fit_sims_irrOR1)
summary(pool(fit_sims_irrOR1))


# get AIC
AIC_simsirrOR1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = unit_quality_factor ~ as.factor(ecohet_max), data = dt_MI_water[[i]], Hess = TRUE, start = c( 2.78, -4.48, -7.93, 5.23)))
mean(unlist((AIC_simsirrOR1[[i]][2])))

#Model 1b: unit quality: ecohet 
fit_sims_irrOR1b = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~ eco_maxval, data = dt_MI_water[[i]], Hess = TRUE, start = c( 2.78, -7.93, 5.23))
#look into
pool(fit_sims_irrOR1b)
summary(pool(fit_sims_irrOR1b))

# getting startvalues for unit quality model
fit_sims_irr_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~as.factor(ecohet_max), data = dt_MI_water[[i]])
#look into
pool(fit_sims_irr_start)
summary(pool(fit_sims_irr_start))
# startvalue for sochet is -1.92
# startvalues for trust are 24.50 and 12.58
# startvalues for ecohet are 2.78 and -4.48


#Model 2: balance:  ecohet
fit_sims_irrOR3 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~as.factor(ecohet_max), data = dt_MI_water[[i]], Hess = TRUE, start = c( -0.95, -2.14, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR3)
summary(pool(fit_sims_irrOR3))

#Model 2b: balance:  ecohet
fit_sims_irrOR3b = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~eco_maxval, data = dt_MI_water[[i]], Hess = TRUE, start = c( -0.95, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR3b)
summary(pool(fit_sims_irrOR3b))



#Model 3: unit quality: sochet ecohet 
fit_sims_irrOR3 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR3)
summary(pool(fit_sims_irrOR3))

# get AIC
AIC_simsirrOR3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]]))
#look into)
mean(unlist((AIC_simsirrOR3[[i]][2])))

#Model 3b: unit quality: sochet ecohet 
fit_sims_irrOR3b = 
  foreach(i = 1:n.sims) %do%
  polr(formula = unit_quality_factor ~sochet_max + eco_maxval, data = dt_MI_water[[i]], Hess = TRUE, start = c(-1.92, 2.78, -7.93, 5.23))
#look into
pool(fit_sims_irrOR3b)
summary(pool(fit_sims_irrOR3b))

#Model 3c: unit quality: trust
fit_sims_irrOR3c = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~trust_factor, data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR3c)
summary(pool(fit_sims_irrOR3c))

# get AIC
AIC_simsirrOR3c = 
  foreach(i = 1:n.sims) %do%
  extractAIC( bayespolr(formula = unit_quality_factor ~trust_factor, data = dt_MI_water[[i]]))
#look into)
mean(unlist((AIC_simsirrOR3c[[i]][2])))

#Model 4: balance: sochet ecohet 
fit_sims_irrOR4 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]], Hess = TRUE, start = c(-0.82, -0.95, -2.14, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR4)
summary(pool(fit_sims_irrOR4))

# get AIC
AIC_simsirrOR4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]], Hess = TRUE, start = c(-0.82, -0.95, -2.14, -2.94, -1.22, 1.60)))
#look into)
mean(unlist((AIC_simsirrOR4[[i]][2])))

#Model 4c: balance: sochet ecohet 
fit_sims_irrOR4c = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~sochet_max + eco_maxval, data = dt_MI_water[[i]], Hess = TRUE, start = c(-0.82, -0.95, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR4c)
summary(pool(fit_sims_irrOR4c))

#Model 4b: balance: trust
fit_sims_irrOR4 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~trust_factor, data = dt_MI_water[[i]], Hess = TRUE, start = c(3.32, 4.74, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR4)
summary(pool(fit_sims_irrOR4))

# get AIC
AIC_simsirrOR4b = 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = balance_factor ~trust_factor, data = dt_MI_water[[i]], Hess = TRUE, start = c(3.32, 4.74, -2.94, -1.22, 1.60)))
#look into)
mean(unlist((AIC_simsirrOR4b[[i]][2])))

# getting startvalues for balace models 
fit_sims_irr_start = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ trust_factor, data = dt_MI_water[[i]], Hess = TRUE, start = c(0.61, 1.44, 3.84, 6.89))
#look into
pool(fit_sims_irr_start)
summary(pool(fit_sims_irr_start))
# startvalue for sochet is -0.82
# startvalues for trust are 3.32 and 4.74
# startvalues for ecohet are -0.95 and -2.14

#Model 5: unit quality: sochet ecohet trust --> solved with bayespolr
fit_sims_irrOR5 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max) + trust_factor, data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR5)
summary(pool(fit_sims_irrOR5)) 

# get AIC
AIC_simsirrOR5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max) + trust_factor, data = dt_MI_water[[i]]))
#look into)
mean(unlist((AIC_simsirrOR5[[i]][2])))


#Model 6: balance: sochet ecohet trust
fit_sims_irrOR6 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_water[[i]], Hess = TRUE,start = c(-0.82, -0.95, -2.14, 3.32, 4.74, -2.94, -1.22, 1.60))
#look into
pool(fit_sims_irrOR6)
summary(pool(fit_sims_irrOR6))

# get AIC
AIC_simsirrOR6 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_water[[i]], Hess = TRUE,start = c(-0.82, -0.95, -2.14, 3.32, 4.74, -2.94, -1.22, 1.60)))
#look into)
mean(unlist((AIC_simsirrOR6[[i]][2])))


#Model 7a: trust: sochet ecohet 
fit_sims_irrOR7 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = trust_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR7)
summary(pool(fit_sims_irrOR7))


# get AIC
AIC_simsirrOR7 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = trust_factor ~sochet_max + as.factor(ecohet_max), data = dt_MI_water[[i]]))
#look into)
mean(unlist((AIC_simsirrOR7[[i]][2])))



#Model 7a: trust: sochet 
fit_sims_irrOR7 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = trust_factor ~sochet_max, data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR7)
summary(pool(fit_sims_irrOR7))

#Model 7b: trust: ecohet
fit_sims_irrOR7 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = trust_factor ~as.factor(ecohet_max), data = dt_MI_water[[i]])
#look into
pool(fit_sims_irrOR7)
summary(pool(fit_sims_irrOR7))

dt_MI[[i]]$trust_factor
dt_MI[[i]]$trust

hist(dt_MI[[i]]$trust)




#--------------------------------------------------------------------------------------------------------------------------
# Robustness checks: means of economic and sociocultural heterogeneity instead of max
#--------------------------------------------------------------------------------------------------------------------------


#Model 3: Unit quality: irrigation, sochet, ecohet
fit_sims3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ water + sochet_mean + ecohet_mean  ,data = dt_MI[[i]] )
# look into
pool(fit_sims3)
summary(pool(fit_sims3))
pooled.r.squared.homemade.adj(fit_sims3)

#AIC
AIC_sims3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ water + sochet_mean + ecohet_mean,data = dt_MI[[i]] ))
mean(unlist((AIC_sims3[[i]][2])))


#model 4: Balance: irrigation, sochet, ecohet 
fit_sims4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + sochet_mean + ecohet_mean ,data = dt_MI[[i]] )
# look into
pool(fit_sims4)
summary(pool(fit_sims4))
pooled.r.squared.homemade.adj(fit_sims4)

#AIC
AIC_sims4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + sochet_mean + ecohet_mean,data = dt_MI[[i]] ))
mean(unlist((AIC_sims4[[i]][2])))


#Model 5: Unit quality: irrigation, sochet, ecohet, trust
fit_sims5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ water + trust + sochet_mean + ecohet_mean ,data = dt_MI[[i]] )
# look into
pool(fit_sims5)
summary(pool(fit_sims5))
pooled.r.squared.homemade.adj(fit_sims5)

#AIC
AIC_sims5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ water + trust + sochet_mean + ecohet_mean,data = dt_MI[[i]] ))
mean(unlist((AIC_sims5[[i]][2])))


#Model 6: Balance: irrigation, sochet, ecohet, trust
fit_sims6 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + trust+ sochet_mean + ecohet_mean ,data = dt_MI[[i]] )
# look into
pool(fit_sims6)
summary(pool(fit_sims6))
pooled.r.squared.homemade.adj(fit_sims6)

#AIC
AIC_sims6 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + trust+ sochet_mean + ecohet_mean,data = dt_MI[[i]] ))
mean(unlist((AIC_sims6[[i]][2])))


#Model 7: Unit quality: irrigation x trust, sochet, ecohet 
fit_sims7 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + water + ecohet_mean + sochet_mean + water*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims7)
summary(pool(fit_sims7))
pooled.r.squared.homemade.adj(fit_sims7)

#for fishery
fit_sims7 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + fishery + ecohet_mean + sochet_mean + fishery*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims7)
summary(pool(fit_sims7))
pooled.r.squared.homemade.adj(fit_sims7)


#AIC
AIC_sims7 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + water + ecohet_mean + sochet_mean + water*trust,data = dt_MI[[i]] ))
mean(unlist((AIC_sims7[[i]][2])))



#Model 8: Balance: irrigation x trust, sochet, ecohet
fit_sims8 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ water + trust+ sochet_mean+ ecohet_mean + water*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims8)
summary(pool(fit_sims8))
pooled.r.squared.homemade.adj(fit_sims8)

#for fishery
fit_sims8_fish = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ fishery + trust+ sochet_mean+ ecohet_mean + fishery*trust,data = dt_MI[[i]] )
# look into
pool(fit_sims8_fish)
summary(pool(fit_sims8_fish))
pooled.r.squared.homemade.adj(fit_sims8_fish)


#AIC
AIC_sims8 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ water + trust+ sochet_mean+ ecohet_mean + water*trust,data = dt_MI[[i]] ))
mean(unlist((AIC_sims8[[i]][2])))

#Model 9: Trust: ecohet, sochet 
fit_sims9 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~ sochet_mean + ecohet_mean,data = dt_MI[[i]] )
# look into
pool(fit_sims9)
summary(pool(fit_sims9))
pooled.r.squared.homemade.adj(fit_sims9)

#AIC
AIC_sims9 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~ sochet_mean + ecohet_mean,data = dt_MI[[i]] ))
mean(unlist((AIC_sims9[[i]][2])))


#Model 10: Trust: ecohet, sochet , water 
fit_sims10 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~ sochet_mean + ecohet_mean + water ,data = dt_MI[[i]] )
# look into
pool(fit_sims10)
summary(pool(fit_sims10))
pooled.r.squared.homemade.adj(fit_sims9)

#AIC
AIC_sims10 = 
  foreach(i = 1:n.sims) %do%
  extractAIC( lm(formula = trust ~ sochet_mean + ecohet_mean + water ,data = dt_MI[[i]] ))
mean(unlist((AIC_sims10[[i]][2])))

#----------------- Separate sample results: Fishery ---------------------------------------------------------------------


#Model 1: unit quality: sochet, ecohet 
fit_sims_fish1 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish1)
summary(pool(fit_sims_fish1))
pooled.r.squared.homemade.adj(fit_sims_fish1)

#AIC
AIC_sims_fish1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish1[[i]][2])))



#Model 2: balance: sochet, ecohet
fit_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish2)
summary(pool(fit_sims_fish2))
pooled.r.squared.homemade.adj(fit_sims_fish2)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))


#Model 3: unit quality: trust, sochet, ecohet 
fit_sims_fish3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish3)
summary(pool(fit_sims_fish3))
pooled.r.squared.homemade.adj(fit_sims_fish3)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))



#Model 4: balance: trust, sochet, ecohet 
fit_sims_fish4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish4)
summary(pool(fit_sims_fish4))
pooled.r.squared.homemade.adj(fit_sims_fish4)

#AIC
AIC_sims_fish2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish2[[i]][2])))


#Model 5: Trust: sochet, ecohet 
fit_sims_fish5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~  sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] )
# look into
pool(fit_sims_fish5)
summary(pool(fit_sims_fish5))
pooled.r.squared.homemade.adj(fit_sims_fish5)

#AIC
AIC_sims_fish5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~  sochet_mean + ecohet_mean ,data = dt_MI_fishery[[i]] ))
mean(unlist((AIC_sims_fish5[[i]][2])))


#----------------- Separate sample results: Irrigation ---------------------------------------------------------------------

#Model 1: unit quality: sochet, ecohet 
fit_sims_water1 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water1)
summary(pool(fit_sims_water1))
pooled.r.squared.homemade.adj(fit_sims_water1)

#AIC
AIC_sims_water1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water1[[i]][2])))



#Model 2: balance: sochet, ecohet
fit_sims_water2 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water2)
summary(pool(fit_sims_water2))
pooled.r.squared.homemade.adj(fit_sims_water2)

#AIC
AIC_sims_water2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water2[[i]][2])))


#Model 3: unit quality: trust, sochet, ecohet 
fit_sims_water3 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = unit_quality_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water3)
summary(pool(fit_sims_water3))
pooled.r.squared.homemade.adj(fit_sims_water3)

#AIC
AIC_sims_water3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = unit_quality_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water3[[i]][2])))


#Model 4: balance: trust, sochet, ecohet 
fit_sims_water4 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = balance_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water4)
summary(pool(fit_sims_water4))
pooled.r.squared.homemade.adj(fit_sims_water4)

#AIC
AIC_sims_water4 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = balance_cont ~ trust + sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water4[[i]][2])))


#Model 5: Trust: sochet, ecohet 
fit_sims_water5 = 
  foreach(i = 1:n.sims) %do%
  lm(formula = trust ~  sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] )
# look into
pool(fit_sims_water5)
summary(pool(fit_sims_water5))
pooled.r.squared.homemade.adj(fit_sims_water5)

#AIC
AIC_sims_water5 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(lm(formula = trust ~  sochet_mean + ecohet_mean ,data = dt_MI_water[[i]] ))
mean(unlist((AIC_sims_water5[[i]][2])))




#------------------------------------------------------------------------------------------------------------------------------
# Robustness checks: Indirect effect tests
#------------------------------------------------------------------------------------------------------------------------------

# WHOLE SAMPLE MAIN TRUST = FISHERY

# Ecohet --> trust and trust --> unit quality for whole sample (main trust for fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.363,0.106) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.543, 0.129) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.011
p_value_Z #0.009



# Ecohet --> trust and trust --> balance for whole sample (main trust fishery)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.364,0.106) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.444, 0.218) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.093
p_value_Z #0.089



# WHOLE SAMPLE MAIN TRUST = IRRIGATION




# Ecohet --> trust and trust --> balance for whole sample (main trust irrigation)
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.363,0.106) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.605, 0.199) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 90),pt(q = Z_sims,df = 90))
p_value_t #0.029
p_value_Z #0.026



#-----------------SEPARATE SAMPLE: FISHERIES----------------------------------------------

# Ecohet --> trust and trust --> unit quality 
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.288,0.164) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.528, 0.174) # trust --> unit quality 
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.154
p_value_Z #0.146



# Ecohet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.288,0.164) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.532, 0.251) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 38),pt(q = Z_sims,df = 38))
p_value_t #0.212
p_value_Z #0.204



#--------------SEPARATE SAMPLE: IRRIGATION------------------------------------------------------------------------------------



# Ecohet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.372,0.148) #ecohet --> trust 
B_sims = rnorm(n.sims.indirect,0.570, 0.195) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.071
p_value_Z #0.065

# Sochet --> trust and trust --> balance
n.sims.indirect = 100000
A_sims = rnorm(n.sims.indirect,-0.593,0.351) #sochet --> trust 
B_sims = rnorm(n.sims.indirect, 0.570, 0.195) # trust --> balance
C_sims = A_sims * B_sims
sd_C_sims = sd(C_sims)
mean_C_sims = mean(C_sims)

Z_sims = mean_C_sims/sd_C_sims
p_value_Z = 2*min(1-pnorm(q = Z_sims),pnorm(q = Z_sims))
library(stats)
p_value_t = 2*min(1-pt(q = Z_sims,df = 50),pt(q = Z_sims,df = 50))
p_value_t #0.167
p_value_Z #0.161




#--------------------------------------------------------------------------------------------------------------------------
# Available case analysis
#--------------------------------------------------------------------------------------------------------------------------

#Load available case data
setwd("G:/Nuffield DPhil/Paper 1/Data")
load(file = "dt_temp_clean_MI.RData")


dt_available = dt_temp_clean_MI
head(dt_available)
names(dt_available)[which(names(dt_available)=="opleveloplevid")]="sliceid"


#available case sochet

  names(dt_available)[which(names(dt_available)=="oplevelraceid")]="raceid"
  names(dt_available)[which(names(dt_available)=="oplevelrelanims")]="relanims"
  names(dt_available)[which(names(dt_available)=="oplevelethncid")]="ethncid"
  names(dt_available)[which(names(dt_available)=="oplevelclanid")]="clanid"
  names(dt_available)[which(names(dt_available)=="oplevelcommlang")]="commlang"
  names(dt_available)[which(names(dt_available)=="oplevelsex")]="sex"
  names(dt_available)[which(names(dt_available)=="oplevelsocstrat")]="casteid"
  names(dt_available)[which(names(dt_available)=="oplevelcultvwr")]="cultvwr"

  dt_available$raceid_cont = as.numeric(as.character(unlist(dt_available$raceid)))
  dt_available$relanims_cont = as.numeric(as.character(unlist(dt_available$relanims)))
  dt_available$ethncid_cont = as.numeric(as.character(unlist(dt_available$ethncid)))
  dt_available$commlang_cont = as.numeric(as.character(unlist(dt_available$commlang)))
  dt_available$sex_cont = as.numeric(as.character(unlist(dt_available$sex)))
  dt_available$clanid_cont = as.numeric(as.character(unlist(dt_available$clanid)))
  dt_available$casteid_cont = as.numeric(as.character(unlist(dt_available$casteid)))
  dt_available$cultvwr_cont = as.numeric(as.character(unlist(dt_available$cultvwr)))

rm_inf = function(x){ifelse(x=="-Inf",NA,x)}
dt_available$sochet_max_temp = rm_inf(apply(dt_available[,c("raceid_cont","relanims_cont","ethncid_cont","commlang_cont","sex_cont","clanid_cont", "casteid_cont")],1,max,na.rm=TRUE)) #make sochet max 

# by sliceid
library(data.table)
dt_available = as.data.table(dt_available)[,sochet_max:=max(sochet_max_temp),by=sliceid]
dt_available = as.data.frame(dt_available)
table(dt_available$sochet_max)

dt_available$sochet_max=(dt_available[,which(names(dt_available)=="sochet_max")])-1 #set min to 0 

#available case sochet as factor
dt_available$sochet_factor = as.factor(as.character(unlist(dt_available$sochet_max)))



#create available case sochet as mean of all seven categories 
for(i in 1:n.sims){
  dt_available$sochet_mean = apply(dt_available[,c("raceid_cont","relanims_cont","ethncid_cont","commlang_cont","sex_cont","clanid_cont", "casteid_cont")],1,mean)
}




#available case unit quality and balance
names(dt_available)[which(names(dt_available)=="oplevelendblnc")]="balance_factor" #make factor balance
levels(dt_available$balance_factor)=c("0", "1", "2", "3") #set min to 0 
dt_available$balance_cont = as.numeric(as.character(unlist(dt_available$balance_factor))) #make continuous balance



dt_available$unit_quality_cont=5-as.numeric(as.character(unlist(dt_available[,which(names(dt_available)=="oplevelendqual")]))) #make continuous unit quality with min to 0 
dt_available$unit_quality_factor = as.factor(as.character(unlist(dt_available$unit_quality_cont))) #make factor unit quality 

#available case trust 
dt_available$trust=3-as.numeric(as.character(unlist(dt_available[,which(names(dt_available)=="oplevelendtrust")])))
dt_available$trust_factor = as.factor(as.character(unlist(dt_available$trust)))


#available case water and fishery
  dt_available$fishery = ifelse(dt_available$resourcesector1==1,1,0)
  dt_available$water = ifelse(dt_available$resourcesector1==5,1,0)

table(dt_available$water)
table(dt_available$fishery)
table(dt_available$balance_cont)
table(dt_available$unit_quality_cont)
table(dt_available$trust)

#remove duplicates of sliceid 
dt_available = dt_available[-which(duplicated(dt_available$sliceid)),]
#remove if not water not fishery
dt_available = dt_available[-which(dt_available$fishery==0 & dt_available$water==0),]


#Model 1: unit quality: water
model1 = lm(formula = unit_quality_cont ~water, data = dt_available)
dt_available
summary(model1)
length(model1$fitted.values)

table(dt_available$water)
#OR
model1OR = polr(formula = unit_quality_factor ~water, data = dt_available)
summary(model1OR)
(ctable <- coef(summary(model1OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p
dim(model1OR$fitted.values)

#AIC
extractAIC(model1)

#Model 2: balance: water
model2 = lm(formula = balance_cont ~water, data = dt_available)
summary(model2)
length(model2$fitted.values)


model2OR = polr(formula = balance_factor ~water, data = dt_available)
summary(model2OR)
(ctable <- coef(summary(model2OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model2OR$fitted.values)

#AIC
extractAIC(model2)

#Model 3: unit quality: water sochet 
model3 = lm(formula = unit_quality_cont ~ sochet_max + water, data = dt_available)
summary(model3)
length(model3$fitted.values)

model3OR = polr(formula = unit_quality_factor ~ sochet_max+ water, data = dt_available)
summary(model3OR)
(ctable <- coef(summary(model3OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model3OR$fitted.values)

#AIC

extractAIC(model3)



#Model 4: balance: water sochet 
model4 = lm(formula = balance_cont ~sochet_max + water, data = dt_available)
summary(model4)
length(model4$fitted.values)

model4OR = polr(formula = balance_factor ~sochet_max + water, data = dt_available)
summary(model4OR)
(ctable <- coef(summary(model4OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model4OR$fitted.values)

#AIC
extractAIC(model4)

#Model 5: Unit quality: water sochet trust
model5 = lm(formula = unit_quality_cont ~ trust + sochet_max + water, data = dt_available)
summary(model5)
length(model5$fitted.values)

model5OR = polr(formula = unit_quality_factor ~ trust + sochet_max+ water, data = dt_available)
summary(model5OR)
(ctable <- coef(summary(model5OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model5OR$fitted.values)

#AIC
extractAIC(model5)


#Model 6: balance: water sochet trust
model6 = lm(formula = balance_cont ~ trust + sochet_max + water, data = dt_available)
summary(model6)
length(model6$fitted.values)

model6OR = polr(formula = balance_factor ~trust + sochet_max + water, data = dt_available)
summary(model6OR)
(ctable <- coef(summary(model6OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model6OR$fitted.values)

#AIC
extractAIC(model6)

#Model 7: unit quality: water*trust sochet 
model7 = lm(formula = unit_quality_cont ~ trust*water + trust + sochet_max + water, data = dt_available)
summary(model7)
length(model7$fitted.values)

model7OR = polr(formula = unit_quality_factor ~ trust*water + trust + sochet_max+ water, data = dt_available)
summary(model7OR)
(ctable <- coef(summary(model7OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model7OR$fitted.values)

#AIC
extractAIC(model7)


#model 8 : balance: water*trust sochet 
model8 = lm(formula = balance_cont ~ trust*water + trust + sochet_max + water, data = dt_available)
summary(model8)
length(model8$fitted.values)

model8OR = polr(formula = balance_factor ~trust*water + trust + sochet_max + water, data = dt_available)
summary(model8OR)
(ctable <- coef(summary(model8OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model8OR$fitted.values)

#AIC
extractAIC(model8)

#model 9: trust: sochet
model9 = lm(formula = trust ~ sochet_max, data = dt_available)
summary(model9)
length(model9$fitted.values)

model9OR = polr(formula = trust_factor ~sochet_max, data = dt_available)
summary(model9OR)
(ctable <- coef(summary(model9OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p

dim(model9OR$fitted.values)

#AIC
extractAIC(model9)

#model 10: sochet: water 
model10 = lm(formula = sochet_max ~water, data = dt_available)
summary(model10)
length(model10$fitted.values)

model10OR = polr(formula = sochet_factor ~water, data = dt_available)
summary(model10OR)
(ctable <- coef(summary(model10OR)))
p = pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
p


