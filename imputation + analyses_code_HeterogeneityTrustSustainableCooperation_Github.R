rm(list=ls())
options(scipen=999)
library(foreign)
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
library(gmodels)
library(data.table)
library(mediation)
library(Hmisc)
n.sims = 100

#-----------------------------------------------------------------------------------------------------------------------------
# Prepare data for Multiple Imputation procedure
#-----------------------------------------------------------------------------------------------------------------------------

dt = read.table(file = 'file',sep = '\t',header=TRUE)
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
setwd("directory")
save(dt_temp_clean_MI, file = "dt_temp_clean_MI.RData", compress = TRUE)

#Save imputed data
setwd("directory")
save(dt_MI,file = "dt_MI.RData",compress=TRUE)

#Load imputed and available data to create variables 
rm(list=ls())
setwd("directory")
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
# Descriptive stats of control vars irrigation vs fishery
#--------------------------------------------------------------------------------------------------------------------------


# Get average (mode/mean) data over 100 dataframes

fun.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
rowMode = function(x){apply(x,1,fun.mode)}


# averaging datasets
dt_av = data.table(array(NA,dim(dt_MI[[1]])))
colnames(dt_av) = colnames(dt_MI[[1]])
for(j in 1:dim(dt_MI[[1]])[2]) {
  if( class(dt_MI[[1]][,j]) =="numeric" & length(unique(dt_MI[[1]][,j]))>10 |  class(dt_MI[[1]][,j]) =="integer" & length(unique(dt_MI[[1]][,j]))>10
  ){
    temp = rowMeans(foreach(i = 1:100,.combine = 'cbind') %do% as.numeric(as.character(unlist(dt_MI[[i]][,j]))))
  }else{
    temp = 
      if(class(dt_MI[[1]][,j]) == "factor") {
        as.factor(rowMode(foreach(i = 1:100,.combine = 'cbind') %do% as.character(unlist(dt_MI[[i]][,j]))))
      }else{
        as.numeric(rowMode(foreach(i = 1:100,.combine = 'cbind') %do% as.character(unlist(dt_MI[[i]][,j]))))
      }
  }
  dt_av[,j] = temp
}



CrossTable(dt_av$resourcevarspace, dt_av$water)
CrossTable(dt_av$closedaccess, dt_av$water)
CrossTable(dt_av$exit_prop, dt_av$water)


CrossTable(dt_av$resourcevaryear, dt_av$water)
CrossTable(dt_av$resourcevarotime, dt_av$water)


CrossTable(dt_av$resourcepredvar1, dt_av$water)
CrossTable(dt_av$resourcepredvar2, dt_av$water)
CrossTable(dt_av$resourcepredvar3, dt_av$water)
CrossTable(dt_av$resourcequalbetr, dt_av$water)



#--------------------------------------------------------------------------------------------------------------------------
# Correlations on imputed datasets
#--------------------------------------------------------------------------------------------------------------------------

#Combined sample
#EH Trust
cor_cet = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$trust, type = c("spearman"))$r[1,2]
p_cet = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_cet)
sd(cor_cet)
mean(p_cet)

#SH Trust
cor_cst = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$trust, type = c("spearman"))$r[1,2]
p_cst = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_cst)
sd(cor_cst)
mean(p_cst)


#EH UQ
cor_ceu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_ceu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_ceu)
sd(cor_ceu)
mean(p_ceu)

#EH B
cor_ceb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_ceb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$ecohet_max, dt_MI[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_ceb)
sd(cor_ceb)
mean(p_ceb)

#SH UQ
cor_csu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_csu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_csu)
sd(cor_csu)
mean(p_csu)

#SH B
cor_csb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_csb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$sochet_max, dt_MI[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_csb)
sd(cor_csb)
mean(p_csb)


#Trust UQ
cor_ctu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_ctu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_ctu)
sd(cor_ctu)
mean(p_ctu)

#Trust B
cor_ctb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_ctb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI[[i]]$trust, dt_MI[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_ctb)
sd(cor_ctb)
mean(p_ctb)

#Fishing grounds

#EH Trust
cor_fet = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$trust, type = c("spearman"))$r[1,2]
p_fet = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_fet)
sd(cor_fet)
mean(p_fet)

#SH Trust
cor_fst = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$trust, type = c("spearman"))$r[1,2]
p_fst = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_fst)
sd(cor_fst)
mean(p_fst)

#EH UQ
cor_feu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_feu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_feu)
sd(cor_feu)
mean(p_feu)

#EH B
cor_feb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_feb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$ecohet_max, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_feb)
sd(cor_feb)
mean(p_feb)


#SH UQ
cor_fsu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_fsu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_fsu)
sd(cor_fsu)
mean(p_fsu)

#SH B
cor_fsb= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_fsb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$sochet_max, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_fsb)
sd(cor_fsb)
mean(p_fsb)


#Trust UQ
cor_ftu= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_ftu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_ftu)
sd(cor_ftu)
mean(p_ftu)



#Trust B
cor_ftb= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_ftb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_fishery[[i]]$trust, dt_MI_fishery[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_ftb)
sd(cor_ftb)
mean(p_ftb)

# Irrigation systems

#EH Trust
cor_iet= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$trust, type = c("spearman"))$r[1,2]
p_iet = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_iet)
sd(cor_iet)
mean(p_iet)

#SH Trust
cor_ist= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$trust, type = c("spearman"))$r[1,2]
p_ist = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$trust, type = c("spearman"))$P[1,2]
mean(cor_ist)
sd(cor_ist)
mean(p_ist)

#EH UQ
cor_ieu= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_ieu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_ieu)
sd(cor_ieu)
mean(p_ieu)

#EH B
cor_ieb= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_ieb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$ecohet_max, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_ieb)
sd(cor_ieb)
mean(p_ieb)


#SH UQ
cor_isu= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_isu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_isu)
sd(cor_isu)
mean(p_isu)

#SH B
cor_isb= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_isb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$sochet_max, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_isb)
sd(cor_isb)
mean(p_isb)

#Trust UQ
cor_itu= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$r[1,2]
p_itu = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$unit_quality_cont, type = c("spearman"))$P[1,2]
mean(cor_itu)
sd(cor_itu)
mean(p_itu)

#Trust B
cor_itb= foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$r[1,2]
p_itb = foreach(i = 1:n.sims,.combine = 'c') %do% rcorr(dt_MI_water[[i]]$trust, dt_MI_water[[i]]$balance_cont, type = c("spearman"))$P[1,2]
mean(cor_itb)
sd(cor_itb)
mean(p_itb)




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
  bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max) + water,data = dt_MI[[i]] )
# look into
pool(fit_simsOR10)
summary(pool(fit_simsOR10))

# get AIC
AIC_simsOR10 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max) + water,data = dt_MI[[i]] ))
mean(unlist((AIC_simsOR10[[i]][2])))


#------------------- Ordinal logistic regressions: FISHERIES ----------------------------------------------------------------------------

# Make factor variable from factor variable trust 
for (i in 1:n.sims) {
  dt_MI_fishery[[i]]$trust_factor = as.factor(as.character(unlist(dt_MI_fishery[[i]]$trust)))
}

#Model 1: unit quality: sochet ecohet
fit_sims_fishOR1 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~ sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR1)
summary(pool(fit_sims_fishOR1))


# get AIC
AIC_simsfishOR1 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~ sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR1[[i]][2])))

#Model 2: balance: sochet ecohet
fit_sims_fishOR2 = 
  foreach(i = 1:n.sims) %do%
  polr(formula = balance_factor ~ sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR2)
summary(pool(fit_sims_fishOR2))

# get AIC
AIC_simsfishOR2 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(polr(formula = balance_factor ~ sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR2[[i]][2])))

#Model 3: unit quality: sochet ecohet trust 
fit_sims_fishOR3 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = unit_quality_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR3)
summary(pool(fit_sims_fishOR3))

# get AIC
AIC_simsfishOR3 = 
  foreach(i = 1:n.sims) %do%
  extractAIC(bayespolr(formula = unit_quality_factor ~as.factor(sochet_max) + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]]))
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
  polr(formula = balance_factor ~ sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR__6)
summary(pool(fit_sims_fishOR__6))

# get AIC
AIC_simsfishOR6= 
  foreach(i = 1:n.sims) %do%
  extractAIC( polr(formula = balance_factor ~sochet_max + as.factor(ecohet_max) + as.factor(trust), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR6[[i]][2])))


#Model 7: trust: sochet ecohet (with Bayespolr)
fit_sims_fishOR7 = 
  foreach(i = 1:n.sims) %do%
  bayespolr(formula = trust_factor ~ sochet_max + as.factor(ecohet_max), data = dt_MI_fishery[[i]])
#look into
pool(fit_sims_fishOR7)
summary(pool(fit_sims_fishOR7))

# get AIC
AIC_simsfishOR7= 
  foreach(i = 1:n.sims) %do%
  extractAIC( bayespolr(formula = trust_factor ~ as.factor(sochet_max) + as.factor(ecohet_max), data = dt_MI_fishery[[i]]))
mean(unlist((AIC_simsfishOR7[[i]][2])))



#------------------------- Ordinal logistic regressions: IRRIGATION ----------------------------------------------------------------------------
# Make factor variables from continuous variable trust 
for (i in 1:n.sims) {
  dt_MI_water[[i]]$trust_factor = as.factor(as.character(unlist(dt_MI_water[[i]]$trust)))
}


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


# startvalue for sochet is -0.82
# startvalues for trust are 3.32 and 4.74
# startvalues for ecohet are -0.95 and -2.14

#Model 5: unit quality: sochet ecohet trust -->  with bayespolr
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





#-----------------------------------------------------------------------
# Moderated Mediation Models 
# https://ademos.people.uic.edu/Chapter15.html
#-----------------------------------------------------------------------

# Get average (mode/mean) data over 100 dataframes

fun.mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
rowMode = function(x){apply(x,1,fun.mode)}


# averaging datasets
dt_av = data.table(array(NA,dim(dt_MI[[1]])))
colnames(dt_av) = colnames(dt_MI[[1]])
for(j in 1:dim(dt_MI[[1]])[2]) {
  if( class(dt_MI[[1]][,j]) =="numeric" & length(unique(dt_MI[[1]][,j]))>10 |  class(dt_MI[[1]][,j]) =="integer" & length(unique(dt_MI[[1]][,j]))>10
  ){
    temp = rowMeans(foreach(i = 1:100,.combine = 'cbind') %do% as.numeric(as.character(unlist(dt_MI[[i]][,j]))))
  }else{
    temp = 
      if(class(dt_MI[[1]][,j]) == "factor") {
        as.factor(rowMode(foreach(i = 1:100,.combine = 'cbind') %do% as.character(unlist(dt_MI[[i]][,j]))))
      }else{
        as.numeric(rowMode(foreach(i = 1:100,.combine = 'cbind') %do% as.character(unlist(dt_MI[[i]][,j]))))
      }
  }
  dt_av[,j] = temp
}



#### Moderated Mediation For unit quality ######
med.fit.u = lm(formula = trust ~ ecohet_max + sochet_max + water , data = dt_av)
out.fit.u = lm(formula = unit_quality_cont ~ trust + trust*water + sochet_max + ecohet_max + water, data = dt_av)

summary(med.fit.u)
summary(out.fit.u)

#water =1
med.test.ui = mediate(med.fit.u, out.fit.u, sims = 1000, treat = "ecohet_max", mediator = "trust", covariates = list(water =1), 
                      boot = FALSE)
summary(med.test.ui)
plot(med.test.ui, main = "Economic heterogeneity, trust & unit quality (Irr.)")


#water = 0 
med.test.uf = mediate(med.fit.u, out.fit.u, sims = 1000, treat = "ecohet_max", mediator = "trust", covariates = list(water =0),
                      boot = FALSE)
summary(med.test.uf)
plot(med.test.uf, main = "Economic Heterogeneity, Trust & Unit Quality (Fish.)")


#For balance
med.fit.b = lm(formula = trust ~ ecohet_max + sochet_max + water , data = dt_av)
out.fit.b = lm(formula = balance_cont ~ trust + trust*water + sochet_max + ecohet_max + water, data = dt_av)

#water =1
med.test.bi = mediate(med.fit.b, out.fit.b, sims = 1000, treat = "ecohet_max", mediator = "trust", covariates = list(water =1),
                      boot = FALSE)
summary(med.test.bi)
plot(med.test.bi, main = "Economic Heterogeneity, Trust & Balance (Irr.)")

#water = 0 
med.test.bf = mediate(med.fit.b, out.fit.b, sims = 1000, treat = "ecohet_max", mediator = "trust", covariates = list(water =0),
                      boot = FALSE)
summary(med.test.bf)
plot(med.test.bf)


#for sochet in balance

#water =1
med.test.bis = mediate(med.fit.b, out.fit.b, sims = 1000, treat = "sochet_max", mediator = "trust", covariates = list(water =1),
                       boot = FALSE)
summary(med.test.bis)
plot(med.test.bis, main = "Sociocultural Heterogeneity, Trust & Balance (Irr.)")

#water = 0 
med.test.bfs = mediate(med.fit.b, out.fit.b, sims = 1000, treat = "sochet_max", mediator = "trust", covariates = list(water =0),
                       boot = FALSE)
summary(med.test.bfs)
plot(med.test.bfs)


#sochet in Unit quality 
#water =1
med.test.uis = mediate(med.fit.u, out.fit.u, sims = 1000, treat = "sochet_max", mediator = "trust", covariates = list(water =1),
                       boot = FALSE)
summary(med.test.uis)
plot(med.test.uis)

#water = 0 
med.test.ufs = mediate(med.fit.u, out.fit.u, sims = 1000, treat = "sochet_max", mediator = "trust", covariates = list(water =0),
                       boot = FALSE)
summary(med.test.ufs)
plot(med.test.ufs)



#--------------------------------------------------------------------------------------------------------------------------
# Available case analysis
#--------------------------------------------------------------------------------------------------------------------------

#Load available case data
#Load imputed data
#Set working directory
rm(list=ls())
options(scipen=999)
setwd("directory")
load(file = "dt_temp_clean_MI.RData")

dt_available = dt_temp_clean_MI


mean(is.na(dt_available))


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



# available case max varinc by sliceid = ecohet_max
library(data.table)
dt_available = as.data.table(dt_available)[,ecohet_max:=max(as.numeric(as.character(unlist(subgroupsubvar)))),by=c("sliceid")]
dt_available = as.data.frame(dt_available)

dt_available$ecohet_factor = as.factor(as.character(unlist(dt_available$ecohet_max)))



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


#Available case fishieries
dt_available_water = as.data.frame(dt_available[which(dt_available$water == 1),])
dt_available_fishery = as.data.frame(dt_available[which(dt_available$fishery ==1),])

#--------------------------------------------------------------------------
# Available case correlations
#--------------------------------------------------------------------------


library(Hmisc)
rcorr(dt_available$ecohet_factor, dt_available$trust_factor,  type=c("spearman"))
rcorr(dt_available$sochet_factor, dt_available$trust, type=c("spearman"))
rcorr(dt_available$ecohet_factor, dt_available$unit_quality_factor, type=c("spearman"))
rcorr(dt_available$ecohet_factor, dt_available$balance_factor, type=c("spearman"))
rcorr(dt_available$sochet_factor, dt_available$unit_quality_factor, type=c("spearman"))
rcorr(dt_available$sochet_factor, dt_available$balance_factor, type=c("spearman"))
rcorr(dt_available$trust_factor, dt_available$unit_quality_factor, type=c("spearman"))
rcorr(dt_available$trust_factor, dt_available$balance_factor, type=c("spearman"))


#separate for fisheries and irrigation systems
rcorr(dt_available_water$ecohet_factor, dt_available_water$trust_factor, type=c("spearman"))
rcorr(dt_available_water$sochet_factor, dt_available_water$trust_factor, type=c("spearman"))
rcorr(dt_available_water$ecohet_factor, dt_available_water$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_water$ecohet_factor, dt_available_water$balance_factor, type=c("spearman"))
rcorr(dt_available_water$sochet_factor, dt_available_water$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_water$sochet_factor, dt_available_water$balance_factor, type=c("spearman"))
rcorr(dt_available_water$trust_factor, dt_available_water$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_water$trust_factor, dt_available_water$balance_factor, type=c("spearman"))

rcorr(dt_available_fishery$ecohet_factor, dt_available_fishery$trust_factor, type=c("spearman"))
rcorr(dt_available_fishery$sochet_factor, dt_available_fishery$trust_factor, type=c("spearman"))
rcorr(dt_available_fishery$ecohet_factor, dt_available_fishery$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_fishery$ecohet_factor, dt_available_fishery$balance_factor, type=c("spearman"))
rcorr(dt_available_fishery$sochet_factor, dt_available_fishery$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_fishery$sochet_factor, dt_available_fishery$balance_factor, type=c("spearman"))
rcorr(dt_available_fishery$trust_factor, dt_available_fishery$unit_quality_factor, type=c("spearman"))
rcorr(dt_available_fishery$trust_factor, dt_available_fishery$balance_factor, type=c("spearman"))


