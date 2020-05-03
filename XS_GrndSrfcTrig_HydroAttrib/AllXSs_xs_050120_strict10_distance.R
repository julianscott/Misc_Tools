library(tidyverse)
library(rdist)

# setwd("D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis")
# all input files must have a common XS column that identifies the XS ID
# code written to expect XS in HECRAS format. Station 0 is the far left topo point.

table1 <- read.csv("All_XSs.csv",header=TRUE)
table1 <- table1[,c(1:5)]
# each rating curve must be in ascending order and identified by XS ID
table2 <- read.csv("All_XSs_RCs.csv",header=TRUE)
table3 <- read.csv("All_sites_EP30.csv",header=TRUE)
# table of stages of the Q500 at each XS
table4 <- read.csv("All_XSs_z500.csv",header=TRUE)
df1 <- as.data.frame(table1)
df.rc <- as.data.frame(table2)
df.z500 <- as.data.frame(table4)

# remove exceptions
df.z500 <- filter(df.z500,XS != "Elkhead Creek 2_XSR106.3")
df.z500 <- filter(df.z500,XS != "Rock Creek_XS96")

# Build table of XS endpoint for Encampment (see line ~500)
df.endp <- filter(table1,Site == "Encampment")
endpts <- c(which(df.endp$hec_station == 0),c(which(df.endp$hec_station == 0)[-1] - 1,nrow(df.endp)))
endpts <- endpts[order(endpts)]

df.endp <- df.endp[endpts,]
df.endp$pt <- rep(c("start","end"),11)

# x is row number 
slope.fn <- function(x) { 
  if (x == nrow(df1)) { # if x is the last row in the table, NA
    "NA"
  } else if (df1[x + 1,"hec_station"] == df1[x,"hec_station"]) { 
    # if the hec stations are equal (vertical jump w/no horiz change), add small number and compute
    (df1[x,"hec_topo"] - df1[x + 1,"hec_topo"])/(df1[x,"hec_station"] - df1[x + 1,"hec_station"] + 0.00001)
    # if the 2nd XS does not equal the 1st XS, NA   
  } else if (df1[x + 1,"XS"] != df1[x,"XS"]) { 
    "NA"
  } else {
    # compute slope
    (df1[x,"hec_topo"] - df1[x + 1,"hec_topo"])/(df1[x,"hec_station"] - df1[x + 1,"hec_station"]) 
  }
}

slope <- sapply(1:nrow(df1),function(x) slope.fn(x)) 
slope <- as.numeric(slope)

df3 <- data.frame(df1,"slope" = slope) # intermediate df

# calculates intercept for the cum. distance & z pairs
intercept <- sapply(1:nrow(df3),function(x) df3[x,"hec_topo"]-df3[x,"slope"]*df3[x,"hec_station"])
df4 <- data.frame(df3,intercept,PtType="survey") # intermediate df
 
#Distance function 
D.fn <- function(x){
  if (x == nrow(df4)) { # if x is the last row in the table, NA
    "NA"
  } else if (df4[x + 1,"XS"] != df4[x,"XS"]) { # if the 2nd XS does not equal the 1st XS, NA  
    "NA"
  } else {
    # compute ground distance between surveyed points
    ((df4[x,"hec_station"] - df4[x + 1,"hec_station"])^2 + (df4[x,"hec_topo"] - df4[x + 1,"hec_topo"])^2)^0.5 
  }
}

# for every row in df4, use D.fn to calculate the ground distance between each point
# last point in a Xs gets a "NA"
df4$ground.d <- sapply(1:nrow(df4),function(x) D.fn(x))
df4$ground.d <- as.numeric(df4$ground.d)
head(df4)

# Get cumulative ground distance for each XS
df4 <- df4 %>%
  group_by(XS) %>%
  mutate(cum.d = cumsum(ground.d)) %>% 
  as.data.frame()

x=1
i="Encampment_XS1"
XS.i = filter(df4,XS == i)
# create final dataframe to contain all XS, all sites
df.complete <- data.frame(Site=character(0),
                          XS=character(0),
                          VegXS=character(0),
                          hec_station=numeric(0),
                          hec_topo=numeric(0),
                          slope=numeric(0),
                          intercept=numeric(0),
                          PtType=character(0),
                          ground.d=numeric(0),
                          cum.d=numeric(0),
                          int.ground.d=numeric(0),
                          int.cum.d=numeric(0))
for (i in unique(df4$XS)) {
  # define the XS.i
  XS.i = filter(df4,XS == i)
  
  # Critical! No NA values allowed in XS.i$cum.d besides the final value!
  test <- length(XS.i$cum.d) == which(is.na(XS.i$cum.d))
  if(length(test) == 1 & any(test == TRUE) == FALSE) {
    print(paste0("Script ended because XS ",i," contains NA value(s) at stations besides other than the final station"))
    break
  }
  
  # Get total ground distance for XS.i
  Di_total = sum(XS.i$ground.d,na.rm = T)
  
  # Set desired interval. Code optimized for 10 cm. Caution if other interval used.
  interval = 0.1
  
  N.full = Di_total/interval # number of 10 cm segments in the XS
  #n.complete = Di_total %/% interval # number of complete 10 cm segments in the XS
  #n.remainder = Di_total %% interval # portion of remaining 10 cm segment
  N.rows = floor(N.full) + 1 # number of rows to accomodate n.complete + n.remainder
   # set a temp dataframe to be filled interatively and dumped into df.complete
  df.built <- data.frame(Site=character(0),
                         XS=character(0),
                         VegXS=character(0),
                         hec_station=numeric(0),
                         hec_topo=numeric(0),
                         slope=numeric(0),
                         intercept=numeric(0),
                         PtType=character(0),
                         ground.d=numeric(0),
                         cum.d=numeric(0),
                         int.ground.d=numeric(0),
                         int.cum.d=numeric(0))
  
  for (x in 1:N.rows) {
    # Get XS without the last value, which is always NA (see 'test' above)
    XS.i2 <- XS.i[-which(is.na(XS.i$cum.d)),"cum.d"]
    # Find in the vector of cumulative distances for XS.i2 
    # the interval that contains decile x. 
    # findInterval()returns the index of the left boundary of the interval.
    # findInterval()+1 yields the right boundary.
    # Because vector XS.i2 is a column in XS.i, surveyPT.i 
    # corresponds to row number in XS.i. 
    surveyPT.i <- findInterval(x*interval,XS.i2) + 1
    
    # ground distance between consecutive survey points that bound decile x
    surveyPT.D <- XS.i[surveyPT.i,"ground.d"]
    
    # Initialize the first row of the df.build dataframe for x = 1 and begin trig calcs
    if (x == 1) {
      df.built <- rbind(df.built,data.frame(XS.i[1,],int.ground.d = 0.1,int.cum.d = 0.1))
    } else {
      # if the cumulative ground distance of point i (xi*interval) is within 
      # the same two SurveyPts as the previous (e.g. xi-1), then the triangles are similar
      if (findInterval(x*interval,XS.i2) + 1 == findInterval((x-1)*interval,XS.i2) + 1) {
        # ground distance (hypotenuse) of small inset triangle
        hyp <- x*interval - (XS.i[surveyPT.i,"cum.d"] - XS.i[surveyPT.i,"ground.d"])
        # cosin of large triangle formed by consecutive survey points
        surveyPT.cos <- ((XS.i[surveyPT.i,"hec_station"] - XS.i[surveyPT.i+1,"hec_station"])/XS.i[surveyPT.i,"ground.d"])
        # Calc new hec_station using similar cosine ratios for the two triangles
        h.seq = XS.i[surveyPT.i,"hec_station"] - hyp*surveyPT.cos 
        v.seq = XS.i[surveyPT.i,"slope"]*h.seq + XS.i[surveyPT.i,"intercept"]
        df.built <- rbind(df.built, data.frame(Site = XS.i[surveyPT.i,"Site"],
                                              XS = XS.i[surveyPT.i,"XS"],
                                              VegXS = XS.i[surveyPT.i,"VegXS"],
                                              hec_station = h.seq,
                                              hec_topo =  v.seq,
                                              slope = XS.i[surveyPT.i,"slope"],
                                              intercept = XS.i[surveyPT.i,"intercept"],
                                              PtType = "interpolated",
                                              ground.d = XS.i[surveyPT.i,"ground.d"],
                                              cum.d = XS.i[surveyPT.i,"cum.d"],
                                              int.ground.d = interval,
                                              int.cum.d = interval + df.built[x-1,"int.cum.d"]))
        
        # Else if the cumulative ground distance of point i (xi*interval) is 
        # not within the same two SurveyPts as the previous (e.g. xi-1)
      } else if ((findInterval(x*interval,XS.i2) + 1 > findInterval((x-1)*interval,XS.i2) + 1)) {
        # calc hyp1, the first portion of the interval, which is along the line 
        # defined by the previous SurveyPt interval
        surveyPT.i <- findInterval((x-1)*interval,XS.i2) + 1
        hyp1 <- interval - (x*interval - XS.i[surveyPT.i,"cum.d"])
        # calc hyp2, the second portion of the interval, which is along the line
        #  defined by the next SurveyPt interval
        surveyPT.i <- findInterval((x)*interval,XS.i2) + 1
        hyp2 <- interval - hyp1
        # find the hec_staion along the line defined by the next SurveyPt interval
        #  given hypotenuse and the cos ratio of the big-to-inset triangle
    
        #  cosin of large triangle formed consecutive survey points
        surveyPT.cos <- ((XS.i[surveyPT.i,"hec_station"] - XS.i[surveyPT.i+1,"hec_station"])/XS.i[surveyPT.i,"ground.d"])
        
        # Calc the new hec_station using similar cosine ratios for the two triangles
        h.seq = XS.i[surveyPT.i,"hec_station"] - hyp2*surveyPT.cos 
        v.seq = XS.i[surveyPT.i,"slope"]*h.seq + XS.i[surveyPT.i,"intercept"]
        df.built <- rbind(df.built, data.frame(Site = XS.i[surveyPT.i,"Site"],
                                               XS = XS.i[surveyPT.i,"XS"],
                                               VegXS = XS.i[surveyPT.i,"VegXS"],
                                               hec_station = h.seq,
                                               hec_topo =  v.seq,
                                               slope = XS.i[surveyPT.i,"slope"],
                                               intercept = XS.i[surveyPT.i,"intercept"],
                                               PtType = "2xinterpolated",
                                               ground.d = XS.i[surveyPT.i,"ground.d"],
                                               cum.d = XS.i[surveyPT.i,"cum.d"],
                                               int.ground.d = interval,
                                               int.cum.d = interval + df.built[x-1,"int.cum.d"]))
      }
    }
  }
  # Bundle up results
  df.complete <- rbind(df.complete,df.built)
}

tmp1 <- filter(df.complete,Site=="Encampment" & XS == "Encampment_XS8")
tmp2 <- filter(df4,Site=="Encampment" & XS == "Encampment_XS8")

ggplot()+
  geom_point(data=tmp1,aes(x=hec_station,y=hec_topo,color = PtType))+
  geom_point(data=tmp2,aes(x=hec_station,y=hec_topo,color = PtType),shape=1)+
  xlim(23,27)
  


df.final <- df.complete
df.final$oo <- 1:nrow(df.final)
# remove NA values (those where topo is NA at the final segment of each XS.)

df.final <- df.final[which(!is.na(df.final$hec_topo)),]
write.csv(df.final,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_complete.csv")

df.final <- as.data.frame(read.csv("D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_complete.csv"))
df.final <- df.final[,-1]

# identify the topo points above the Q500 stage
list.i <- list()
for (i in 1:nrow(df.final)) {
 XS.i = as.character(df.final[i,"XS"])
 z500.i <- filter(df.z500,XS == XS.i)$z500
 if (df.final[i,"hec_topo"] > z500.i) {
   list.i[i] <- i
 }
}

df.final <- df.final %>% 
  left_join(df.z500[c("XS","z500")],by = "XS") %>% 
  mutate(z500.exclude = ifelse(hec_topo > z500,"TRUE","FALSE"))


# these are the row numbers of those points with elevation greater than the stage of the Q500 for each XS
z500.exclude <- unlist(list.i)

# identify min Q in the gage record for each XS
mQ.site <- list()
for (i in unique(table3$site)) {
  # min Q for the gage
  Qmin.i = min(filter(table3,site == i)$CMS) 
  mQ.site[[i]] = Qmin.i
}

mQ.site.z500 <- sapply(1:nrow(df.z500),function(x) mQ.site[[as.character(df.z500[x,"Site"])]]) 
df.z500$minQgage <-  mQ.site.z500


# identify min Q in the rating curve for each XS
mQ.rc <- list()
for (i in unique(df.rc$XS)) {
  # min Q for the rc
  Qmin.rc = min(filter(df.rc,XS == i)$Q)
  mQ.rc[[i]] = Qmin.rc
}
mQ.rc.z500 <- sapply(1:nrow(df.z500),function(x) mQ.rc[[as.character(df.z500[x,"XS"])]])
df.z500$minQrc <-  unlist(mQ.rc.z500)

# identify Z associated with the min Q in the rating curve for each XS
mZ.rc <- list()
for (i in unique(df.rc$XS)) {
  Zmin.rc = min(filter(df.rc,XS == i)$Elev)
  mZ.rc[[i]] = Zmin.rc
}
mZ.rc.z500 <- sapply(1:nrow(df.z500),function(x) mZ.rc[[as.character(df.z500[x,"XS"])]])
df.z500$minZrc <-  mZ.rc.z500

## calc stage for each minQgage value in df.z500

# rename for work
t1.sub <- df.z500

########### For fitting power curves to z-Q rating curves for each XS in rating curve table
# empty list of power functions
power.list <- list()

# loop through each XS in rc table and fit a power curve. Store power curve in power.list
for (i in unique(df.rc$XS)) {
  df.sub = filter(df.rc,XS == i)
  y = df.sub$Elev
  x = df.sub$Q
  # a = min y; b = coef
  power.i <- tryCatch({nls(y~y0+a*I(x^b),start=list(y0=min(y),a=0.5,b=0.5),trace=T)
  }, error=function(e){cat("ERROR on i =",i,":",conditionMessage(e), "\n")})
  power.list[[i]] <- power.i
}

# for creating a plot of all the functions, with elevations standardized to 0
dev.off()
for (n in names(power.list)) {
  par(new=TRUE) 
  power.fun = power.list[[n]]
  x = seq(0,max(df.rc$Q),1)
  y = coef(power.fun)[["y0"]] + coef(power.fun)[["a"]] * (x^coef(power.fun)[["b"]])
  ynorm = y - min(y)
  if (max(ynorm) > 9) { # note that Rock Creek XS83 is an anomaly. Coef y0 is meters smaller than the smallest elevation   
    print(n)            # used to fit the model. The end result appears to be a well fitting model, but, when Q=0, 
    x = x[-1]           # the y0 coef becomes y, which creates a large gap from Q=0 to Q=1. As long as predictions
    y=y[-1]             # aren't made in that range, which there shouldn't be, b/c we only predict for Zs > minQ in record, 
    ynorm = y - min(y)  #there is no problem.
  }
  plot(x,ynorm,type="n",ylim=c(0,10),ylab = "Stage (meters above y0)",xlab = "Q (cms)",main = "Power Function Rating Curves for all XSs")
  lines(x,ynorm)
}


###########################################################

i=101
# register each Q-XS pair with the correct XS-based rating curve
Q.register <- function(i){
  # for use with following sapply command
  # where i equal row number in dataframe of plot ID, plot elevations, and plot XS
  XS.i = as.character(t1.sub[i,"XS"])
  power.i = power.list[[XS.i]] # from Final_InunQ_and_EP_for_[plots/XSs].R
  Q.i = t1.sub[i,"minQgage"]
  coef(power.i)[["y0"]] + coef(power.i)[["a"]] * (Q.i^coef(power.i)[["b"]])
}

# calculate stage for each minQgage for each row in the table
t1.sub$Z <- unlist(sapply(1:nrow(t1.sub),function(i) Q.register(i)))
df.z500$minZgage <- as.numeric(t1.sub$Z)
df.z500$final_minZ <- df.z500$minZgage

# TRUE means the topo point is less than the z500 cutoff, so it should be maintained as a valid prediction point
topo.exce.z500 <- sapply(1:nrow(df.final),function(x) df.final[x,"hec_topo"] < df.z500[df.z500$XS == as.character(df.final[x,"XS"]),]$z500)
# TRUE means the topo point is greater than the zMin cutoff, so it should be maintained as a valid prediction point
topo.less.zMin <- sapply(1:nrow(df.final),function(x) df.final[x,"hec_topo"] > df.z500[df.z500$XS == as.character(df.final[x,"XS"]),]$final_minZ)

df.final$topo.exce.z500 <- topo.exce.z500
df.final$topo.less.zMin <- topo.less.zMin

write.csv(df.final,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_complete.csv")

df.final_trimmed <- filter(df.final,topo.exce.z500 == "TRUE" & topo.less.zMin == "TRUE")
#df.final_trimmed <- filter(df.final,topo.less.zMin == "TRUE")

i="Encampment_XS1"

for (i in unique(df.final$XS)){
  #plot XS i for QAQC
  xs.plot <- ggplot() +
    geom_point(data = filter(df.final,XS == i),aes(x = hec_station,y = hec_topo),color = "black") +
    geom_point(data = filter(df.final_trimmed,XS == i),aes(x = hec_station,y = hec_topo),color = "red") +  
    labs(x = "XS Station",y = "Elevation") +
    ggtitle(i) +
    theme_bw(base_size = 20)+
    coord_fixed(ratio = 4)
  #plot(xs.plot)
  ggsave(filename=paste("XS_surfaces/",i,".jpg",sep=""),plot= xs.plot,device="jpeg",width=12, height = 4, units = "in")
}


write.csv(df.final_trimmed,"XS_strict10_ZminZ500_trimmed.csv")

###### calc inund Q for each point in table
# create a subset of the environmental matrix with just relevent columns
t1.sub <- data.frame(Site = df.final_trimmed$Site, oo = df.final_trimmed$oo,Elevation = df.final_trimmed$hec_topo,XS = df.final_trimmed$XS)
t1.sub <- df.final_trimmed
t1.sub <- df.final

sapply(1:nrow(df.final),function(x) df.final[x,"hec_topo"] > df.z500[df.z500$XS == as.character(df.final[x,"XS"]),]$final_minZ)

i=74156
# register each elevation-XS pair with the correct XS-based rating curve
Elev.register <- function(i){
  # for use with following sapply command
  # where i equal row number in dataframe of plot ID, plot elevations, and plot XS
  XS.i = as.character(t1.sub[i,"XS"])
  power.i = power.list[[XS.i]]
  Elev.i = t1.sub[i,"hec_topo"]
  # solve for InunQ.i algabraicly from the power equation 'Elev=y0+a*InunQ^b'
  InunQ.i = ((Elev.i - coef(power.i)[["y0"]]) / coef(power.i)[["a"]])^(1 / coef(power.i)[["b"]])
  #InunQ.i = ((Elev.i - 98.6869) / coef(power.i)[["a"]])^(1 / coef(power.i)[["b"]])
  #Elev.i = coef(power.i)[["y0"]] + coef(power.i)[["a"]] * (0.104776155^coef(power.i)[["b"]])
}

# calculate inundating Q for each row in the table
t1.sub$InunQ <- unlist(sapply(1:nrow(t1.sub),function(i) Elev.register(i)))


# the following points have topo too small to be calc using power curve for their XS. They are moved out of the 
# modeling surface (to be with the points with topo > z500 and <Zmin)
t1.sub <- t1.sub[-which(is.na(t1.sub$InunQ)),]

t1.sub <- left_join(t1.sub,df.z500[,c("XS","final_minZ","minQgage")],by="XS")
which(tmp$InunQ < tmp$minQgage)

write.csv(t1.sub,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_ZminZ500trim_wInunQ.csv")

#df.final <- as.data.frame(read.csv("XS_Final_ZminTrim.csv"))
df.final <- t1.sub


#####################
### Next, calculate the EP for every InunQ, based on the conditional statement:

# Given that RQ = the range of discharges in the flow duration table and
# Qi = an inundating discharge
# 
# -	If Qi is within RQ, then linearly interpolate EP between the closest two discharges in 
# RQ.
# 
# -	Else if Qi > max(RQ), then regress EP from linear model fit to the largest 10
# discharges in RQ.  
# 
# -	Else if Qi < min(RQ), then set EP to 1. 

# set source to the function file
source("D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/GuildedFloraFunction_042717.R")
source("D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/GuildedFloraFunction_080217.R")


# define Q-EP table of 30 yr daily average flows for each gage
EP_table <- data.frame(read.table("All_sites_EP30.csv", header=T, sep=","))
#EP_table <- data.frame(read.table("All_sites_EP30_S4.csv", header=T, sep=","))
#EP_table <- data.frame(read.table("Flow_Scenarios/mgmt.flows.m1.csv", header=T, sep=","))


# take a Q-EP table for QS1 scenario, and calculate EP for the modeling surface
#InunQ_table <- as.data.frame(read.csv("XS_Final_trimmed_062117.csv",header = T))
#InunQ_table <- t1.sub
#colnames(InunQ_table)[1] <- "site"
## exclude sites as desired
#site_list <- c("Encampment","Elkhead Creek 1","Lower North Brush","N. Brush","Rock Creek","Elkhead Creek 2","Hog Park Site 1","Hog Park Site 2", "Hog Park Gorge")

#InunQ_table <- filter(InunQ_table,site %in% site_list)
#InunQ_table <- InunQ_table[,c("site","XS","oo","InunQ")]

#EP_table <- Q.scenarios
colnames(EP_table)
Q.column.name <- "CMS"
#Q.column.name <- "QS7..Peak.shifted.to.April.15th"
#i="Encampment"

InunQ.df <- t1.sub

colnames(InunQ.df)[1] <- "site"

# check for NA values
which(is.na(InunQ.df$InunQ))

# check colnames to ID the proper Q column header in EP_table. Corresponds to the scenario
colnames(EP_table)
# the inundating Q column should be titled "InunQ"
colnames(InunQ.df)
site.list <- as.character(unique(t1.sub$Site))
Scenario_name = "EP_Gage"

## send EP table and InunQ tables off to the interpolation function (in GuildedFloraFunction_042717.R) 
# to perform interpolation. Returns values of 1 for Qi < min(RQ) and 0.0001 for Qi > max(RQ).
EP.interp.df <- .getEPInterp(EP_table,InunQ.df,Q.column.name)
EP.interp.df <- .getModel_Surface(EP_table = EP_table,
                    InunQ_table = InunQ.df,
                    Q.column.name = "CMS",
                    site_list = site.list,
                    Scenario_name)

##############

t1.sub$LnEP30 <- EP.interp.df

which(is.na(t1.sub$LnEP30))
max(t1.sub$LnEP30)
min(t1.sub$LnEP30)

write.csv(t1.sub,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_ZminZ500trim_wGageEP.csv")

##
# filter out MI River, Elkhead Creek 2 River R XS, and the Hog Park Sites to make a Modeling Scenario
t2.sub <- filter(t1.sub,Site %in% c("Encampment","Elkhead Creek 1","Lower North Brush","N. Brush","Rock Creek","Elkhead Creek 2"))
t2.sub <- filter(t2.sub,!grepl('Elkhead Creek 2_XSR', XS))
write.csv(t2.sub,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_ZminZ500trim_wGageEP_MS.csv")

# filter out all but the Hog Park Sites to make a reservoir Modeling Scenario
t2.sub <-  filter(t1.sub,grepl('Hog', Site))
write.csv(t2.sub,"D:/EnvironmentalFlowsProject/Medicine_Bow_Route_NF/Veg_Data/Analysis/_Final_Modeling_Surface_Files/XS_strict10_ZminZ500trim_wGageEP_HogMS.csv")


# add local Northing and Easting data to Encampment subset

Enc.trimmed.topo <- as.data.frame(read.csv("")) #XS_strict10_ZminZ500_trimmed.csv
Enc.trimmed.topo <- filter(left_join(XS.sub4,XS.df[,c("oo","hec_station")],by = "oo"),site == "Encampment")
Enc.trimmed.topo$XS <- droplevels(Enc.trimmed.topo$XS)
Enc.trimmed.topo$E <- 0
Enc.trimmed.topo$N <- 0

df.endp$XS <- droplevels(df.endp$XS)

for (i in 1:nrow(Enc.trimmed.topo)){
  XS.i = Enc.trimmed.topo[i,"XS"]
  x1 = filter(df.endp,XS == XS.i & pt == "start")$Local_x
  x2 = filter(df.endp,XS == XS.i & pt == "end")$Local_x
  y1 = filter(df.endp,XS == XS.i & pt == "start")$Local_y
  y2 = filter(df.endp,XS == XS.i & pt == "end")$Local_y
  A = x2-x1
  B = y2-y1
  C =  (A^2 + B^2)^0.5
  c = Enc.trimmed.topo[i,"hec_station"]
  a = c*(A/C)
  b = c*(B/C)
  Enc.trimmed.topo[i,"E"] = x1+a
  Enc.trimmed.topo[i,"N"] = y1+b
}

Enc.complete.topo <- as.data.frame(read.csv("XS_strict10_complete.csv"))
#Enc.complete.topo <- XS.df
Enc.complete.topo <- filter(Enc.complete.topo,Site == "Encampment")
Enc.complete.topo$XS <- droplevels(Enc.complete.topo$XS)
Enc.complete.topo <- Enc.complete.topo[which(!Enc.complete.topo$oo %in% tmp$oo),]
Enc.complete.topo$E <- 0
Enc.complete.topo$N <- 0
#Enc.complete.topo$upl_lol <- "NA"

for (i in 1:nrow(Enc.complete.topo)){
  XS.i = Enc.complete.topo[i,"XS"]
  x1 = filter(df.endp,XS == XS.i & pt == "start")$Local_x
  x2 = filter(df.endp,XS == XS.i & pt == "end")$Local_x
  y1 = filter(df.endp,XS == XS.i & pt == "start")$Local_y
  y2 = filter(df.endp,XS == XS.i & pt == "end")$Local_y
  A = x2-x1
  B = y2-y1
  C =  (A^2 + B^2)^0.5
  c = Enc.complete.topo[i,"hec_station"]
  a = c*(A/C)
  b = c*(B/C)
  Enc.complete.topo[i,"E"] = x1+a
  Enc.complete.topo[i,"N"] = y1+b
  Enc.complete.topo[i,"upl_lol"] = ifelse(Enc.complete.topo[i,"hec_topo"] < min(filter(Enc.trimmed.topo,XS == XS.i)$hec_topo),"low","up")
}

ggplot(filter(df.final,Site == "Encampment"),aes(x=hec_station, y=hec_topo))+
  geom_point() +
  facet_wrap(~XS)

plot(x=Enc.trimmed.topo$E,y=Enc.trimmed.topo$N)


# function to sum the meters (from the Distance column in df.final) of ground surface in EP5 decile bins
# if the InunQ is > the max flow in the record, given EP5 of 0. If InunQ < min flow, given a 1. 

df.final.ex <- XS.prob




d.total <- sum(df.final.ex$Distance)
myfun <- function(x) sum(subset(df.final.ex,df.final.ex$EP5>x & df.final.ex$EP5<=x+0.1)$Distance/d.total*100)
r1 <- seq(from=0.000001,to=1,0.1)
EP5.bin <- sapply(r1, myfun)

tp <- ggplot()+
  geom_bar(aes(x=factor(round(log10(r1),1)),y=EP5.bin),stat="identity")+
  labs(x="Exceedence Probability (EP5)",y="Available Habitat (%)")+
  ggtitle("All XS")+
  theme_grey(base_size=20)
plot(tp)


ggplot(data=data.frame(EP5=EP5.bin,site="Encampment"),aes(x=site,y=EP5)) +
  geom_jitter(alpha = .9) +
  geom_violin(alpha = .50)
 



cross_sections <- filter(XS.df,site == "Encampment")
distance_var <- "ground.d" 
summary_var <- "EP30"
bin_by <- 0.1
bin_low <- 0
bin_hi <- 1
d_norm <- 1#sum(cross_sections[[distance_var]])




write.csv(df.final.ex,file="Encampment_XS_EP5_S2Q.csv")
write.csv(df.final,file="Encampment_XS_complete_topo.csv")


for (i in unique(df.final.ex$XS)){
  
  test <- subset(df.final.ex,XS==i)
  myfun <- function(x) sum(subset(test,EP5>x & EP5<=x+0.1)$Distance)
  r1 <- seq(from=0,to=1,0.1)
  EP5.bin <- sapply(r1, myfun)
  tp <- ggplot()+
    geom_bar(aes(x=factor(r1),y=EP5.bin),stat="identity")+
    labs(x="Exceedence Probability (EP5)",y="Meters")+
    ggtitle(unique(test$XS))+
    theme_grey(base_size=20)
  plot(tp)
}

ggplot()+
  geom_bar(aes(x=factor(r1),y=EP5.bin),stat="identity")+
  labs(x="Exceedence Probability (EP5)",y="Meters")+
  theme_grey(base_size=20)












