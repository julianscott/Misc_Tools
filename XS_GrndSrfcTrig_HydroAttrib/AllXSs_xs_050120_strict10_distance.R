library(tidyverse)
library(rdist)

# all input files must have a common XS column that identifies the XS ID
# code written to expect XS in HECRAS format. Station 0 is the far left topo point.
# See csv files for format and expected headers

table1 <- read.csv("All_XSs.csv",header=TRUE)
df1 <- as.data.frame(table1)

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

# Option to set interators for testing
# x <- 1
# x <- 28
# i <- "Hog Park Gorge_XS132.3"
# XS.i <- filter(df4,XS == i)

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
# Loop through each XS and slice into 10 cm sections along the ground surface
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
  # Get XS without the last value, which is always NA (see 'test' above)
  XS.i2 <- XS.i[-which(is.na(XS.i$cum.d)),"cum.d"]
  # N.rows
  for (x in 1:N.rows) {
    # Find in the vector of cumulative distances for XS.i2 
    # the interval that contains decile x. 
    # findInterval()returns the index of the left boundary of the interval.
    # findInterval()+1 yields the right boundary.
    # Because vector XS.i2 is a column in XS.i, surveyPT.i 
    # corresponds to row number in XS.i. 
    surveyPT.i <- findInterval(x*interval,XS.i2) + 1
    
    # ground distance between consecutive survey points that bound decile x
    surveyPT.D <- XS.i[surveyPT.i,"ground.d"]
    
    # if the cumulative ground distance of point i (xi*interval) is within 
    # the same two SurveyPts as the previous (e.g. xi-1), then the triangles are similar
    if (findInterval(x*interval,XS.i2) + 1 == findInterval((x-1)*interval,XS.i2) + 1) {
      # ground distance (hypotenuse) of small inset triangle
      hyp <- x*interval - (XS.i[surveyPT.i,"cum.d"] - XS.i[surveyPT.i,"ground.d"])
      # hyp <- interval
      # cosin of large triangle formed by consecutive survey points
      surveyPT.cos <- ((XS.i[surveyPT.i,"hec_station"] - XS.i[surveyPT.i+1,"hec_station"])/XS.i[surveyPT.i,"ground.d"])
      # Calc new hec_station using similar cosine ratios for the two triangles
      h.seq <- XS.i[surveyPT.i,"hec_station"] - hyp*surveyPT.cos 
      v.seq <- XS.i[surveyPT.i,"slope"]*h.seq + XS.i[surveyPT.i,"intercept"]
      
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
                                            int.cum.d = ifelse(x == 1,
                                                                interval,
                                                                interval + df.built[x-1,]$int.cum.d)))
      
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
  # Bundle up results
  df.complete <- rbind(df.complete,df.built)
  print(paste0("Cross section: ",i, "complete."))
}

# Plot one cross section
tmp1 <- filter(df.complete,XS == unique(df4$XS)[1])
tmp2 <- filter(df4, XS == unique(df4$XS)[1])


tmp1 <- filter(df.complete,XS == "Elkhead Creek 1_XS21.9")
tmp2 <- filter(df4, XS == "Elkhead Creek 1_XS21.9")

library(plotly)

ggplot()+
  geom_point(data=tmp1,aes(x=hec_station,y=hec_topo,color = PtType),size = 1)+
  geom_point(data=tmp2,aes(x=hec_station,y=hec_topo,color = PtType),shape=1) +
  coord_equal()

ggplotly()

df.final <- df.complete
df.final$oo <- 1:nrow(df.final)

# remove NA values (those where topo is NA at the final segment of each XS.)
df.final <- df.final[which(!is.na(df.final$hec_topo)),]
# write.csv(df.final,"XS_strict10_complete.csv")







