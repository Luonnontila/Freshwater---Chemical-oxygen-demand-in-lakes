### Data obtaining functions

get_VESLA_sampling_sites <- function(x) {
  
Q1 <- sapply(x$V_Vesimuodostuma_Id, function(x) paste0("http://rajapinnat.ymparisto.fi/api/vesla/2.0/odata/Paikka?$filter=V_Vesimuodostuma_Id%20eq%20",
                                    x,
                                    "&$select=Paikka_Id,Nimi,J_Jarvi_Id,V_Vesimuodostuma_Id,H_Kunta_Id,KoordErTmIta,KoordErTmPohj,Koordsto,Lisatieto,Tarkkuus"))
Q1 <- lapply(Q1, function(x) httr::GET(paste0(x)))
Q1 <- lapply(Q1, function(x) jsonlite::fromJSON(httr::content(x, as="text"))[[2]])
Q1 <- dplyr::bind_rows(Q1)
#
Q1 <- dplyr::left_join(Q1, x[,c("Jarvi_nimi","Jarvi_tyyppi","Jarvi_tyyppi3","Maant_jako","V_Vesimuodostuma_Id"
)], by=c("V_Vesimuodostuma_Id"))

return(Q1)
}


get_VESLA_all_data <- function(x) {
  # Create queries for each sampling site
  Q1 <- sapply(x, function(x) paste0("https://rajapinnat.ymparisto.fi/api/vesla/2.0/odata/VedenlNayte?$filter=Naytteenotto/Paikka_Id%20eq%20",
                                     x,
                                     "&$expand=VedenlTulos,Naytteenotto"))
  
  # Obtain the data
  Q1 <- lapply(Q1, function(x) httr::GET(paste0(x)))
  
  # Transform the JSON data into list objects
  Q1 <- lapply(Q1, function(x) jsonlite::fromJSON(httr::content(x, as="text")))
  
  # Withdraw only the data necessary
  Q1 <- lapply(Q1, function(x) {x[[2]]})
  
  return(Q1)
}


### Shaping functions

shape_VESLA_data_events <- function(x,t) {
  S1 <- lapply(x, function (x) {x[2]$Naytteenotto[,c(1,2,3,4,6)]})
  #S2 <- lapply(x, function (x) {x[5]$SyvyysYla})
  #S3 <- lapply(x, function (x) {x[6]$SyvyysAla})
  S1 <- dplyr::bind_rows(S1)
  #S2 <- unlist(S2)
  #S3 <- unlist(S3)
  #S1 <- cbind(S1,SyvyysYla=S2,SyvyysAla=S3)
  
  S1$Aika.e <- as.POSIXct(S1$Aika, tz = "", "%Y-%m-%dT%H:%M:%OS")
  S1 <- transform(S1, npvm = as.numeric(Aika.e),
                  na  = as.numeric(format(Aika.e, '%Y')),
                  nkk = as.numeric(format(Aika.e, '%m')),
                  pvma    = as.numeric(format(Aika.e, '%j')))
  
  t <- as.POSIXct(t, tz="","%Y-%m-%d")
  
  S1 <- S1[which(S1$Aika.e>t),]
  
  S1 <- unique(S1)
  
  return(S1)
}


shape_VESLA_data_samples <- function(x,y) {
  S1 <- lapply(x, function (x) {x[3:8]})
  S1 <- dplyr::bind_rows(S1)
  S1$SyvyysYla <- as.numeric(S1$SyvyysYla)
  S1$SyvyysAla <- as.numeric(S1$SyvyysAla)
  S1 <- S1[which(S1$SyvyysYla<=max(y)),]
  S1 <- S1[which(S1$SyvyysYla>=min(y)),]
  
  return(S1)
}


shape_VESLA_data_results <- function(x,y) {
  
  S1 <- lapply(x, function (x) {dplyr::bind_rows(x[1]$VedenlTulos)})
  S1 <- dplyr::bind_rows(S1)
  S1 <- dplyr::filter(S1, is.na(Lippu_Id))
  S1 <- S1[which(S1$Maaritys_Id %in% y),]
  
  return(S1)
}


shape_VESLA_data_combine <- function(x,y,z,d) {
  
  D1 <- dplyr::inner_join(x, y, by="VedenlNayte_Id")
  D1 <- dplyr::inner_join(D1, z, by="Naytteenotto_Id")
  D1 <- dplyr::inner_join(D1, d, by="Paikka_Id")
  D1 <- D1[which(D1$Arvo!=0),]
  #D1 <- D1[which(D1$SyvyysYla<=max(g)),]
  #D1 <- D1[which(D1$SyvyysYla>=min(g)),]
  
  D1$Paikka_Id <- as.factor(D1$Paikka_Id)
  D1$V_Vesimuodostuma_Id <- as.factor(D1$V_Vesimuodostuma_Id)
  D1$J_Jarvi_Id <- as.factor(D1$J_Jarvi_Id)
  D1$Jarvi_tyyppi <- as.factor(D1$Jarvi_tyyppi)
  
  #D1$Aika.e <- as.POSIXct(D1$Aika, tz = "", "%Y-%m-%dT%H:%M:%OS")
  #D1 <- transform(D1, npvm = as.numeric(Aika.e),
  #                na  = as.numeric(format(Aika.e, '%Y')),
  #                nkk = as.numeric(format(Aika.e, '%m')),
  #                pvma    = as.numeric(format(Aika.e, '%j')))
  
  #t <- as.POSIXct(t, tz="","%Y-%m-%d")
  
  #D1 <- D1[which(D1$Aika.e>t),]
  
  return(D1)
  
}


shape_VESLA_data_comparison <- function(x,y,z) {
  
  # x = obtained query object
  # y = first time window (c(start.year, end.year))
  # z = second time window (c(star.year, end.year))
  
  xy <- x[which(x$na>=min(y)&x$na<=max(y)),]
  xy$f_na <- as.factor(paste0("f",xy$na))
  xy$period <- "beginning"
  
  xz <- x[which(x$na>=min(z)&x$na<=max(z)),]
  xz$f_na <- as.factor(paste0("f",xz$na))
  xz$period <- "end"
  
  xe <- rbind(xy,xz)
  
  if(nrow(xe)>0){
    
    return(xe)}
  
  else{
    print("No observations")
  }
}

shape_VESLA_data_time_region <- function(x,y) {
  
  # x = data list
  # y = time window (c(start.year, end.year))
  
  
  xe <- x[which(x$na>=min(y)&x$na<=max(y)),]
  xe$f_na <- as.factor(paste0("f",xe$na))
  
  if(nrow(xe)>0){
    
    return(xe)}
  
  else{
    print("No observations")
  }
  
}


### Summarizing functions

summarize_VESLA_lake_info <- function(x) {
  
  D1  <- x %>% 
    group_by(Jarvi_tyyppi3, Maant_jako, Jarvi_nimi, J_Jarvi_Id, V_Vesimuodostuma_Id) %>% 
    summarize(n_samples = n(),
              n_sites = length(unique(Paikka_Id)),
              n_waterbodies = length(unique(V_Vesimuodostuma_Id)),
              n_years = length(unique(na)),
              started = min(Aika.e),
              ended = max(Aika.e))
  
  return(D1)
}


summarize_VESLA_data <- function(x,y) {
  
  if(y=="J_Jarvi_Id"){
  
  D1 <- x %>% 
    group_by(Jarvi_tyyppi, Jarvi_tyyppi3, Maant_jako, Jarvi_nimi, na, J_Jarvi_Id) %>% 
    summarize(median_CODmn = median(Arvo, na.rm=TRUE),
              mean_CODmn = mean(Arvo, na.rm=TRUE),
              n_CODmn = n(),
              sd_CODmn = sd(Arvo, na.rm=TRUE))
  return(D1)
  }
  
  if(y=="Paikka_Id"){
  
  D1 <- x %>% 
      group_by(Jarvi_tyyppi, Jarvi_tyyppi3, Maant_jako, Jarvi_nimi, na, J_Jarvi_Id, V_Vesimuodostuma_Id, Paikka_Id) %>% 
      summarize(median_CODmn = median(Arvo, na.rm=TRUE),
                mean_CODmn = mean(Arvo, na.rm=TRUE),
                n_CODmn = n(),
                sd_CODmn = sd(Arvo, na.rm=TRUE))
  return(D1)
  }
  
  else {
    print("Summarizing level is neither 'J_Jarvi_Id' nor 'Paikka_Id'")
  }
  
  
  
}



### Analysis functions

analyze_VESLA_minimum_MA <- function(x,y) {
  # x = data used in the analyses
  # y = time window size
  
  # Calculate moving averages for each lake
  moving_ave <- x %>% group_by(J_Jarvi_Id) %>% mutate(ma5=zoo::rollapply(median_CODmn, y, mean, align="center", fill=NA))
  
  # Obtain minimum moving averages for each lake
  mamins <- moving_ave %>% group_by(J_Jarvi_Id) %>% summarize(min_MA5 = min(ma5, na.rm=TRUE))
  mamins <- mamins[mamins$min_MA5!=Inf,]
  mamins <- mamins[!is.na(mamins$J_Jarvi_Id),]
  
  # Shape into a list with J_Jarvi_Id as the grouping
  shapelist <- lapply(mamins$J_Jarvi_Id, function(x){
    moving_ave[which(moving_ave$J_Jarvi_Id==x),]})
  
  # Select years corresponding to minimums
  maminyearslist <- lapply(seq(from=1, to=length(shapelist)), function(x){
    shapelist[[x]][which(shapelist[[x]]$ma5==mamins$min_MA5[x]), c("Jarvi_tyyppi3", "Maant_jako", "Jarvi_nimi","na","J_Jarvi_Id")]
  })
  
  # Shape list to data.frame
  maminyearsdf <- data.table::rbindlist(maminyearslist)
  
  uniqyears <- unique(maminyearsdf$na)
  yearsmode <- uniqyears[which.max(tabulate(match(maminyearsdf$na, uniqyears)))]
  
  timeframe <- c((yearsmode-floor(y/2)), (yearsmode+floor(y/2)))
  
  return(timeframe)
  
  }
  


analyze_VESLA_status <- function(x){
  # x = shaped indicator data
  if(any(colnames(x)=="Paikka_Id")){
    print("Lowest level of grouping: sampling site")
    
  M1 <- brm(bf(median_CODmn ~ 0 + period + Jarvi_tyyppi3 +
                 (1|Jarvi_tyyppi3:J_Jarvi_Id) +
                 (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id) +
                 (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id:Paikka_Id)),
            data=x,
            family = lognormal(),
            prior = prior(normal(0,5), class="b"),
            cores = 4,
            chains = 4,
            #init = 0.1,
            #seed = 123,
            sample_prior = TRUE,
            iter = 10000, warmup = 4000, thin = 5,
            control = list(adapt_delta = 0.999, max_treedepth = 11))
  
  
  return(M1)
  
  } else {
    
    print("Lowest level of grouping: lake")
    M1 <-  brm(bf(median_CODmn ~ period + Jarvi_tyyppi3
                    (1|Jarvi_tyyppi3:J_Jarvi_Id) + 
                    (1|f_na)),
               data=x,
               family = lognormal(),
               cores = 4, 
               seed = 123,
               iter = 5000, warmup = 1000, thin = 5,
               control = list(adapt_delta = 0.99,
                              max_treedepth = 12))
    
    return(M1)}
  
  
  }


analyze_VESLA_trend <- function(x) {
  if(any(colnames(x)=="Paikka_Id")){
    print("Lowest level of grouping: sampling site")
    M1 <- brm(bf(median_CODmn ~ na + Jarvi_tyyppi3 +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id) +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id) +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id:Paikka_Id)),
              prior =  prior(normal(0,5), class = b),
              data=x,
              family = lognormal(),
              cores = 4, 
              seed = 123,
              iter = 6000, warmup = 1000, thin = 5,
              control = list(adapt_delta = 0.999, max_treedepth=11))
    
  }
  
  else{
    print("Lowest level of grouping: lake")
    M1 <-  brm(bf(log(median_CODmn) ~ na + JArvi_tyyppi3 +
                    (1|Jarvi_tyyppi3:J_Jarvi_Id) + 
                    (1|f_na)),
               data=x,
               family = lognormal(),
               cores = 4, 
               seed = 123,
               iter = 5000, warmup = 1000, thin = 5,
               control = list(adapt_delta = 0.99,
                              max_treedepth = 12))
    
    return(M1)}
}

analyze_VESLA_GAM <- function(x) {
  
  if(any(colnames(x)=="Paikka_Id")){
    print("Lowest level of grouping: sampling site")
    M1 <- brm(bf(median_CODmn ~ Jarvi_tyyppi3 +
                   s(na, by=Jarvi_tyyppi3) +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id) +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id) +
                   (1|Jarvi_tyyppi3:J_Jarvi_Id:V_Vesimuodostuma_Id:Paikka_Id) +
                   (1|na)),
              data=x,
              family = lognormal(),
              cores = 4, 
              seed = 123,
              iter = 6000, warmup = 1000, thin = 5,
              control = list(adapt_delta = 0.99, max_treedepth=11))
    
  }
  
  else{
    print("Lowest level of grouping: lake")
    M1 <-  brm(bf(log(median_CODmn) ~ Jarvi_tyyppi + 
                  s(na, by=Jarvi_tyyppi) +
                  (1|Jarvi_tyyppi:J_Jarvi_Id) + 
                    (1|na)),
             data=x,
             family = gaussian(),
             cores = 4, 
             seed = 123,
             iter = 5000, warmup = 1000, thin = 5,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 12))
  
  return(M1)}
  
}


### Output functions

output_table_full_VESLA_GAM <- function(x) {
  ot1 <- conditional_effects(x)
  ot1df <- as.data.frame(ot1$`na:Jarvi_tyyppi`)
  
  write.csv(ot1df, file="Output/COD_GAM_fit.csv", row.names=FALSE)
  
  return(ot1df)
}

output_table_VESLA_status_comparison <- function(x,y){
  
  # x = brmsfit object1 produced with the function analyze_VESLA_indicator_status
  # y = data1 used in the fitting
  
  # Create a set of null hypotheses to which compare against. 
  # The first hypothesis tests for whether a difference between the compared time periods exists at all. 
  # The ones that follow after the first hypothesis are for specific defined status classes. 
  hyp <- c("periodend=periodbeginning", # 1 # This is for testing if there is a difference between the time periods (and obtaining an estimate for that)
           "(exp(periodend)/exp(periodbeginning))-1=0", # 2
           "exp(periodend)/exp(periodbeginning)>1.5", # 3  # This is the null hypothesis for state being "very bad" (i.e. x > 1.5)
           # The opposite H1 is thus "less than very bad" (i.e. x < 1.5)
           "exp(periodend)/exp(periodbeginning)>1.2", # 4 # This is the null hypothesis of state being "bad or worse" (i.e. x > 1.2)
           # The opposite H1 is thus "less than bad" (i.e. x < 1.2)
           "exp(periodend)/exp(periodbeginning)>1.05", # 5 # This is the null hypothesis of state being "satisfactory or worse" (i.e. x > 1.05)
           # The opposite H1 is thus "less than satisfactory" (i.e. x < 1.05)
           "exp(periodend)/exp(periodbeginning)>0.80", # 6 # This is the null hypothesis of state being "good or worse" (i.e. x > 0.8)
           # The opposite H1 is thus "very good" (i.e. x < 0.8)
           "exp(periodend)/exp(periodbeginning)<0.80") # 7 # This is the null hypothesis of state being "very good" (i.e. x > 0.5)
           # The opposite H1 is thus "worse than very good" (i.e. x < 0.5)
  
  # Perform hypotehsis testing using brms
  h95 <- hypothesis(x, hyp, alpha=0.05)
  h90 <- hypothesis(x, hyp, alpha=0.1)
  
  
  # Transform null hypothesis tests to BF10 Bayes factors
  bf10 <- h95$hypothesis$Evid.Ratio[3:7]
  
  # Obtain difference estimate
  pe <- h95$hypothesis$Estimate[2]*100
  ci95 <- c(h95$hypothesis$CI.Lower[2]*100,h95$hypothesis$CI.Upper[2]*100)
  ci90 <- c(h90$hypothesis$CI.Lower[2]*100,h90$hypothesis$CI.Upper[2]*100)
  
  # Obtain error estimates
  err <- round(h95$hypothesis$Est.Error[2]*100,1)
  err.prec <- round((h95$hypothesis$Est.Error[2]/abs(h95$hypothesis$Estimate[2]))*100,1)
  
  # Posterior draws
  status_samples <- ((exp(as_draws_df(x)$b_periodend)/exp(as_draws_df(x)$b_periodbeginning))-1)*100
  n_samples <- length(status_samples)
  
  samples_very_bad <- status_samples[which(status_samples>50)]
  prec_very_bad <- length(samples_very_bad)/n_samples * 100
  
  samples_bad <- status_samples[which(status_samples<50)]
  samples_bad <- samples_bad[which(samples_bad>20)]
  prec_bad <- length(samples_bad)/n_samples * 100
  
  samples_satisfactory <- status_samples[which(status_samples<20)]
  samples_satisfactory <- samples_satisfactory[which(samples_satisfactory>5)]
  prec_satisfactory <- length(samples_satisfactory)/n_samples *100
  
  samples_good <- status_samples[which(status_samples<5)]
  samples_good <- samples_good[which(samples_good>(-20))]
  prec_good <- length(samples_good)/n_samples * 100
  
  samples_excellent <- status_samples[which(status_samples<(-20))]
  prec_excellent <- length(samples_excellent)/n_samples * 100
  
  
  qL <- quantile(((exp(as_draws_df(x)$b_periodend)/exp(as_draws_df(x)$b_periodbeginning))-1)*100, probs=0.10)
  qR <- quantile(((exp(as_draws_df(x)$b_periodend)/exp(as_draws_df(x)$b_periodbeginning))-1)*100, probs=0.90)
  
  # If the credible interval does not include zero, then proceed to checking minimum decline or minimum rise
  if((qL*qR)>0){
    # If the difference seems to be negative (a change to worse), then obtain the minimum decline (with 90% certainty)
    # If the difference seems to be positive (a change to better), then obtain the minimum rise (with 90% certainty)
    if(qR<0){
      e90 <- qR
    } else{
      e90 <- qL
    }
    # If the 80% credible interval includes zero then qL*qR < 0. 
    # Then we will proceed to checking the minimum decline (worst case scenario; with 90% certainty).
  }else{
    e90 <- "no evidence of change"
  }
  
  
  
  status <-data.frame(status_class = c("excellent", "good", "satisfacory", "poor", "very poor"),
                      status_rank = c(1,2,3,4,5),
                      prec_posterior = c(prec_excellent, prec_good, prec_satisfactory, prec_bad, prec_very_bad),
                      bf10 = rev(bf10))
  
  
  # Type in interpretations of different BF10 Bayes factors
  conf_class <- data.frame(class=c("extreme supporting evidence",
                                   "very strong supporting evidence", 
                                   "strong supporting evidence", 
                                   "moderate supporting evidence", 
                                   "anecdotal countering evidence",
                                   "no evidence",
                                   "anecdotal countering evidence",
                                   "moderate countering evidence",
                                   "strong countering evidence",
                                   "very strong countering evidence",
                                   "extreme countering evidence"),
                           bf10_u=c(Inf,100,30,10,3,1,1,1/3,0.1,1/30,0.01),
                           bf10_l=c(100,30,10,3,1,1,1/3,0.1,1/30,0.01,-Inf))
  
  # Withdraw how much support each status class gets based on BF10 Bayes factors
  status$support <- sapply(seq(1:nrow(status)), 
                           function(i){
                             conf_class$class[which(conf_class$bf10_l<=status$bf10[i]&conf_class$bf10_u>=status$bf10[i])]
                           })
  
  if(length(grep("extreme supporting evidence",status$support))>0){
    selected_status <- status[grep("extreme supporting evidence",status$support),]
    selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
  }else{
    if(length(grep("very strong supporting evidence",status$support))>0){
      selected_status <- status[grep("very strong supporting evidence",status$support),]
      selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
    }else{
      if(length(grep("strong supporting evidence",status$support))>0){
        selected_status <- status[grep("strong supporting evidence",status$support),]
        selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
      }else{
        if(length(grep("moderate supporting evidence",status$support))>0){
          selected_status <- status[grep("moderate supporting evidence",status$support),]
          selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
        }else{
          selected_status <- data.frame(status_class="undetermined", support="no evidence")
        }
      }
    }
  }
  
  
  
  
  state <- data.frame(indicator = paste0("VESLA_COD"),
                      reference_time_range = paste0(min(y$na[which(y$period=="beginning")]),
                                                    "-",
                                                    max(y$na[which(y$period=="beginning")])),
                      comparison_time_range = paste0(min(y$na[which(y$period=="end")]),
                                                     "-",
                                                     max(y$na[which(y$period=="end")])),
                      difference.cred.I.95 = paste0("[",round(ci95[1],3),";",round(ci95[2],3),"]"),
                      difference.cred.I.90 = paste0("[",round(ci90[1],3),";",round(ci90[2],3),"]"),
                      difference.point.est = pe,
                      difference.error = err,
                      difference.error.prec = err.prec,
                      estimate_90 = e90,
                      prec.excellent = prec_excellent,
                      prec.good = prec_good,
                      prec.satisfactory = prec_satisfactory,
                      prec.poor = prec_bad,
                      prec.very.poor = prec_very_bad,
                      status = selected_status$status_class,
                      status_support = selected_status$support)
  
  return(state)
}

output_table_VESLA_trend <- function(x,y) {

# x = fitted model
# y = data used in the fitting


#ye <- cbind(y, sd_log = (log(y$upper)-log(y$mean)))
ye <- y

ot <- brms::fixef(x, probs=c(0.025,0.05,0.95,0.975))
ot <- as.data.frame(ot)
ot <- ot[!row.names(ot) %in% "Intercept",]
ot <- ot[row.names(ot) %in% "na",]

hyp <- c("na=0", # This is for testing if there is a trend (and obtaining an estimate for that)
         "na<(-0.05)",  # This is the hypothesis of trend being "rapidly improving" (i.e. x < -0.05)
         # The opposite H0 is thus "less than rapidly improving" (i.e. x > -0.05)
         "na<0", # This is the hypothesis of trend being "at least improving" (i.e. x < 0)
         # The opposite H0 is thus trend being "less than improving" (i.e. decreasing or stable; x < 0)
         "na>0", # This is the null hypothesis of trend being "at least decreasing" (i.e. x > -0.05)
         # The opposite H0 is thus trend being "rapidly decreasing" (i.e. x < -0.05; see below)
         "na>0.05") # This is the hypothesis of state being "rapidly decreasing" (i.e. x < 0.5), i.e. state being "very poor"



# Perform hypotehsis testing using brms
h <- hypothesis(x, hyp)

# Transform null hypothesis tests to BF10 Bayes factors
bf10 <- h$hypothesis$Evid.Ratio[2:5]

status <-data.frame(status_class = c("rapdily improving", "minimum improving", "minimum decreasing", "rapidly decreasing"),
                    status_rank = c(1,2,3,4),
                    bf10 = bf10)

# Type in interpretations of different BF10 Bayes factors
conf_class <- data.frame(class=c("extreme supporting evidence",
                                 "very strong supporting evidence", 
                                 "strong supporting evidence", 
                                 "moderate supporting evidence", 
                                 "anecdotal countering evidence",
                                 "no evidence",
                                 "anecdotal countering evidence",
                                 "moderate countering evidence",
                                 "strong countering evidence",
                                 "very strong countering evidence",
                                 "extreme countering evidence"),
                         bf10_u=c(Inf,100,30,10,3,1,1,1/3,0.1,1/30,0.01),
                         bf10_l=c(100,30,10,3,1,1,1/3,0.1,1/30,0.01,-Inf))

# Withdraw how much support each status class gets based on BF10 Bayes factors
status$support <- sapply(seq(1:nrow(status)), 
                         function(i){
                           conf_class$class[which(conf_class$bf10_l<=status$bf10[i]&conf_class$bf10_u>=status$bf10[i])]
                         })

if(length(grep("extreme supporting evidence",status$support))>0){
  selected_status <- status[grep("extreme supporting evidence",status$support),]
  selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
}else{
  if(length(grep("very strong supporting evidence",status$support))>0){
    selected_status <- status[grep("very strong supporting evidence",status$support),]
    selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
  }else{
    if(length(grep("strong supporting evidence",status$support))>0){
      selected_status <- status[grep("strong supporting evidence",status$support),]
      selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
    }else{
      if(length(grep("moderate supporting evidence",status$support))>0){
        selected_status <- status[grep("strong supporting evidence",status$support),]
        selected_status <- selected_status[which(selected_status$status_rank == min(selected_status$status_rank)),]
      }else{
        selected_status <- data.frame(status_class="undetermined", support="no evidence")
      }
    }
  }
}


if(ot$Q2.5[1]*ot$Q97.5[1]>0){
  # If the difference seems to be negative (a change to worse), then obtain the minimum decline (with 90% certainty)
  # If the difference seems to be positive (a change to better), then obtain the minimum rise (with 90% certainty)
  if(ot$Q97.5[1]<0){
    e90 <- quantile((as_draws_df(x)$b_na), probs=0.9)
  } else{
    e90 <- quantile((as_draws_df(x)$b_na), probs=0.1)
  }
  # If the credible interval includes zero then ci95[1]*ci95[2]<0. 
  # Then we will proceed to checking the minimum decline (worst case scenario; with 90% certainty).
}else{
  e90 <- quantile((as_draws_df(x)$b_na), probs=0.9)
}


ote <- as.data.frame(cbind(trend = ot$Estimate, 
                           trend.se = ot$Est.Error,
                           trend.prob90 = e90,
                           trend.lower95 = ot$Q2.5,
                           trend.upper95 = ot$Q97.5,
                           trend.lower90 = ot$Q5,
                           trend.upper90 = ot$Q95))

ote <- dplyr::mutate_at(ote, .vars = c("trend", "trend.se", "trend.prob90", "trend.lower95", "trend.upper95", "trend.lower90", "trend.upper90"), as.numeric)


ote <- cbind(indicator = paste0("VESLA_COD_trend"),
             s.year = rep(min(y$na),nrow(ot)),
             e.year = rep(max(y$na),nrow(ot)),
             ote,
             trend.r = exp(ote$trend),
             trend.r.prob90 = exp(ote$trend.prob90),
             trend.r.lower95 = exp(ote$trend.lower95),
             trend.r.upper95 = exp(ote$trend.upper95),
             trend.r.lower90 = exp(ote$trend.lower90),
             trend.r.upper90 = exp(ote$trend.upper90),
             trend.interpretation = selected_status$status_class,
             trend.support = selected_status$support)


return(ote)


}

output_table_tidy_VESLA_GAM <- function(x) {
 
 resolution <- (max(x$data$na)-min(x$data$na))+1
  
 ot <- conditional_effects(x, resolution = resolution)
 otdf <- as.data.frame(ot$`na:Jarvi_tyyppi`)
  
 otF <- otdf[,c("na","Jarvi_tyyppi3","estimate__","lower__","upper__")]
 otE <- otdf[,c("na","Jarvi_tyyppi3","estimate__","lower__","upper__")]
 otS <- otdf[,c("na","Jarvi_tyyppi3","estimate__","lower__","upper__")]
 
 otE$Jarvi_tyyppi3 <- as.character(otE$Jarvi_tyyppi3)
 otE$Jarvi_tyyppi3[which(otE$Jarvi_tyyppi3=="humuksinen")] <- "humic"
 otE$Jarvi_tyyppi3[which(otE$Jarvi_tyyppi3=="vähähumuksinen")] <- "non-humic"
 otE$Jarvi_tyyppi3 <- as.factor(otE$Jarvi_tyyppi3)
 
 
 colnames(otF) <- c("vuosi", "järvityyppi", "kemiallinen hapenkulutus (mg/l)","95% luottamusvälin alaraja", "95% luottamusvälin yläraja")
 colnames(otE) <- c("year", "lake type", "chemical oxygen demand (mg/l)","lower boundary of the 95% credible interval", "upper boundary of the 95% credible interval")
 colnames(otS) <- c("År", "sjö typ", "kemisk syreförbrukning (mg/l)","nedre gränsen för det 95% trovärdiga intervallet", "övre gränsen för det 95 % trovärdiga intervallet")

 
  write.table(otF, file="Output/COD_GAM_fit_tidy_FIN.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=",", quote = FALSE)
  write.table(otE, file="Output/COD_GAM_fit_tidy_ENG.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=".", quote = FALSE)
  write.table(otS, file="Output/COD_GAM_fit_tidy_SWE.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=".", quote = FALSE)
  
  
  return(otE)

}

output_table_wide_VESLA_GAM <- function(x) {
  
  ot1 <- x[,c("na","Jarvi_tyyppi3","estimate__")]
  ot2 <- x[,c("na","Jarvi_tyyppi3","lower__")]
  ot3 <- x[,c("na","Jarvi_tyyppi3","upper__")]
  
  ote1 <- tidyr::pivot_wider(ot1, names_from = "Jarvi_tyyppi3",
                     values_from = "estimate__")
  ote1.1 <- ote1[,c("na","humuksinen")]
  colnames(ote1.1) <- c("vuosi","humusjärven kem. hapenkulutus")
  ote1.2 <- ote1[,c("vähähumuksinen")]
  colnames(ote1.2) <- c("vähähumuksisen järven kem. hapenkulutus")  
  
  ote2 <- tidyr::pivot_wider(ot2, names_from = "Jarvi_tyyppi3",
                             values_from = "lower__")
  ote2.1 <- ote2[,c("humuksinen")]
  colnames(ote2.1) <- c("humusjärven kem. hapenkulutuksen 95% luottamusvälin alaraja")
  ote2.2 <- ote2[,c("vähähumuksinen")]
  colnames(ote2.2) <- c("vähähumuksisen järven kem. hapenkulutuksen 95% luottamusvälin alaraja") 
  
  ote3 <- tidyr::pivot_wider(ot3, names_from = "Jarvi_tyyppi3",
                             values_from = "upper__")
  ote3.1 <- ote3[,c("humuksinen")]
  colnames(ote3.1) <- c("humusjärven kem. hapenkulutuksen 95% luottamusvälin yläraja")
  ote3.2 <- ote3[,c("vähähumuksinen")]
  colnames(ote3.2) <- c("vähähumuksisen järven kem. hapenkulutuksen 95% luottamusvälin yläraja")
  
  ot.dfF <- cbind(ote1.1, ote2.1, ote3.1, ote1.2, ote2.2, ote3.2)
  ot.dfS <- ot.dfF
  colnames(ot.dfS) <- c("år", "kemisk syreförbrukning i humusrika sjöar (mg/l)",
                        "nedre gränsen för det 95% trovärdiga intervallet för kemisk syreförbrukning i humusrika sjöar", 
                        "övre gränsen för det 95 % trovärdiga intervallet för kemisk syreförbrukning i humusrika sjöar",
                        "kemisk syreförbrukning i sjöar med låg humus koncentration (mg/l)",
                        "nedre gränsen för det 95% trovärdiga intervallet för kemisk syreförbrukning i sjöar med låg humus koncentration",
                        "övre gränsen för det 95 % trovärdiga intervallet för kemisk syreförbrukning i sjöar med låg humus koncentration")
  
  ot.dfE <- ot.dfF
  colnames(ot.dfE) <- c("year", "chemical oxygen demand in naturally humic lakes (mg/l)",
                        "lower boundary of 95% credible interval for chemical oxygen demand in naturally humic lakes",
                        "upper boundary of 95% credible interval for chemical oxygen demand in naturally humic lakes",
                        "chemical oxygen demand in naturally non-humic lakes (mg/l)",
                        "lower boundary of 95% credible interval for chemical oxygen demand in naturally non-humic lakes",
                        "upper boundary of 95% credible interval for chemical oxygen demand in naturally non-humic lakes")
  
  write.table(ot.dfF, file="Output/COD_GAM_fit_wide_FIN.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=",")
  write.table(ot.dfE, file="Output/COD_GAM_fit_wide_ENG.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=".")
  write.table(ot.dfS, file="Output/COD_GAM_fit_wide_SWE.csv", fileEncoding = "UTF-8", row.names=FALSE, sep=";", dec=",")
  
  return(ot1df)
}
