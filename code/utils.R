reload_source <- function(){
                library(tidyverse)
                library(lubridate)
                
                library(cowplot)
                
                library(RColorBrewer)
                library(viridis)
                library(flextable)
                library(ggpubr)
                #library(officer)
                
                library(randomForest)
                library(cvAUC)
                
                library(icenReg)
                library(survival)
                library(rmcorr)
                
                
                library(rpart)
                library(rpart.plot)
                
                
                
                source("code/utils.R")
}


#Figure 1

pControl<- function(data,title) {
                data %>%
                                gather(biomarker,value,-c(sample,dos)) %>%
                                ggplot(aes(x=dos,y=value#,col=biomarker
                                ))+
                                geom_jitter(alpha=0.4,
                                            height = 0,
                                            size=0.8,
                                            shape=16) +
                                ylab(title)+
                                facet_grid(biomarker~.,switch="y")+ theme_cowplot()+
                                theme(axis.title.x = element_blank(),
                                      axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank(),
                                      panel.grid.major.x = element_blank(),
                                      panel.grid.minor.x = element_blank()
                                )
                
}



pCase<- function(data) {
                data %>%
                                ggplot(aes(x=dos,y=value))+
                                geom_point(#aes(col=sample %in% olddata$sample),
                                           # alpha=0.4,
                                           # size=1,
                                           alpha=0.4,
                                           size=0.8,
                                           shape=16) + 
                                geom_smooth(se=FALSE,alpha=0.3,
                                            method="loess",
                                            fullrange = TRUE,
                                            xseq = seq(0,100, length=200)
                                            
                                            )+
                                geom_line(aes(group=id),alpha=0.05,size=0.3)+
                                xlab("Days since symptom onset")+
                                scale_x_continuous(breaks = seq(0,126,14)
                                                   )+
                                theme_cowplot()+
                                theme(
                                                panel.grid.minor.x = element_blank(),
                                                #axis.title.x = element_blank(),
                                                axis.title.y = element_blank(),
                                                axis.text.y = element_blank(),
                                                axis.ticks.y = element_blank(),
                                )
}



cv_titer_repeated_measure_auc <- function(dat,k,biomarker,col_of_interest){
                library(cvAUC)
                
                titers <- dat[,biomarker] %>% unlist
                truths <- dat[,col_of_interest] %>% unlist
                ids <- dat[,'id'] %>% unlist
                
                folds <- cv_folds_rm(dat,k)
                preds <- vector("list",length=k)
                
                
                
                for (k in 1:length(folds)){
                                perf <- ROCR::performance(ROCR::prediction(titers[-folds[[k]]],truths[-folds[[k]]]), "sens", "spec")
                                df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
                                cutoff <- df[which.max(df$sens + df$spec), "cut"]
                                
                                preds[[k]] <- ifelse(titers[folds[[k]]]>=cutoff,1,0) %>% unlist %>% unname
                }
                
                out <- ci.pooled.cvAUC(predictions=preds %>% unlist,
                                       labels=truths[unlist(folds)],
                                       ids=ids[unlist(folds)],
                                       folds=rep(1:k,map(folds,length)),
                                       confidence=0.95)
                
                return(out)
}


cv_folds_rm <- function(dat, k){
                
                ## get number of observations per person
                unique_ids <- unique(dat$id) ## unique ids in the dataset
                obs_indices <- lapply(1:length(unique_ids),function(x) which(dat$id == unique_ids[x])) ## list of observations per person
                ## get warnings for not having equal values in splits
                person_id_folds <- suppressWarnings(split(unique_ids,rep(1:k,length(unique_ids))))
                obs_id_folds <- lapply(person_id_folds,
                                       function(my_fold){which(dat$id %in% my_fold)})
                
                return(obs_id_folds)
                
}






## cross validated AUC for RF models
cv_rf_model <- function(dat,my_formula,k=20,ntrees=2000,...){
                
                ## creating folds
                folds <- cv_folds_rm(dat,k)
                
                my_preds <- my_truths <- numeric()
                #window <- str_extract(my_formula[2] %>% as.character,"inf_[0-9_]+")
                tmp_pred <- var_importance <- vector("list",length=k)
                
                
                list(folds)
                
                for (f in 1:length(folds)){
                                #cat(sprintf("Cross-validating AUC, fold %s \n",f))
                                my_forest <- randomForest(my_formula,
                                                          ntree=ntrees,
                                                          data=dat[-c(folds[[f]]),],
                                                          importance=TRUE,...)

                                my_truths <- c(my_truths,dat[folds[[f]],"case"] %>% unlist)
                                tmp_pred[[f]] <- predict(my_forest,newdata=dat[c(folds[[f]]),],predict.all=TRUE)
                                var_importance[[f]] <- importance(my_forest)

                                ## get probability from trees
                                my_preds <- c(my_preds,
                                              apply(tmp_pred[[f]][[2]],1,function(x) mean(x==1)) %>%
                                                              as.numeric
                                )
                                
                                
                }
                

                out <- vector("list",length=4)

                out[[1]] <- ci.pooled.cvAUC(predictions=my_preds,
                                            labels=my_truths,
                                            ids=dat[unlist(folds),'id'] %>% unlist,
                                            folds=rep(1:k,map(folds,length)),
                                            confidence=0.95) %>%
                                unlist %>% t %>%
                                data.frame

                out[[1]] <- out[[1]] #%>% mutate(time_window=window)

                out[[2]] <- bind_cols(truth=my_truths,pred=my_preds,fold=rep(1:k,map(folds,length)))

                out[[3]] <- tmp_pred

                out[[4]] <- var_importance

                names(out) <- c("perf_summary","outs","preds_full","var_imp")
                return(out)
}


## ROCs for cross validated results (multiple)
make_rocs <- function(preds,truths,k_folds,ribbon=TRUE,title="",annot_auc=TRUE){
                votes <- sapply(preds,function(pred)
                                apply(pred$individual,1,function(x) mean(x==1)),simplify=FALSE)
                
                perf_auc <- sapply(1:k_folds,function(x)
                                prediction(votes[[x]],truths[[x]]) %>% performance(.,"auc"))
                
                auc <- sapply(perf_auc,function(x) x@y.values[[1]])
                
                perf_roc <- sapply(1:k_folds,function(x) {
                                prediction(votes[[x]],truths[[x]]) %>% performance(.,"tpr","fpr")},
                                simplify = F)
                
                df <- data.frame(fold = rep(1:k_folds,sapply(perf_roc,function(x) x@x.values %>% unlist %>% length)),
                                 fpr=sapply(perf_roc,function(x) x@x.values[[1]],simplify=FALSE) %>% unlist,
                                 tpr=sapply(perf_roc,function(x) x@y.values[[1]],simplify=FALSE) %>% unlist)
                
                
                # cat(sprintf("AUC Range: %.2f-%.2f \n",range(auc)[1],range(auc)[2]))
                # cat(sprintf("AUC IQR: %.2f-%.2f \n",quantile(auc,.25),quantile(auc,.75)))
                # cat(sprintf("AUC Median and Mean: %.2f, %.2f \n",median(auc),mean(auc)))
                
                df <- df %>% dplyr::mutate(fpr_r = round(fpr,2),tpr_r=round(tpr,2))
                ci_ribbon <- df %>% group_by(fpr_r) %>% dplyr::summarize(med=median(tpr),low_ci=quantile(tpr,.025),upper_ci=quantile(tpr,.975))
                
                if(ribbon){
                                
                                ggplot(data=ci_ribbon) +
                                                geom_line(aes(x=fpr_r,y=med),col='royalblue',lwd=1.2)+
                                                geom_ribbon(data=ci_ribbon,aes(x=fpr_r,ymin=low_ci,ymax=upper_ci),alpha=.5,fill='royalblue') +
                                                xlab('False Positive Rate') + ylab('True Positive Rate') -> gg
                                
                } else {
                                
                                ggplot(data=df,aes(x=fpr,y=tpr)) +
                                                geom_line(aes(group=fold),color=AddAlpha('#2c7fb8',0.9)) +
                                                xlab('False Positive Rate') + ylab('True Positive Rate') +
                                                annotate("text", x = 0.5, y = 0.15, size=2,label = title) -> gg
                                
                                if(annot_auc) gg <- gg + annotate("text", x = 0.5, y = 0.1,size=2, label = paste0("AUC Range=",round(min(auc),3),"-",round(max(auc),3)))
                                
                }
                
                return(list(gg,df,auc))
                
}


make_roc_varimp_plot <- function(preds,truths,imps,my_title,panel_label="",sub=FALSE,ribbon=TRUE,...){
                
                my_roc <- make_rocs(preds,truths,length(preds),title=my_title,ribbon=ribbon,...)
                
                importance_df <- do.call('rbind',imps) %>% as.data.frame()

                importance_df$variable <- rownames(importance_df) %>% str_remove("\\.[0-9]") #%>% pretty_antibody
                
                levels_vars <- importance_df %>%
                                group_by(variable) %>%
                                dplyr::summarize(med=median(MeanDecreaseAccuracy)) %>%
                                select(med) %>% unlist() %>% order()

                impplot <- importance_df %>%
                                dplyr::mutate(variable1=factor(variable,
                                                               levels=sort(unique(importance_df$variable))[levels_vars])) %>%
                                ggplot() +
                                geom_bar(aes(x=variable1,y=MeanDecreaseAccuracy,fill=variable),stat='identity',alpha=0.5) +
                                coord_flip() + scale_fill_discrete(guide=F) + xlab('') + ylab('relative importance') +
                                theme_classic() + theme(axis.title=element_text(size=8),
                                                        axis.text.x = element_text(size=6))

                if(sub){
                                gg <- my_roc[[1]] + annotate("text",x=0,y=0,label=my_title,hjust=0) +
                                                annotate("text",x=0,y=1,label=panel_label,hjust=0) +
                                                annotation_custom(ggplotGrob(impplot), xmin=0.35, xmax=1, ymin=0, ymax=.4)

                                return(gg)

                } else{

                                return(multiplot(impplot,my_roc[[1]],cols=2))

                }
}

##' Adds alpha to a set of colors
##' @title
##' @param COLORS
##' @param ALPHA
##' @return
AddAlpha <- function(COLORS, ALPHA){
                if(missing(ALPHA)) stop("provide a value for alpha between 0 and 1")
                RGB <- col2rgb(COLORS, alpha=TRUE)
                RGB[4,] <- round(RGB[4,]*ALPHA)
                NEW.COLORS <- rgb(RGB[1,], RGB[2,], RGB[3,], RGB[4,], maxColorValue = 255)
                return(NEW.COLORS)
}

## some helper functions
sens_func <- function(x,truths){
                mean((truths==x)[which(truths==1)])
}

spec_func <- function(x,truths){
                mean((truths==x)[which(truths==0)])
}

youden_func <- function(x,truths,w_spec=.5){
                (1-w_spec)*sens_func(x,truths) + w_spec*spec_func(x,truths) - 1
}

##' Uses a cross validation approach to get the opimimal cutpoints
##' then tests performance on hold-out. Only uses one observation per
##' person
##' @param dat
##' @param my_threshold threshold for postivie/negative (if known), otehrwise will compute based on weighted youden index
##' @param w_spec #desired minimum specificity
##' @param biomarker
##' @param col_of_interest
##' @param nsims
##' @param print_me
##' @param training_prop
##' @return
##' @author
test_thresholds <- function(dat,
                            my_threshold=NA,
                            w_spec=1,
                            biomarker,
                            col_of_interest="case", #time period
                            nsims=1000,
                            print_me=FALSE,
                            training_prop=0.7){
                
                unique_titers <- dat[,biomarker] %>% unlist() %>% signif(.,3) %>% unique %>% sort
                
                # if(biomarker=='vib'){
                #                 ## get unique values of the titer
                #                 unique_titers <- c(dat$vibinab,dat$vibogaw) %>% unique %>% sort
                #
                #                 ## add new column to dat for max_vib
                #                 dat <- dat %>% dplyr::mutate(vib=pmax(vibogaw,vibinab))
                #
                # } else {
                #                 if(!biomarker %in% colnames(dat)) stop("biomarker is not a valid name of a column in dat. Only non-column name allowed is 'vib,' which is for max vibriocidal")
                #
                #                 unique_titers <- round(dat[,biomarker]) %>% unlist %>% unique %>% sort
                # }

                unique_ids <- unique(dat$id) ## unique ids in the dataset
                obs_indices <- sapply(1:length(unique_ids),
                                      function(x) which(dat$id == unique_ids[x])) ## list of observations per person
               
                ## number of people in training set
                n_train <- (length(unique_ids) * training_prop) %>% ceiling

                truths <- preds <- vector(mode='list', length=nsims)
                testing_perf <- matrix(ncol=5,nrow=nsims) %>% data.frame()
                colnames(testing_perf) <- c('threshold','sens_test','spec_test','sens_train','spec_train')

                for (i in 1:nsims){

                                cat("Simulation:", i, "of", nsims,"\n")

                                train_ids <- sample(length(unique_ids),n_train,replace=F) ## get ids of those who we want to sample for training
                                test_ids <- setdiff(1:length(unique_ids),train_ids) ## then those for testing

                                ## sample one time point for each person
                                train_obs_inds <- sapply(train_ids,function(x) sample(obs_indices[[x]],1))
                                test_obs_inds <- sapply(test_ids,function(x) sample(obs_indices[[x]],1))

                                training_titers <- dat[train_obs_inds,biomarker]
                                training_truths <- dat[train_obs_inds,col_of_interest]

                                testing_titers <- dat[test_obs_inds,biomarker]
                                testing_truths <- dat[test_obs_inds,col_of_interest]
                                
                                

                                # ## now get optimal cut-point (if needed) from the
                                ## training set
                                if(is.na(my_threshold)){
                                                ## returns a matrix with columns for each possible cutpoint and rows for observations
                                                preds_training <- sapply(unique_titers,function(x) training_titers>=x)

                                                sens_max_id <-  which.max(apply(preds_training,2,function(x) sens_func(x,truths=training_truths)*
                                                                                (spec_func(x,truths=training_truths)>=w_spec)))

                                                sens_spec_thresh <- c(unique_titers[sens_max_id],
                                                                      apply(preds_training,2,function(x) sens_func(x,truths=training_truths))[sens_max_id],
                                                                      apply(preds_training,2,function(x) spec_func(x,truths=training_truths))[sens_max_id])



                                                my_thresh <- sens_spec_thresh[1]
                                                sens_train <- sens_spec_thresh[2]
                                                spec_train <- sens_spec_thresh[3]

                                } else {

                                                my_thresh <- my_threshold
                                                preds_training <- training_titers>=my_thresh

                                                sens_train <- sens_func(preds_training,truths=training_truths)
                                                spec_train <- spec_func(preds_training,truths=training_truths)

                                }

                                ## now using the testing data, what is the sens and spec from these cut-offs?
                                testing_preds <- data.frame(pred=testing_titers >= my_thresh,
                                                            truth=testing_truths)

                                testing_perf[i,] <- c(threshold=my_thresh,
                                                      sens_test=sens_func(testing_preds[,1],testing_preds[,2]),
                                                      spec_test=spec_func(testing_preds[,1],testing_preds[,2]),
                                                      sens_train=sens_train,
                                                      spec_train=spec_train)

                }

                ## don't want to keep this loaded
                library(reshape2)
                thresh_df <- melt(testing_perf,id.vars='threshold')
                detach('package:reshape2')

                if (print_me){
                                thresh_df %>% ggplot() +
                                                geom_histogram(aes(value,fill=factor(threshold)),position='identity') +
                                                facet_wrap(~variable) +
                                                labs(title=col_of_interest)  -> gg
                                print(gg)
                }

                return(thresh_df)
}


##' Uses a cross validation approach to get the opimimal cutpoints
##' then tests performance on hold-out. Only uses one observation per
##' person
##' @param dat
##' @param my_threshold threshold for postivie/negative (if known), otehrwise will compute based on weighted youden index
##' @param w_spec #desired minimum specificity
##' @param biomarker
##' @param col_of_interest
##' @param nsims
##' @param print_me
##' @param training_prop
##' @return
##' @author
test_thresholdRF <- function(dat,
                            my_threshold=NA,
                            w_spec=1,
                            col_of_interest="case", 
                            nsims=1000,
                            print_me=FALSE,
                            training_prop=0.7){
                
                ntree=1000
                unique_titers <- (1:ntree)/ntree # tree values 
                #unique_titers <- (1:round(ntree*0.7))/ntree # tree values 
                
                
                unique_ids <- unique(dat$id) ## unique ids in the dataset
                obs_indices <- sapply(1:length(unique_ids),
                                      function(x) which(dat$id == unique_ids[x])) ## list of observations per person
                
                ## number of people in training set
                n_train <- (length(unique_ids) * training_prop) %>% ceiling
                
                truths <- preds <- vector(mode='list', length=nsims)
                testing_perf <- matrix(ncol=5,nrow=nsims) %>% data.frame()
                colnames(testing_perf) <- c('threshold','sens_test','spec_test','sens_train','spec_train')
                
                for (i in 1:nsims){

                                cat("Simulation:", i, "of", nsims,"\n")

                                train_ids <- sample(length(unique_ids),n_train,replace=F) ## get ids of those who we want to sample for training
                                test_ids <- setdiff(1:length(unique_ids),train_ids) ## then those for testing
                                                
                                ## sample one time point for each person
                                train_obs_inds <- sapply(train_ids,function(x) sample(obs_indices[[x]],1))
                                test_obs_inds <- sapply(test_ids,function(x) sample(obs_indices[[x]],1))

                                training_set <- dat[train_obs_inds,c(col_of_interest,"IgA","IgG","IgM")]
                                trainingRF<-training_set %>%
                                                randomForest(factor(case)~IgM+IgA+IgG,
                                                             ntree=ntree,
                                                             data=.)
                                training_truths <- dat[train_obs_inds,col_of_interest]


                                testing_set <- dat[test_obs_inds,c(col_of_interest,"IgA","IgG","IgM")]
                                testing_truths <- dat[test_obs_inds,col_of_interest]

                                # ## now get optimal cut-point (if needed) from the
                                ## training set
                                if(is.na(my_threshold)){
                                                ##

                                                ## returns a matrix with columns for each possible cutpoint and rows for observations
                                                preds_training <- sapply(unique_titers,function(x) trainingRF$votes[,2]>=x)

                                                sens_max_id <-  which.max(apply(preds_training,2,function(x)
                                                                sens_func(x,truths=training_truths)*
                                                                                (spec_func(x,truths=training_truths)>=w_spec)))

                                                sens_spec_thresh <- c(unique_titers[sens_max_id],
                                                                      apply(preds_training,2,function(x) sens_func(x,truths=training_truths))[sens_max_id],
                                                                      apply(preds_training,2,function(x) spec_func(x,truths=training_truths))[sens_max_id])

                                                my_thresh <- sens_spec_thresh[1]
                                                sens_train <- sens_spec_thresh[2]
                                                spec_train <- sens_spec_thresh[3]
                                                

                                } else {

                                                my_thresh <- my_threshold
                                                preds_training <- trainingRF$votes[,2]>=my_thresh

                                                sens_train <- sens_func(preds_training,truths=training_truths)
                                                spec_train <- spec_func(preds_training,truths=training_truths)


                                                
                                                }

                                ## now using the testing data, what is the sens and spec from these cut-offs?
                                testing_preds <- data.frame(pred=predict(trainingRF, testing_set, type="prob")[,2]>= my_thresh,
                                                            truth=testing_truths)
                                
                                
                                

                                testing_perf[i,] <- c(threshold=my_thresh,
                                                      sens_test=sens_func(testing_preds[,1],testing_preds[,2]),
                                                      spec_test=spec_func(testing_preds[,1],testing_preds[,2]),
                                                      sens_train=sens_train,
                                                      spec_train=spec_train)
                                
                                

                }

                ## don't want to keep this loaded
                library(reshape2)
                thresh_df <- melt(testing_perf,id.vars='threshold')
                detach('package:reshape2')

                if (print_me){
                                thresh_df %>% ggplot() +
                                                geom_histogram(aes(value,fill=factor(threshold)),position='identity') +
                                                facet_wrap(~variable) +
                                                labs(title=col_of_interest)  -> gg
                                print(gg)
                }

                return(thresh_df)
}


##' takes the cases data frame to generate TL and TR
##' for both seroCONversion and seroREVERSION.
##' assumes that the first serconversion is the true one
##' assumes the first seroreversion is false, observation is removed
##' @param data input the cases you are interested in
##' @param iso isotype desired
##' @return a list 1.seroconversion df, 2. seroreversion df, 3. double seroconverters
convertTLTR <- function(data, iso){
                
                #remove double seroconverting measurements
                if(iso=="IgG"){
                                data <- filter(data, !(id==57 & dos>23))
                               
                }
                if(iso=="IgM"){
                                data <- filter(data, !(id==652 & dos>38))
                                
                }
                if(iso=="IgA"){
                                data <- filter(data, !(id==70 & dos>10))
                                data <- filter(data, !(id==180 & dos>51))
                                data <- filter(data, !(id==610 & dos>51))
                                
                }
                
                
                #seropositivity cutoffs
                sero_cut <- data.frame(IgG=0.57,#0.5,  #IgG
                                       IgM=2.63,#0.68, #IgM
                                       IgA=2.02)#2.02) #IgA
                cutp <- sero_cut[iso] %>% unlist()
                
                #get dataset together
                # data <- filter(data,!is.na(Severity2)) #why is this here???
                init1 <- data[c("id","sample",iso,"dos")]
                names(init1)[3] <- "Ig"
                
                #define seroconversion and seroreversion
                init2 <- init1 %>% 
                                filter(!(is.na(dos) | is.na(Ig)))%>%
                                arrange(id,dos) %>%
                                # serostatus
                                mutate(seroPos=ifelse(Ig>=cutp,1,0)) %>%
                                group_by(id) %>%
                                #cumulative serostatus
                                mutate(cumSeroPos=cumsum(seroPos)) %>%
                                # 0 = always seronegative
                                # 1 = seropositive
                                # 2 = seroreverted
                                mutate(ConRev=ifelse(seroPos==0 & cumSeroPos>0,2,seroPos))%>% 
                                mutate(ConRev2 = lead(ConRev,1)-ConRev)                
                
                # check for double sero-converters 
                neg <- init2  %>% filter(ConRev2<0)
                double <- filter(data,id %in% neg$id) 
                double <- double[c("id","sample",iso,"dos")]
                names(double)[3] <- "Ig"
                double <- double %>% mutate(seroPos=ifelse(Ig>=cutp,1,0))
                
                if(nrow(double)>0) print("CHECK FOR DOUBLE SEROCONVERTERS")
                
                #i think we just remove the double seroconverters and do those ones manually
                
                init3 <- init2 %>% 
                                #remove the double seroconverters 
                                filter(!(id %in% neg$id))%>%
                                # # remove first seroreversion measurement
                                # # for double converters
                                # filter(!(sample %in% neg$sample))%>%
                                # for each sero-state for each individual, find max and min
                                group_by(id,ConRev) %>%
                                summarize(min=min(dos),max=max(dos))
                
                
                
                # dataframe for seroconversion
                convDf <- init3 %>%
                                # remove sero-reverters
                                filter(ConRev!=2) %>%
                                # find TL and TR
                                mutate(day=ifelse(ConRev==1,min,max)) %>%
                                mutate(interval = recode(ConRev,
                                                         `0`="TL",
                                                         `1`="TR")
                                ) %>%
                                select(id,interval,day) %>%
                                spread(interval,day) %>%
                                # add in zeroes and infinity for missing
                                mutate(TL=ifelse(is.na(TL),0,TL),
                                       TR=ifelse(is.na(TR),Inf,TR))

                # summary(convDf$TR-convDf$TL)
                # #this individual was seropositive on the first day
                # hi <- cases %>% filter(id==225)


                # dataframe for seroreversion
                # for those who always negative
                revDf1 <- init3 %>%
                                group_by(id) %>%
                                filter(max(ConRev)==0) %>%
                                # find TL (or other assumption)
                                mutate(TL=max, TR=Inf) %>%
                                select(id,TL,TR)

                #remove initial seronegativity
                revDf2 <- init3 %>%
                                filter(ConRev!=0) %>%
                                # find TL and TR
                                mutate(day=ifelse(ConRev==2,min,max)) %>%
                                mutate(interval = recode(ConRev,
                                                         `1`="TL",
                                                         `2`="TR")
                                ) %>%
                                select(id,interval,day) %>%
                                spread(interval,day) %>%
                                mutate(TR=ifelse(is.na(TR),Inf,TR))

                #combine
                revDf <- bind_rows(revDf1,revDf2)
                # summary(revDf$TR-revDf$TL)

                #bring individual characterisistics
                indy <- data %>% select(id,cohort,gender,age,ageCat,req_ICU,death,immunosuppressed,
                                        Severity,Severity2,Severity3

                                        )%>%
                                distinct(id,cohort,gender,age,ageCat,req_ICU,death,immunosuppressed,
                                         Severity,Severity2,Severity3)


                convDf <-left_join(convDf,indy)
                revDf <-left_join(revDf,indy)

                return(list(convDf,revDf,double))
                # return(init3)
                
                
}

##'  find bootstrap 95% CI of sensitivity for a given cutoff
##' @param dat input the prediction data
##' @param biomarker isotype desired
##' @param threshold threshold of importance
##' @param col_of_interest case vs control
bootstrapSensSpec<- function(dat,
                             nsims=1000,
                             biomarker,
                             threshold,
                             col_of_interest){
                
                sens<-NULL; spec<-NULL
                
                unique_ids <- unique(dat$id) ## unique ids in the dataset
                obs_indices <- sapply(1:length(unique_ids),
                                      ## list of observations per person
                                      function(x) which(dat$id == unique_ids[x])) 
                
                for (i in 1:nsims){
                                if((i %% (nsims/10))==0) cat("Simulation: ",i,"of",nsims,"\n")
                                
                                getIDs<- sample(length(unique_ids),length(unique_ids),replace = TRUE)
                                obs_inds <- sapply(getIDs,function(x) sample(obs_indices[[x]],1))
                                
                                
                                cases<- dat[obs_inds,col_of_interest]
                                titers<- dat[obs_inds,biomarker]
                                
                                sens <- c(sens,mean(titers[cases==1]>=threshold))
                                spec <- c(spec,mean(titers[cases==0]< threshold))
                                }
                return(data.frame(biomarker,threshold,sens,spec))
}


##'  find bootstrap 95% CI of sensitivity for a cmobination of cutoffs
##' @param dat input the prediction data
##' @param biomarker isotype desired
##' @param threshold threshold of importance
##' @param col_of_interest case vs control

double.bootstrapSensSpec<- function(dat,
                                    nsims=1000,
                                    biomarker1,threshold1,
                                    biomarker2,threshold2,
                                    biomarker3=NA,threshold3=NA,
                                    col_of_interest){
                
                sens<-NULL; spec<-NULL
                
                unique_ids <- unique(dat$id) ## unique ids in the dataset
                obs_indices <- sapply(1:length(unique_ids),
                                      function(x) which(dat$id == unique_ids[x])) ## list of observations per person
                for (i in 1:nsims){
                                
                                if((i %% (nsims/10))==0) cat("Simulation: ",i,"of",nsims,"\n")
                                
                                getIDs<- sample(length(unique_ids),length(unique_ids),replace = TRUE)
                                obs_inds <- sapply(getIDs,function(x) sample(obs_indices[[x]],1))
                                
                                cases<- dat[obs_inds,col_of_interest]
                                titers1<- dat[obs_inds,biomarker1]
                                titers2<- dat[obs_inds,biomarker2]
                                
                                #two thresholds
                                if(is.na(biomarker3)){
                                                TP<-titers1[cases==1]>=threshold1 | titers2[cases==1]>=threshold2
                                                TN<-titers1[cases==0]< threshold1 & titers2[cases==0]< threshold2
                                }
                                
                                #three thresholds
                                if(!is.na(biomarker3)){
                                                
                                                titers3<- dat[obs_inds,biomarker3]
                                                TP<-titers1[cases==1]>=threshold1 | 
                                                                titers2[cases==1]>=threshold2 | 
                                                                titers3[cases==1]>=threshold3
                                                TN<-titers1[cases==0]<threshold1  & 
                                                                titers2[cases==0]<threshold2  & 
                                                                titers3[cases==0]<threshold3  
                                }
                                sens <- c(sens,mean(TP))
                                spec <- c(spec,mean(TN))
                }
                
                return(data.frame(biomarker1,threshold1,biomarker2,threshold2,
                                  biomarker3,threshold3,sens,spec))
                
}




##'  calculate bootstrap for parameters of accelerrated failure time
##' @param dat input the prediction data
##' @param isotype isotype desired
##' @param type case vs control
convertBootstrap <- function(data,
                             form=cbind(TL,TR) ~ 1,
                             distribution="weibull",
                             nsim=1000,
                             isotype, # isotype: IgM, IgA, IgG
                             type){ #convert or revert
                
                
                
                # if(isotype=="IgG") distribution="weibull"
                # if(isotype=="IgM") distribution="lnorm"
                # if(isotype=="IgA") distribution="gamma"
                
                
                #get data into interval
                dat <- convertTLTR(data,isotype)
                
                if(type=="convert") dat <- dat[[1]]
                if(type=="revert") dat <- dat[[2]]
                
                #initialize data.frame to store matrix
                storeBoot <- data.frame()
                for (b in 1:nsim){
                                #get bootstrap sample
                                s <- sample(1:nrow(dat),nrow(dat),replace = TRUE)
                                
                                #fit model
                                fit <- ic_par(form, model = "aft",
                                              dist = distribution,
                                              data = dat[s,])
                                
                                
                                
                                xtra <- data.frame(t(fit$coefficients),
                                                   sample=b)
                                
                                if(is.null(fit$xlevels)==FALSE){
                                                xtra <- data.frame(t(fit$coefficients),
                                                                   sample=b,
                                                                   ref=fit$xlevels[[1]][1])
                                }
                                
                                
                                #store parameters
                                storeBoot <- bind_rows(storeBoot,xtra)
                                
                                
                                
                                #     #use surv CIs
                                #     survCIs(fit_wICU,p=c(.025,.25,.5,.75,.975),newdata=newDataRI)$cis[[1]] %>%
                                # as.data.frame() %>% mutate(Pop="No ICU"),
                                # survCIs(fit_wICU,p=c(.025,.25,.5,.75,.975),newdata=newDataRI)$cis[[2]] %>%
                                # as.data.frame() %>% mutate(Pop="ICU"),
                                
                                if(b %% (nsim/10) ==0) cat(b,"bootstrap samples out of",nsim, "completed\n")
                                
                }
                
                storeBoot %>% mutate(isotype=isotype,type=type)
                
                
                
                
                
                
}

#assess different models by log likelihood
assessFit<- function(data){
                fitW <- data %>%
                                ic_par(cbind(TL,TR) ~1, model = "aft",dist = "weibull",data =.) 
                fitG <- data %>%
                                ic_par(cbind(TL,TR) ~1, model = "aft",dist = "gamma",data =.) 
                fitLN <- data %>%
                                ic_par(cbind(TL,TR) ~1, model = "aft",dist = "lnorm",data =.) 
                
                cat("weibull",fitW$llk,"\n")
                cat("gamma",fitG$llk,"\n")
                cat("lognormal",fitLN$llk,"\n")
                
                cat("MAX LL:",max(fitW$llk,fitLN$llk,fitG$llk))
                
}



# fits of the mle for both models of interval censored data
# (not the bootstraps)
getFits<- function(data,iso){
                
                #non-parametric
                np_fit <- ic_np(cbind(TL,TR) ~0, data)
                
                np_fit_df<-data.frame(surv=cumsum(np_fit$p_hat),
                                      t1=np_fit$T_bull_Intervals[1,],
                                      t2=np_fit$T_bull_Intervals[2,]
                ) %>%
                                mutate(surv2=lag(surv,1)) %>%
                                mutate(surv2=ifelse(is.na(surv2),0,surv2)) %>%
                                gather(type,surv,-c(t1,t2)) %>%
                                mutate(isotype=iso)%>%
                                arrange(-surv,-t1) #%>% # THIS IS EVERYTHING
                
                
                
                #parametric
                p_fit<- ic_par(cbind(TL,TR) ~1, 
                               model = "aft",dist = "lnorm",data=data)
                
                p_fit_df <- data.frame(
                                x=days,
                                survp=1-plnorm(days,
                                               (p_fit$coefficients[1]),
                                               exp(p_fit$coefficients[2]))) %>%
                                mutate(isotype=iso)
                
                
                return(list(np_fit_df,p_fit_df))
                
                
}


#easy graphs to compare subgroups
compareTopics <- function(data,topic){
                
                
                subG<- data[c("sample","IgG", "dos","id",topic)] 
                colnames(subG)[5] <-"topic"
                pG <- subG %>%        select(sample,IgG, dos,id,topic) %>%
                                gather(biomarker,value,-c(sample,dos,id,topic))%>%
                                ggplot(aes(x=dos,y=value,col=factor(topic)))+
                                geom_point(size=0.8,alpha=0.5) + 
                                geom_smooth(se=FALSE,alpha=0.3,
                                            method="loess",
                                            fullrange = TRUE,
                                            xseq = seq(0,100, length=200))+
                                geom_line(aes(group=id),alpha=0.05,size=0.3)+
                                xlab("Weeks since symptom onset")+
                                scale_x_continuous(breaks = seq(0,119,7),
                                                   labels = 0:17
                                )+
                                scale_y_continuous(breaks = c(0.1,10,1000),
                                                   labels = c(0.1,10,1000),
                                                   trans="log10",limits = c(0.01,1200))+
                                scale_color_brewer(topic,palette = "Dark2")+
                                ylab("IgG")+
                                geom_hline(yintercept = 0.57,lty=2,size=0.4)+
                                theme_cowplot()
                
                subM<- data[c("sample","IgM", "dos","id",topic)] 
                colnames(subM)[5] <-"topic"
                pM <- subM %>%        select(sample,IgM, dos,id,topic) %>%
                                gather(biomarker,value,-c(sample,dos,id,topic))%>%
                                ggplot(aes(x=dos,y=value,col=factor(topic)))+
                                geom_point(size=0.8,alpha=0.5) + 
                                geom_smooth(se=FALSE,alpha=0.3,
                                            method="loess",
                                            fullrange = TRUE,
                                            xseq = seq(0,100, length=200))+
                                geom_line(aes(group=id),alpha=0.05,size=0.3)+
                                xlab("Weeks since symptom onset")+
                                scale_x_continuous(breaks = seq(0,119,7),
                                                   labels = 0:17
                                )+
                                scale_y_continuous(breaks = c(0.1,10,1000),
                                                   labels = c(0.1,10,1000),
                                                   trans="log10",limits = c(0.01,1200))+
                                scale_color_brewer(topic,palette = "Dark2")+
                                ylab("IgM")+
                                geom_hline(yintercept = 2.63,lty=2,size=0.4)+
                                theme_cowplot()
                
                
                subA<- data[c("sample","IgA", "dos","id",topic)] 
                colnames(subA)[5] <-"topic"
                pA <- subA %>%        select(sample,IgA, dos,id,topic) %>%
                                gather(biomarker,value,-c(sample,dos,id,topic))%>%
                                ggplot(aes(x=dos,y=value,col=factor(topic)))+
                                geom_point(size=0.8,alpha=0.5) + 
                                geom_smooth(se=FALSE,alpha=0.3,
                                            method="loess",
                                            fullrange = TRUE,
                                            xseq = seq(0,100, length=200))+
                                geom_line(aes(group=id),alpha=0.05,size=0.3)+
                                xlab("Weeks since symptom onset")+
                                scale_x_continuous(breaks = seq(0,119,7),
                                                   labels = 0:17
                                )+
                                scale_y_continuous(breaks = c(0.1,10,1000),
                                                   labels = c(0.1,10,1000),
                                                   trans="log10",limits = c(0.01,1200))+
                                scale_color_brewer(topic,palette = "Dark2")+
                                ylab("IgA")+
                                geom_hline(yintercept = 2.02,lty=2,size=0.4)+
                                theme_cowplot()
                
                
                
                
                plot_grid(pG,pA,pM,ncol=1)
                
                
                
}


#use for Table S2
getCIcovariate <- function(boots,p){
                
                boots %>%
                                gather(variable,paramvalue,-c(isotype,type,sample,mu,log_s,ref)) %>%
                                
                                
                                #calculate using AFT
                                mutate(baseline=qlnorm(p,mu,exp(log_s)), #qlnorm(0.05,(mu),exp(log_s)
                                       covariate=baseline*exp(paramvalue),
                                       diff=covariate-baseline
                                ) %>%
                                group_by(isotype,type,variable,ref) %>%
                                summarize(
                                                #estimate for reference category
                                                Baseline=mean(baseline),
                                                lowerB=quantile(baseline,0.025),
                                                upperB=quantile(baseline,0.975),
                                                #estimate for non-reference category
                                                Covariate=mean(covariate),
                                                lowerCov=quantile(covariate,0.025),
                                                upperCov=quantile(covariate,0.975),
                                                #estunate difference
                                                Diff=mean(diff),
                                                lowerD=quantile(diff,0.025),
                                                upperD=quantile(diff,0.975)
                                ) %>%
                                mutate(
                                                #baseline
                                                Baseline=paste0(sprintf("%.1f",round(Baseline,1)),
                                                                " (",sprintf("%.1f",round(lowerB,1)),
                                                                "-",sprintf("%.1f",round(upperB,1)),")"),
                                                #Covariate
                                                Covariate=paste0(sprintf("%.1f",round(Covariate,1)),
                                                                 " (",sprintf("%.1f",round(lowerCov,1)),
                                                                 "-",sprintf("%.1f",round(upperCov,1)),")"),
                                                #Difference
                                                Difference=paste0(sprintf("%.1f",round(Diff,1)),
                                                                  " (",sprintf("%.1f",round(lowerD,1)),
                                                                  "-",sprintf("%.1f",round(upperD,1)),")"),
                                                
                                                
                                ) %>%
                                mutate(percentile=p)%>%
                                select(percentile,isotype,type,ref,variable,Baseline,Covariate,Difference)
                
}


#make confusion matrix plot
confusionmatrix<- function(rf) {
                
                title <-rf$call$formula %>% as.character()
                title <-title[3]
                
                cmat <- data.frame(rf$confusion, trueclass=c("Control","Case")) %>%
                                rename(Control=X0,Case=X1) %>%
                                select(-class.error)%>%
                                gather(predictclass,number,-trueclass) %>%
                                mutate(trueclass=factor(trueclass,levels=c("Control","Case")))%>%
                                mutate(percent=number/sum(number)) %>%
                                mutate(percent=round(100*percent,0))%>%
                                mutate(text=paste0(number,"\n",percent,"%"))
                
                
                erate <- data.frame(OOB=rf$err.rate[rf$ntree,1] ) %>%
                                mutate(OOB=paste0("OOB error rate: ",round(OOB*100,1),"%"))
                
                p1 <-cmat %>% ggplot(aes(x=predictclass,y=trueclass))+
                                geom_tile(aes(fill=percent))+
                                geom_text(aes(label=text))+
                                ggtitle(title)+
                                xlab("Predicted Class")+
                                ylab("True Class")+
                                theme_bw()+
                                scale_fill_distiller("Percent",palette = "Blues",direction = 1,
                                                     limits=c(0,100))+
                                theme(legend.position = "none",
                                plot.title = element_text(size = 10)
                                )
                
                p2 <-erate %>% ggplot()+
                                geom_text(x=0.5,y=0.5,aes(label=OOB))+
                                theme_void()
                
                
                plot_grid(p1,p2,ncol=1,rel_heights = c(5,1))
                
}











##' Uses AU
##'
##' .. content for \details{} ..
##' @title
##' @param dat
##' @param use_all
##' @param biomarker
##' @param col_of_interest
##' @return
##' 
# titer_auc <- function(dat,use_all=FALSE,biomarker,col_of_interest){
#                 
#                 if(use_all){
#                                 
#                                 titers <- dat[,biomarker]
#                                 truths <- dat[,col_of_interest]
#                                 
#                 } else {
#                                 
#                                 unique_ids <- unique(dat$id) ## unique ids in the dataset
#                                 obs_indices <- sapply(1:length(unique_ids),function(x) which(dat$id == unique_ids[x])) ## list of observations per person
#                                 
#                                 selected_obs_inds <- sapply(1:length(unique_ids),function(x) sample(obs_indices[[x]],1))
#                                 
#                                 titers <- dat[selected_obs_inds,biomarker]
#                                 truths <- dat[selected_obs_inds,col_of_interest]
#                                 
#                 }
#                 
#                 rc <- data.frame(titers,truths)
#                 colnames(rc) <- c("titers","truths")
#                 
#                 return(rc)
#                 
# }
# 
# titer_roc_plot <- function(tauc,annot_titers=TRUE){
#                 
#                 pred <- prediction(tauc$titer,tauc$truth)
#                 perf <- performance(pred,"tpr","fpr")
#                 auc <- performance(pred,"auc")
#                 x <- data.frame(x=perf@x.values[[1]],y=perf@y.values[[1]],titers=perf@alpha.values[[1]])
#                 
#                 gg <- ggplot(data=x) + geom_step(aes(x=x,y=y),direction = "hv") +
#                                 geom_point(aes(x=x,y=y),cex=0,alpha=.5) +
#                                 xlab("False Positive Rate") + ylab("True Positive Rate") +
#                                 annotate("text",x=.9,y=.1,label=sprintf("AUC=%.1f",100*auc@y.values[[1]]))
# 
#                 if(annot_titers)
#                                 gg <- gg +
#                                 geom_text(aes(x=x,y=y,label=titers),alpha=.7,nudge_x=.04)
# 
#                 print(gg)
#                 #print(auc)
#                 return(list(auc,x,gg))
# }
# 


