source("code/utils.R")
reload_source()


serodata <- read_csv("data/raw_data/Complied data 9_2_20 v2.csv",
                     col_types = cols(
                             sample_id=col_character(),
                             NT50 = col_character()
                     )) %>%
        rename(id=sample_id,
               timepoint=tp,
               IgA=iga_avg,
               IgM=igm_avg,
               IgG=igg_avg,
               discharge=dispo_at_testing,
               immunosuppressed=immunosuppressed_y_n,
               death = death_covid_related_y_n
        ) %>%
        
        mutate(
                req_ICU = factor(ifelse(is.na(req_ICU),0,1)), 
                gender=ifelse(gender %in% c(1,2),gender,NA),
                gender=factor(gender,labels=c('Male','Female')),
                ageCat=factor(cut(age,c(0,65,100),right = FALSE),
                              labels=c("<65 years","\u2265 65 years")), 
                death = factor(ifelse(is.na(death),0,1)),
                discharge = factor(discharge,labels=c("Discharged", "Hospitalized","ICU" )),
                immunosuppressed = factor(ifelse(is.na(immunosuppressed),0,1)),
                cohort=factor(cohort,levels=c(3,2,1),
                              labels=c("prepandemic","pcrneg","case")),
                
                #week
                #4 categories
                week=factor(cut(dos,c(-1,7,14,28,80),
                                right = TRUE),
                            labels= c("\u22647 days", "8-14 days",
                                      "15-28 days",">28 days")),
                #20 categories
                week2=factor(cut(dos,c(-1,seq(7,140,7)),
                                 labels=1:20)),
                week2=factor(week2,levels=c(1:20)),
                
                #Severity
                dos=ifelse(Severity==6,NA,dos), #asymptomatic folks
                
                
                Severity2=ifelse(Severity==5,NA,Severity),
                Severity2=ifelse(Severity2==6,1,Severity2),
                Severity2=factor(Severity2,
                                 labels=c("Not Hospitalized",
                                          "Hospitalized, no ICU",
                                          "Hospitalized, required ICU",
                                          "Died due to COVID-19")),
                Severity=ifelse(is.na(Severity),5,Severity),
                Severity=ifelse(Severity==6,1,Severity),
                Severity=factor(Severity,
                                labels=c("Not Hospitalized",
                                         "Hospitalized, no ICU",
                                         "Hospitalized, required ICU",
                                         "Died due to COVID-19",
                                         "Missing")),
                
                Severity3=recode(Severity2,
                                 "Not Hospitalized"= "Not Hospitalized / Hospitalized, no ICU",
                                 "Hospitalized, no ICU"="Not Hospitalized / Hospitalized, no ICU"
                )
                
                
                
        ) %>%
        #remove 4 prepandemic controls missing age
        filter(!is.na(age))


#prediction data
predData <- serodata %>% filter(cohort %in% c("prepandemic","case")) %>%
        filter(!(is.na(IgG)|is.na(IgA)|is.na(IgM))) %>%
        mutate(case=ifelse(cohort=="case",1,0))

#case and control data                
cases <- filter(serodata,cohort=="case")
controls <- filter(serodata,cohort=="prepandemic")
negatives <- filter(serodata,cohort=="pcrneg")

#individual characteristics
individual <- serodata %>% select(id,cohort,gender,age,ageCat,Race,
                                  req_ICU,death,immunosuppressed,Severity)%>%
        distinct(id,cohort,gender,age,ageCat,Race,
                 req_ICU,death,immunosuppressed,Severity)
