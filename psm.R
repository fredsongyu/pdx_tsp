
data <- read.delim("./psm_102018.txt")

# priliminary logit model
logit <- glm(tsp~ln_aadt+Oneway+num_lane+median_per
             +bus_routes+maj_freq+sec_freq+sig_den, 
             offset=log(len), family=binomial(link="logit"), data=data)
summary(logit)

# prediction of propensity scores
pred <- predict(logit,type="response")
data$propensity <- pred

write.table(data, file="matching.csv",sep=",",row.names=F)
