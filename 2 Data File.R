
library(dplyr)

####################
## caravan data ##
####################

data <- read.csv("cara_data.csv")
  head(data)
  str(data)
data <- data %>% rename(class = CARAVAN) 

remove_Start <- c("MGODGE",
                  "MRELOV",
                  "MFWEKIND",
                  "MOPLLAAG",
                  "MBERHOOG",
                  "MSKA",
                  "MHKOOP",
                  "MAUT0",
                  "MZPART",
                  "MINK123M",
                  "MINKGEM")
                  

train_data <- data %>% filter(ORIGIN == "train") %>% dplyr::select(-ORIGIN,-remove_Start)  %>% dplyr::select(kvs,class)
test_data <- data %>% filter(ORIGIN == "test") %>% dplyr::select(-ORIGIN,-remove_Start)  %>% dplyr::select(kvs,class)

num_vars <- ncol(train_data) - 1 # for origin
var_names <- colnames(train_data)[1:num_vars] 

# train_data[, num_vars + 1] <- as.numeric(train_data[, num_vars + 1])
# train_data[, num_vars + 2] <- as.numeric(train_data[, num_vars + 1])
# train_data[, num_vars + 1] <- as.numeric(train_data[, num_vars + 1])
# train_data[, num_vars + 2] <- as.numeric(train_data[, num_vars + 1])

rm(data)

###################################

### keep
library(vip)
var_remove <-   summary(model_use)
var_remove <- as.data.frame(cbind(var_remove$var, var_remove$rel.inf))
var_remove_list <- var_remove %>% filter(V2 >1.5) 
kvs <- var_remove_list$V1
str(kvs)

sum_of_distance
var_remove <- as.data.frame(cbind(names(sum_of_distance), as.numeric(sum_of_distance)))
var_remove_list <- var_remove %>% filter(sum_of_distance > 4550) 
kvs <- var_remove_list$V1
str(kvs)


### remove

rml <- c("AWAPART","APLEZIER","PPLEZIER")

train_data <- train_data %>% dplyr::select(-rml)
test_data <- test_data %>% dplyr::select(-rml)

biplot_data <- biplot_data %>% dplyr::select(-rml)
num_vars <- ncol(train_data) - 1 # for origin
var_names <- colnames(train_data)[1:num_vars] 


##################################
vars <- c("MOPLHOOG",
          "MOPLLAAG",
          "MINKM30",
          "MINKGEM",
          "MKOOPKLA",
          "PWAPART",
          "PPERSAUT",
          "PBRAND",
          "APERSAUT",
          "ALEVEN")

data 
apply(data %>% select(MOPLHOOG,MOPLMIDD,MOPLLAAG,), 1,"sum")


#[1] "PPERSAUT" "MOSTYPE"  "PBRAND"   "MOPLHOOG" "APERSAUT" "MKOOPKLA" "MAUT1"    "MINKM30"  "PWAPART" 