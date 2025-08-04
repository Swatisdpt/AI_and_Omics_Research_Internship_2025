#creating subfolders
dir.create("raw_data")
dir.create("clean_data")
dir.create("scripts")
dir.create("results")
dir.create("plots")

#loading csv file
patient_data <- read.csv(file.choose())

#inspecting the structure of the dataset
str(patient_data)

#converting variables to appropriate data types
patient_data$gender <- as.factor(patient_data$gender)
class(patient_data$gender)

patient_data$diagnosis <- as.factor(patient_data$diagnosis)
class(patient_data$diagnosis)

patient_data$smoker <- as.factor(patient_data$smoker)
class(patient_data$smoker)

#creating new variable for smoking status as binary factors 1 for yes and 0 for no
patient_data$smoking_status <- factor(patient_data$smoker,
                                      levels = c("Yes", "No"),
                                      labels = c(1, 0))
str(patient_data)

#save the cleaned dataset as csv file
write.csv(patient_data, file = "clean_data/patient_info_clean.csv")

#saving R script in script folder
