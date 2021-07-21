load("daphnia_DE.RData")

daphnia.tf.real <- read.csv(file = "TF_List_FULL.csv", header = TRUE)
daphnia.ma <- daphnia_fc$counts
colnames(daphnia.ma) <- c("A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4", "C1", "C2", "D1", "D2", "D3", "D4", "D5", "E1", "E2", "E3", "F1") #please double-check this- I generated these by hand.

daphnia.tf <- daphnia.tf.real[30:50,] #replace this line with your actual tf object

daphnia.df <- as.data.frame(daphnia.ma) #converting matrices to data.frames
daphnia.tf.df <- as.data.frame(daphnia.tf)

tf.merge <- merge(daphnia.df, daphnia.tf.df, by=0) #you may get suplicated columns but you can subset them out at your leisure
