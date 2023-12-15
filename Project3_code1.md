# PREPARES WORKSPACE

## loads necessary libraries
```{r}
library(tidyverse)
```

## Reads in pdb file
```{r}
pdbFile <- readLines("C:/Users/Sternfeld/Desktop/ABT785/1ids.pdb")
```


# GENERATES DATAFRAME FOR ATOMS

## uses substring to collect atom elements to add to dataframe
```{r}
columns <- c("Atom_serial_number", "Atom_name", "Residue_sequence_number", "Residue_name", "Chain", "Temperature_factor")
AtomDF <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(AtomDF) <- columns
for (i in 1:length(pdbFile)) {
  if (startsWith(pdbFile[i], "ATOM")) {
    AtomDF[nrow(AtomDF) + 1,] <- c(substr(pdbFile[i], 7, 11), str_trim(substr(pdbFile[i], 13, 16)), 
                                   substr(pdbFile[i], 23, 26), str_trim(substr(pdbFile[i], 18, 20)), 
                                   str_trim(substr(pdbFile[i], 22, 22)), substr(pdbFile[i], 61, 66))
  }
} 
```

## Converts columns to numerical value
```{r}
NumColumns <-c("Atom_serial_number", "Residue_sequence_number", "Temperature_factor")
AtomDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(AtomDF[[x]]))
```


# GENERATES DATAFRAME FOR SECONDARY STRUCTURES

## uses substring to identify secondary structures to add to dataframe
```{r}
columns <- c("Type", "Chain", "ID#", "Start", "End")
secondaryDF <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(secondaryDF) <- columns
for (i in 1:length(pdbFile)) {
  if (startsWith(pdbFile[i], "HELIX")) {
    secondaryDF[nrow(secondaryDF) + 1,] <- c("HELIX", str_trim(substr(pdbFile[i], 20, 20)), substr(pdbFile[i], 8, 10), 
                                             substr(pdbFile[i], 22, 25), substr(pdbFile[i], 34, 37))
  }
  if (startsWith(pdbFile[i], "SHEET")) {
    secondaryDF[nrow(secondaryDF) + 1,] <- c("SHEET", str_trim(substr(pdbFile[i], 22, 22)), substr(pdbFile[i], 8, 10), 
                                             substr(pdbFile[i], 23, 26), substr(pdbFile[i], 34, 37))
  }
}
```

## Converts columns to numerical value
```{r}
NumColumns <-c("ID#", "Start", "End")
secondaryDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(secondaryDF[[x]]))
```

## Sort by Chain and Residue number 
```{r}
secondaryDF <- secondaryDF %>% arrange(Chain, Start)
```


# GENERATES DATAFRAME FOR LOOP ATOMS

## Make empty dataframe for loop atoms
```{r}
LoopAtoms <- data.frame(matrix(nrow = 0, ncol = length(columns)))
```
## Identify Loop # and their associated atoms
```{r}
for (i in LETTERS[1:4]) {
  #subset df by Chain
  temp <- secondaryDF[secondaryDF$Chain == i, ]
  #loop through non-secondary residues to identify loops
  for (n in 1:(nrow(temp)-1)){
    temp2<-subset(AtomDF,
                  AtomDF$Chain == i &
                    AtomDF$Residue_sequence_number > as.numeric(secondaryDF[n, "End"]) &
                    AtomDF$Residue_sequence_number < as.numeric(secondaryDF[n+1, "Start"]))
    #assigns LoopNumber to atoms
    temp2 <- cbind(temp2, data.frame(LoopNumber = rep(n, nrow(temp2))))
    LoopAtoms <- rbind(LoopAtoms, temp2)
  }
}
```
## Adds a Chain-Loop column
```{r}
LoopAtoms$ChainLoop <- paste0(LoopAtoms$Chain, "-", LoopAtoms$LoopNumber)
```


# FINAL MATH AND OUTPUT
## QUESTION 2 
> Prints the average temperature (B) factor value for atoms in each loop
```{r}
for (i in unique(LoopAtoms$ChainLoop)){
  print(paste0("The average temperature for Loop ",unique(LoopAtoms[LoopAtoms$ChainLoop == i, "LoopNumber"]),
               " in Chain ",unique(LoopAtoms[LoopAtoms$ChainLoop == i, "Chain"])," is ",round(mean(LoopAtoms[LoopAtoms$ChainLoop == i, "Temperature_factor"]),3)))
}
```

## QUESTION 3
> Writes a csv file containing the atoms within loop regions with their associated residue number and name, among additional information
```{r}
write.csv(LoopAtoms, "C:/Users/Sternfeld/Desktop/ABT785/LoopAtoms.csv")
```
