# Project 3: Code 1

## Prepares Workspace

### Loads necessary libraries
```{r}
library(tidyverse)
```
### Reads in pdb file
```{r}
pdbFile <- readLines("C:/Users/Sternfeld/Desktop/ABT785/1ids.pdb")
```

## Generates data frame for atoms

### Converts pdb file into data frame of atoms using substring
> As the pdb file is a txt file with multiple information with specified spacing, substring can be used to extract the desired elements. Here, "ATOM" is searched for at the beginning of each line, specifying the presence of an atom and its associated information. A data frame is generated with the column names specified below.
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

### Converts columns to numerical value
> When the data frame above is made, substring does not know to make the column numeric. To do the math required later, these columns must be converted to numeric, as is done below.
```{r}
NumColumns <-c("Atom_serial_number", "Residue_sequence_number", "Temperature_factor")
AtomDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(AtomDF[[x]]))
```


## Generates data frame for secondary structures

### Converts pdb file into data frame of secondary structures using substring
> The pdb file also contains information for the secondary structures alpha helices and beta sheets. This information is placed into a data frame using substring, similarly to how the atom data frame was generated.
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

### Converts columns to numerical value
> When the data frame above is made, substring does not know to make the column numeric. To do the math required later, these columns must be converted to numeric, as is done below.
```{r}
NumColumns <-c("ID#", "Start", "End")
secondaryDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(secondaryDF[[x]]))
```

### Sort by Chain and Residue number 
> For reasons that will be pointed about below, this data frame must be sorted by Chain and Start residue.
```{r}
secondaryDF <- secondaryDF %>% arrange(Chain, Start)
```


## Generates dataframe for loop atoms

### Make empty dataframe for loop atoms
```{r}
LoopAtoms <- data.frame(matrix(nrow = 0, ncol = length(columns)))
```
### Identify Loop # and their associated atoms
> As a loop is identified as amino acids between two secondary structures, for each Chain, with the secondaryDF structured above (sorted by Chain and then Start residue), we can simply find all atoms from residues that fall between the ***End*** residue of one secondary structure (n) and the ***Start*** of the next secondary structure (n+1).
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
### Adds a Chain-Loop column
> To aid in the coding below (Final math and output section), which asks for the average of the individual loops within each chain, I made a "ChainLoop" column for the data frame.
```{r}
LoopAtoms$ChainLoop <- paste0(LoopAtoms$Chain, "-", LoopAtoms$LoopNumber)
```


## Final math and output
### Question #2 
> Prints the average temperature (B) factor value for atoms in each loop
```{r}
for (i in unique(LoopAtoms$ChainLoop)){
  print(paste0("The average temperature for Loop ",unique(LoopAtoms[LoopAtoms$ChainLoop == i, "LoopNumber"]),
               " in Chain ",unique(LoopAtoms[LoopAtoms$ChainLoop == i, "Chain"])," is ",round(mean(LoopAtoms[LoopAtoms$ChainLoop == i, "Temperature_factor"]),3)))
}
```

### Question #3
> Writes a csv file containing the atoms within loop regions with their associated residue number and name, among additional information
```{r}
write.csv(LoopAtoms, "C:/Users/Sternfeld/Desktop/ABT785/LoopAtoms.csv")
```
