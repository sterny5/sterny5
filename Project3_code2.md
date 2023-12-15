# Project 3; Code 2
## PREPARES WORKSPACE

### loads necessary libraries
```{r}
library(tidyverse)
```
### Reads in pdb file
```{r}
pdbFile <- readLines("C:/Users/Sternfeld/Desktop/ABT785/1ids.pdb")
```
## Generated data frame of alpha carbons

### uses substring to collect atom elements to add to dataframe
```{r}
columns <- c("Atom_serial_number", "Atom_name", "Residue_sequence_number", "Residue_name","Chain", "Temperature_factor", 
             "X_orthogonal_coordinate", "Y_orthogonal_coordinate", "Z_orthogonal_coordinate")
AtomDF <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(AtomDF) <- columns
for (i in 1:length(pdbFile)) {
  if (startsWith(pdbFile[i], "ATOM")) {
    AtomDF[nrow(AtomDF) + 1,] <- c(substr(pdbFile[i], 7, 11), str_trim(substr(pdbFile[i], 13, 16)), 
                           substr(pdbFile[i], 23, 26), str_trim(substr(pdbFile[i], 18, 20)), 
                           str_trim(substr(pdbFile[i], 22, 22)), substr(pdbFile[i], 61, 66), substr(pdbFile[i], 31, 38), 
                           substr(pdbFile[i], 39, 46), substr(pdbFile[i], 47, 54))
  }
}
```

### Converts columns to numerical value
```{r}
NumColumns <-c("Atom_serial_number", "Residue_sequence_number", "Temperature_factor", 
               "X_orthogonal_coordinate", "Y_orthogonal_coordinate", "Z_orthogonal_coordinate")
AtomDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(AtomDF[[x]]))
```

### Isolate Alpha Carbons and make a DF with their information
```{r}
CA_atoms <- subset(AtomDF, Atom_name =='CA')
CA_atoms <- CA_atoms[, c("Residue_sequence_number", "Chain", "X_orthogonal_coordinate", "Y_orthogonal_coordinate", "Z_orthogonal_coordinate")]
```

## Generates data frame for secondary structures

### Uses substring to identify secondary structures to add to dataframe
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

### Converts the following columns to numbers
```{r}
NumColumns <-c("ID#", "Start", "End")
secondaryDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(secondaryDF[[x]]))
```

### merges alpha carbon (CA) coords into DF listing first and last residues
> if statement is used in case there are no alpha helices or beta sheets
```{r}
if(length(secondaryDF$Type)>0){
  secondaryDF <- left_join(secondaryDF, CA_atoms, by = c("Start" = "Residue_sequence_number", "Chain"))
  colnames(secondaryDF) <- c("Type", "Chain", "ID#", "Start", "End", "Start_X", "Start_Y", "Start_Z")
  secondaryDF <- left_join(secondaryDF, CA_atoms, by = c("End" = "Residue_sequence_number", "Chain"))
  colnames(secondaryDF) <- c("Type", "Chain", "ID#", "Start", "End", "Start_X", "Start_Y", "Start_Z", "End_X", "End_Y", "End_Z")
  secondaryDF[,3:11] <- secondaryDF[,3:11] %>% mutate_if(is.character, as.numeric)
  
  #Calculates distance between start and end of secondary structures and appends to "secondaryDF"
  secondaryDF <- mutate(
    secondaryDF,
    Distance = sqrt((secondaryDF$End_X - secondaryDF$Start_X)^2 + 
                    (secondaryDF$End_Y - secondaryDF$Start_Y)^2 + 
                    (secondaryDF$End_Z - secondaryDF$Start_Z)^2))
}
```


## Final Output

>Questions 2-4 are answered within this same code. Determines if there are alpha helices and/or beta sheets in protein and then prints distances between their start and end, if applicable
```{r}
if ("HELIX" %in% unique(secondaryDF$Type) & "SHEET" %in% unique(secondaryDF$Type)) {
  print(secondaryDF)
} else if("HELIX" %in% unique(secondaryDF$Type)){
    print("There are no beta sheets")
    print(secondaryDF)
    } else if("SHEET" %in% unique(secondaryDF$Type)){
    print("There are no alpha helices")
    print(secondaryDF)
    } else print("There are no beta sheets or alpha helices")
```
