# Project 3: Code 2

## Prepares Workspace

### Loads necessary libraries
```{r}
library(tidyverse)
```
### Reads in pdb file
```{r}
pdbFile <- readLines("C:/Users/Sternfeld/Desktop/ABT785/1ids.pdb")
```

## Generates data frame for alpha carbon atoms

### Converts pdb file into data frame of atoms using substring
> As the pdb file is a txt file with multiple information with specified spacing, substring can be used to extract the desired elements. Here, "ATOM" is searched for at the beginning of each line, specifying the presence of an atom and its associated information. A data frame is generated with the column names specified below.
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
> When the data frame above is made, substring does not know to make the column numeric. To do the math required later, these columns must be converted to numeric, as is done below.
```{r}
NumColumns <-c("Atom_serial_number", "Residue_sequence_number", "Temperature_factor", 
               "X_orthogonal_coordinate", "Y_orthogonal_coordinate", "Z_orthogonal_coordinate")
AtomDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(AtomDF[[x]]))
```

### Isolate Alpha Carbons and make a DF with their information
> Only alpha carbon atoms will be used for downstream analysis, so they are subset out of the main atom data frame
```{r}
CA_atoms <- subset(AtomDF, Atom_name =='CA')
CA_atoms <- CA_atoms[, c("Residue_sequence_number", "Chain", "X_orthogonal_coordinate", "Y_orthogonal_coordinate", "Z_orthogonal_coordinate")]
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

### Converts the following columns to numbers
> When the data frame above is made, substring does not know to make the column numeric. To do the math required later, these columns must be converted to numeric, as is done below.
```{r}
NumColumns <-c("ID#", "Start", "End")
secondaryDF[, NumColumns] <- lapply(NumColumns, function(x) as.numeric(secondaryDF[[x]]))
```

### Merges alpha carbon (CA) coordinates with secondary data frame
> Here the secondaryDF and CA_atoms data frames are merged so that the ***Start*** residue of the secondary structure pulls in the alpha carbon coordinates from that residue (and Chain). The columns are then labeled Start_X, Start_Y, and Start_Z, because next this resulting data frame is merged with the alpha carbon data frame again, but this time using the ***End*** residue of the secondary stucture (and Chain). These newly pulled in columns are labeled with End_X, End_Y, and End_Z. A new column is then made calculating the distance between the first and last residue of the secondary structure. In case there are no secondary structures, an if statement is used around this entire section of code, such that it will be skipped if none exist.
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

>Questions 2-4 are answered within this same code. An if statement is used to determine if there are alpha helices and/or beta sheets in protein, printing if only one or neither type exist in the protein. And then if one or both types exists, the data frame of secondary structures will be printed, which includes the distances between their first and last alpha carbon.
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
