# Project 3; Code 3

## PREPARES WORKSPACE

### loads necessary libraries
```{r}
library(tidyverse)
library(stringr)
```

### Reads in pdb file
```{r}
ClustalOutput <- readLines("C:/Users/Sternfeld/Desktop/ABT785/ClustalOutput.txt")
```

### Reads in groups
```{r}
GroupLabels <- read.csv("C:/Users/Sternfeld/Desktop/ABT785/uniprotkbOutput.csv")
```
### Make vector of AAs
```{r}
AAs <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-")
```

## Generates data frame of  IDs and sequences

### Loops through ClustalOutput and Generates a data frame with ID and sequence chunks
```{r}
columns <- c("ID", "Sequence")
df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(df) <- columns
for (i in 1:length(ClustalOutput)) {
  if (!str_extract_all(ClustalOutput[i], pattern = "\\|([^\\|]+)\\|")=="character(0)") {
    df[nrow(df) + 1,] <- c(gsub("\\|", "", str_extract_all(ClustalOutput[i], pattern = "\\|([^\\|]+)\\|")),
                           sub("^.{36}(\\S+)\\s(.*)$", "\\1",  ClustalOutput[i]))
  }
}
```
### Group the data frame by ID
```{r}
grouped_df <- group_by(df, ID)
```
### Concatenate Sequences within each group into single line
```{r}
concatenated_df <- summarize(grouped_df, Sequence = paste0(Sequence, collapse = ""))
```
### Convert to dataframe
```{r}
concatenated_df <- as.data.frame(concatenated_df)
```
## Generating a data frame of amino acid usage by residue number

### Identifies the length of the alignment
```{r}
AlignedLength <- nchar(concatenated_df[1,2])
```
### sets up dataframe
```{r}
columns <- AAs
AAusage <- data.frame(matrix(nrow = AlignedLength, ncol = length(columns)))
colnames(AAusage) <- columns
```

### loops through the number of aligned residues and totals the number of sequencings containing the specified amino acid
```{r}
for (i in 1:AlignedLength){
  AAinPos <- c()
  for (n in 1:length(columns)){
    AAinPos <- c(AAinPos, sum(substr(concatenated_df$Sequence, i, i) == colnames(AAusage)[n]))
  }
  AAusage[i,] <- AAinPos
}
```
## Generates data frame of amino acid use percentage by residue number per grouping

### Merge group label with already generated concatenated sequencing data frame
```{r}
concatenated_df <- left_join(concatenated_df, GroupLabels, by = c("ID" = "Entry"))
```
### Identify group names
```{r}
GroupNames <- unique(concatenated_df$Group)
```
### Create an empty list to store data frames
```{r}
GroupedDataframes <- list()
```
### Loop over the list of group names to make a data frame containing the percent contribution of an amino acid to its aligned residue number 
```{r}
columns <- c(AAs, "-")
for (d in GroupNames[1:length(GroupNames)]) {
  columns <- c(AAs, "-")
  temp <- data.frame(matrix(nrow = AlignedLength, ncol = length(columns)))
  colnames(temp) <- columns
  concatenated_df_subset <- subset(concatenated_df, Group == d)
  for (i in 1:AlignedLength){
    AAinPos <- c()
    for (n in 1:length(columns)){
      AAinPos <- c(AAinPos, sum(substr(concatenated_df_subset$Sequence, i, i) == colnames(temp)[n]))
    }
    temp[i,] <- (AAinPos/sum(AAinPos))*100
  }
  GroupedDataframes[[d]] <- temp
}
```

## Final Output

### Writes a csv containing the number of sequences that use the specified amino acid (columns) at each residue number (rows)
```{r}
write.csv(AAusage, "C:/Users/Sternfeld/Desktop/ABT785/AAusage.csv")
```
### Writes a csv for each grouping containing the percentage of sequences within that group use the specified amino acid (columns) at each residue number (rows)
```{r}
for (i in names(GroupedDataframes)) {
  write.csv(GroupedDataframes[[i]], paste0("C:/Users/Sternfeld/Desktop/ABT785/", i,".csv"))
}
```
