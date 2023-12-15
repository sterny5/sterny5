# Project 3: Code 3

## Prepares Workspace

### Loads necessary libraries
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
> To be used later, a vector of amino acids by single letter codes and a gap ("-") is generated below.
```{r}
AAs <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "-")
```

## Generates data frame of  IDs and sequences

### Loops through ClustalOutput and Generates a data frame with ID and sequence chunks
> As the ClustalOutput file is a txt file with  information denoted between specialized characters or spacing, we can extract the information needed below (ID and Sequence) as specified below.
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
> Due to the nature of the ClustalOutput file, only a certain amount of sequence for each ID/organism will be each line, such that if 20 ID/organisms are tested, the first chunk of sequence will appear for ID# 1 on line one and then the second chunk of sequence for ID #1 will appear on line #21. To get all of the sequence for each ID/organism on one line, they are first grouped and then the sequences are concatenated and place in a new line of a data frame, as shown below.
```{r}
grouped_df <- group_by(df, ID)
concatenated_df <- summarize(grouped_df, Sequence = paste0(Sequence, collapse = ""))
concatenated_df <- as.data.frame(concatenated_df)
```
## Generating a data frame of amino acid usage by residue number
> To answer Question 2, my goal was to generate a data frame that had the total number of aligned residues as rows and the possible amino acids (and gap) in columns. Then, for at each residue (row), for all IDs/organisms tested, the number of sequences that use each amino acids (columns) were counted and reported.

### Identifies the length of the alignment
> As the alignment length was the same for all IDs, any row could have been investigated, here, the sequence length (character count in Column 2) from the ID/organism in Row #1 was tested.
```{r}
AlignedLength <- nchar(concatenated_df[1,2])
```
### Sets up data frame
```{r}
columns <- AAs
AAusage <- data.frame(matrix(nrow = AlignedLength, ncol = length(columns)))
colnames(AAusage) <- columns
```

### Creates amino acid by residue counts table/data frame for all organisms/IDs 
> This code uses two for loops. It loops through the number of aligned residues (i: 1-209 in this case), and for each residue, it loops through each amino acid (n), totalling the number of sequencings containing the specified amino acid and building them into a data frame (AAusage).
```{r}
for (i in 1:AlignedLength){
  AAinPos <- c()
  for (n in 1:length(columns)){
    AAinPos <- c(AAinPos, sum(substr(concatenated_df$Sequence, i, i) == colnames(AAusage)[n]))
  }
  AAusage[i,] <- AAinPos
}
```
## Generates data frame of amino acid use across all organisms/IDs
> This is very similar code to above, however, here, instead of obtaining the ***total number*** of each type of amino acid at a residue position for ***all*** organisms/IDs, it calculates the ***percentage*** of each amino acid found at each residue ***within groups*** of organisms/IDs.

### Merge group label with already generated concatenated sequencing data frame
> In order to group the sequences appropriately, the concatenated_df is merged with the GroupLabels.
```{r}
concatenated_df <- left_join(concatenated_df, GroupLabels, by = c("ID" = "Entry"))
```
### Identify group names
> This list will be used to loop through the data.
```{r}
GroupNames <- unique(concatenated_df$Group)
```
### Create an empty list to store data frames
```{r}
GroupedDataframes <- list()
```
### Generates data frame of amino acid use percentage by residue number per grouping
> First a temporary data frame is generated to subset by group. With each grouping separated, as done above, the total number of residues containing a specific amino acid is calculated. To obtain a percentage, this number is then divided by the total number of residues analyzed at that position and multiplied by 100.
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

### Writes a csv containing the ***total number*** of sequences that use the specified amino acid (columns) at each residue number (rows)
```{r}
write.csv(AAusage, "C:/Users/Sternfeld/Desktop/ABT785/AAusage.csv")
```
### Writes a csv for ***each grouping*** containing the ***percentage*** of sequences within that group use the specified amino acid (columns) at each residue number (rows)
```{r}
for (i in names(GroupedDataframes)) {
  write.csv(GroupedDataframes[[i]], paste0("C:/Users/Sternfeld/Desktop/ABT785/", i,".csv"))
}
```
