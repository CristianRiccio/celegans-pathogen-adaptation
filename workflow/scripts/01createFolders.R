library(yaml)
foldersListOfLists <- yaml::yaml.load_file('input/folders.yaml')

flattenList = function(listOfLists){
  if (class(listOfLists) == 'list')
    lapply(listOfLists, flattenList)
  else if (class(listOfLists) != 'list')
    TRUE
}

longNames <- names(unlist(lapply(foldersListOfLists, flattenList)))
folders <- gsub('.', '/', longNames, fixed = TRUE)

# http://hydroecology.net/iterating-through-lists-of-lists-of-lists/
for (folder in folders) {
  if (!file.exists(folder)) dir.create(folder, recursive = TRUE)
}
