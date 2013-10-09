setwd("/Users/rohosenb/Documents/UCDavis/STA-250/HW0/HW0_Q3")      #setting directory so that new files would be saved in desired location
string = "Hello, my name is Bob. I am a statistician. I like statistics very much."
l = nchar(string)           #finding out length of my string

#creating new fles
for(i in 1:l) write(substring(string,i,i),paste("char",i,".txt"))

#reading the files to get back the original string
new = NULL
for(i in 1:l){
  char = scan(paste("char",i,".txt"),what=character(),sep="$")
  new = paste(new,char,sep="")
}
