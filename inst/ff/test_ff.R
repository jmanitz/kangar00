##################################
##  test ff on simulation data  ##
##################################


library("ff")
library("data.table")


path.sim <- "Q:/1_PROJEKTE/kangar00-paket/kangar00/inst/data/"
             

### data:
           
   geno   <- fread(paste(path.sim,"sim2.geno.txt",sep=""), sep='auto', 
             header=FALSE, data.table=FALSE)  # reads 200 rows + 116417 columns
   rownames(geno) <- geno[,1]
   geno <- geno[,-c(1,116417)] #has column of NAs at the end, I don't know why.
                               # = 116415 'real' columns of genotypes 
   
              
### dataframe to ff-dataframe:

   geno.ffdf <- as.ffdf(geno[1:100,1:100]) #yes, works (on windows and linux)

   geno.ffdf <- as.ffdf("geno = some large dataframe") #don't run!!!
   #windows: 
   # is not finished after 24 h. It fills the whole memory, however the process 
   # does not seem to stay running this long since the processor is only only at 
   # 1-2% capacity. Memory is not released when R is terminated. The last time I
   # tested this I could not even start the task manager to stop the process 
   # and free the memory. Only solution was to turn off the whole computer.
   #linux:
   # does not work, gives following error: 
   # 'Error in ff(initdata = initdata, length = length, levels = levels, ordered
   #  = ordered, : unable to open [alternativ auch] write error
   # 

   #save(geno.ffdf, file="geno.ffdf.RData")
   test.1 <- get(load("geno.ffdf.RData"))  
   #windows + linux:
   # can be saved and loaded, even after I copied it to a different folder
 
     
     

### matrix to ff - (function as.ff() does not work on a dataframe):

   geno.ff  <- as.ff(as.matrix(geno[1:100,1:100])) #works fine (windows + linux)

   geno.ff  <- as.ff(as.matrix(geno))
   #windows:
   # data cannot be transformed into matrix due to out of memory error when 
   # creating the matrix. 
   #linux:
   # works. 

   #save(geno.ff, file="geno.ff.RData")
   test.2 <- get(load("geno.ff.RData"))  
   #windows + linux:
   # can be saved and loaded, even after I copied it to a different folder
   
    
    