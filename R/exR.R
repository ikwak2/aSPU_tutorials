library(aSPU)
g <- "AGRN"

## load variance among phenotypes matrix
varPhe <- as.matrix(read.table("exvarPhe.txt"))

tof <- paste("subdat",g,".txt",sep ="")
cmd2run <- paste("awk 'NR==FNR {h[$1] = $1; next} h[$1] !=\"\" {print $0}' exGenes/",g,".txt exdat.txt >", tof, sep="")
system(cmd2run)

## laod Summary statistics data for gene AGRN
sdat <- read.table(file=tof, colClasses=c("character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "character", "character")  )
system(paste("rm",tof) )


## get corresponding data from reference panel
snpids <- sdat[,1]
tempf = paste("tempids",g,".txt" , sep="")
write.table(snpids,  file = tempf, row.names=FALSE, col.names=FALSE, quote=FALSE )
tof <- paste("tempdata",g,".txt", sep="")

system( paste("awk 'NR==FNR {h[$1] = $0; next} h[$1] !=\"\" {print h[$1]}' SubsetCh1.txt ", tempf," > ", tof, sep="") )

locdat <- read.table(file=tof)
datn <- t(matrix(as.numeric(as.matrix(locdat[,3:381])), dim(locdat)[1],379))

colnames(datn) <- locdat[,1]


## Save used allele for 0,1,2 coding
usedallele <- as.character(locdat[,2])
system(paste("rm",tof))

snpids1 <- snpids[(snpids %in% colnames(datn))]


## Erase column that are all 0's 1's and 2's

a = max(apply( datn == 0, 2, sum )) == 379
b = max(apply( datn == 1, 2, sum )) == 379
c = max(apply( datn == 2, 2, sum )) == 379

if( a | b | c ) {
    eraseit = 0
    eraseit <- c(eraseit, which( apply( datn == 0, 2, sum ) == 112 ),
                 which(apply( datn == 1, 2, sum ) == 112),
                 which(apply( datn == 2, 2, sum ) == 112) )

    datn <- datn[, -eraseit]
    snpids1 <- snpids1[-eraseit]
    usedallele <- usedallele[-eraseit]
}


## Alter 0 and 2 if different alleles are used for 0,1,2 coding
rownames(sdat) <- sdat[,1]
allele <- toupper(sdat[snpids1, 11])
difal <- ifelse( allele == usedallele, 1, -1)
nsdat <- sdat[snpids1,]

subdat2 <- as.matrix(datn)
temp <- subdat2[, difal == -1 ]
temp[ temp == 0] = 3
temp[ temp == 2] = 0
temp[ temp == 3] = 2
subdat2[, difal == -1 ]  <- temp


## Calculate Correlation
datn <- as.matrix(subdat2)
datn <- as.matrix(datn)
corSNP <- cor(as.matrix(datn), use="pairwise.complete.obs")

Ps <- as.matrix(nsdat[,c(4:9)])


## Erase SNPs with correlation > .95 and < -.95
while( sum( corSNP < -0.95 ) > 0 ) {
    toerase <- which( apply( corSNP < -.95, 1, sum) == max( apply( corSNP < -.95, 1, sum) ) )
    corSNP <- corSNP[-toerase[1], -toerase[1]]
    Ps <- Ps[-toerase[1],]
    if( is.matrix(corSNP) == FALSE) break
}

corSNP = as.matrix(corSNP)
while( sum( corSNP > 0.95 ) > ncol(corSNP) ) {
    toerase <- which( apply( corSNP > 0.95, 1, sum) == max( apply( corSNP > 0.95, 1, sum) ) )
    corSNP <- corSNP[-toerase[1], -toerase[1]]
    Ps <- Ps[-toerase[1],]
    if( is.matrix(corSNP) == FALSE) break
}

Ps = as.matrix(Ps)

## Perform MTaSPUsSet test.
## Use MTaSPUsSet() function if n.perm > 10^6 otherwise MTaSPUsSetC() is faster.

##    out <- MTaSPUsSet(Ps, corSNP=corSNP, corPhe = varPhe, pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=10^3, Ps=TRUE)
out <- MTaSPUsSetC(Ps, corSNP=corSNP, corPhe = varPhe, pow=c(1,2,4,8),  pow2 = c(1,2,4,8), n.perm=10^3, Ps=TRUE)

## Perform aSPUs function for each trait.
Out <-NULL
for(i in 1:6) {
    Out[[i]] <- aSPUs(Ps[,i], corrSNP = corSNP, pow=c(1,2,4,8), n.perm = 100, Ps=TRUE)$pvs
}
Out
