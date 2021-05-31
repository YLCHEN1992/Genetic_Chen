# Need packages
library(corrplot)
library(ggplot2)
library(formattable)
library(reshape2)

# Public functions

# Matrix for linkage calculation
Dmdata=c(1,0,0,0,0.5,0,0.5,0,0,0,0.5,0.5,0,0,0,1,
0,0.5,0,0.5,0,1,0,0,0.5,0,0,0.5,0.25,0.25,0.25,0.25,0,0.5,0.5,0)
Dm=matrix(Dmdata,nrow=9,ncol=4,byrow=T)
colnames(Dm)=c("g11","g22","g12","g21")
rownames(Dm)=c("T1111","T1112","T1122","T2211","T2212","T2222","T1211","T1212","T1222")

# Linkage calculation
LDm=function(x){
LD=c()
if(("-" %in% as.character(x[,2]))|("-" %in% as.character(x[,3]))){
for(i in 2:3){LD=c(LD,which(as.character(x[,i])=="-"))}
x=x[-unique(LD),]}
x[x=="21"]="12"
x[x!="11"&x!="12"&x!="22"]="12"
T4=c()
for(i in 1:nrow(x)){
t=as.numeric(Dm[which(rownames(Dm)==paste("T",as.character(x[i,2]),as.character(x[i,3]),sep="")),])
T4=c(T4,t)}
T4M=matrix(T4,ncol=4,byrow=T)
g11=sum(T4M[,1])/nrow(T4M)
g22=sum(T4M[,2])/nrow(T4M)
g12=sum(T4M[,3])/nrow(T4M)
g21=sum(T4M[,4])/nrow(T4M)
D=g11*g22-g12*g21    # CORE EQUATION
q1=genebf(as.character(x[,2]))[4]
p1=genebf(as.character(x[,2]))[5]
q2=genebf(as.character(x[,3]))[4]
p2=genebf(as.character(x[,3]))[5]
if(p1&p2&q1&q2==0){
R2=0
ZD=D}else{
if(D>0){ZD=D/min(q1*p2,q2*p1)
}else if(D<0){ZD=D/min(q1*p1,q2*p2)}else{ZD=D}
R2=D^2/(p1*p2*q1*q2)}
BDR=c(abs(ZD),R2)
BDR}

# Genetic linkage map display
LMAP=function(x,LD="D'"){
library(ggplot2)
library(reshape2)
cormat=round(x,2)
cormat[lower.tri(cormat)]=NA
melted_cormat= melt(cormat,na.rm=TRUE)
ggheatmap=ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
geom_tile(color = "white")+
scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
midpoint =0.5, limit = c(0,1), space = "Lab", 
name=paste("Linkage Disequilibrium\n",as.character(LD),sep="")) +
theme_minimal()+
theme(axis.text.x = element_text(angle = 45,vjust = 1,size = 12,hjust = 1))+
coord_fixed()
map=ggheatmap+
geom_text(aes(Var2, Var1, label = value), color = "black", size = 10)+
theme(
axis.title.x = element_blank(),
axis.title.y = element_blank(),
panel.grid.major = element_blank(),
panel.border = element_blank(),
panel.background = element_blank(),
axis.ticks = element_blank(),
legend.justification = c(1, 0),
legend.position = c(0.35, 0.8),
legend.direction = "horizontal")+
guides(fill = guide_colorbar(barwidth = 8, barheight = 1,
title.position = "top", title.hjust = 0.5))
map
}

# Values calculated and formated
genebf=function(x){
sx=as.character(x)[x!="-"]
sy=as.character(unlist(strsplit(sx,"")))
msx=table(sx)
msy=table(sy)
ntox=sum(msx)
ntoy=sum(msy)
if(is.na(as.numeric(msx["11"]))){msx["11"]=0}
if(is.na(as.numeric(msx["12"]))){msx["12"]=0}
if(is.na(as.numeric(msx["22"]))){msx["22"]=0}
if(is.na(as.numeric(msy["1"]))){msy["1"]=0}
if(is.na(as.numeric(msy["2"]))){msy["2"]=0}
nf11=as.numeric(msx["11"])
nf12=as.numeric(msx["12"])
nf22=as.numeric(msx["22"])
nf1=as.numeric(msy["1"])
nf2=as.numeric(msy["2"])
f11=nf11/ntox
f12=nf12/ntox
f22=nf22/ntox
f1=nf1/ntoy
f2=nf2/ntoy
acf11=f11*ntox
acf22=f22*ntox
acf12=f12*ntox
exf11=ntox*f1^2
exf22=ntox*f2^2
exf12=ntox*2*f1*f2
Ho=sum(nf11,nf22)/ntox              # CORE EQUATION
He=nf12/ntox                        # CORE EQUATION
Ne=1/Ho                             # CORE EQUATION
PIC=1-f1^2-f2^2-2*f2^2*f1^2         # CORE EQUATION
XF=((abs(acf11-exf11)-0.5)^2/(exf11))+
((abs(acf12-exf12)-0.5)^2/(exf12))+
((abs(acf22-exf22)-0.5)^2/(exf22))  # Chi-square value
PV=1-pchisq(XF,1)                   # Chi-square P value
BC=c(f11,f12,f22,f1,f2,XF,PV,nf11,nf12,nf22,nf1,nf2,ntox,ntoy,Ho,He,Ne,PIC)
BC}


# Main function 1 !!!

Genetic_Chen=function(gebd){
address=getwd()

	# Read and normalize file 
x=read.csv(deparse(substitute(gebd)))
sites=ncol(x)-1
frenq=c()
x[x=="21"]="12"
x[x!="11"&x!="12"&x!="22"]="12"

	# Multiple rule
for(i in 1:sites){
assign(paste("site",i,sep=""),as.character(x[,i+1]))
getn=paste("site",i,sep="")
frenq=c(frenq,genebf(get(getn)))} 
Mfrenq=matrix(frenq,nrow=sites,ncol=18,byrow=T) 

	# Values extract
gefbm=data.frame(Sitenames=as.character(colnames(x)[-1]),
FrequenceOF11=Mfrenq[,1],NumberOF11=Mfrenq[,8],
FrequenceOF12=Mfrenq[,2],NumberOF12=Mfrenq[,9],
FrequenceOF22=Mfrenq[,3],NumberOF22=Mfrenq[,10],
FrequenceOF1=Mfrenq[,4],NumberOF1=Mfrenq[,11],
FrequenceOF2=Mfrenq[,5],NumberOF2=Mfrenq[,12],
NumberOFSample=Mfrenq[,13],
NumberOFGene=Mfrenq[,14],
X_Statistics=Mfrenq[,6],P_Value=Mfrenq[,7],
Homozygosity=Mfrenq[,15],Heterozygosity=Mfrenq[,16],Ne=Mfrenq[,17],PIC=Mfrenq[,18])

	# Files save
if (file.exists("./Rgenetics")==TRUE){cat("阁下目标文件夹 Rgenetics 已存在\n")}else{
dir.create("./Rgenetics", recursive=TRUE)
cat("目标文件夹 Rgenetics 已为阁下创建\n")}
setwd("./Rgenetics")
NAME=paste("阁下遗传统计已计算完成",gsub(":","_",Sys.time()),".csv")
write.csv(gefbm,NAME,row.names=FALSE)
cat("阁下基础遗传数据分析已完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)

	# Web display format 
Wgefbm=cbind(Sitenames=gefbm[,1],round(gefbm[,-1],3))
pp= formatter("span", 
style = x ~ style(
font.weight = "bold", 
color = ifelse(x > 0.05, "Green", ifelse(x < 0.05, "Red", "black"))))
nn=formatter("span", 
style =~style(
color ="grey",font.weight = "bold")) 
FQ=color_tile("MediumAquamarine","MediumAquamarine")
webtable=formattable(Wgefbm, align =c("l",rep("c",17)),list('P_Value' =pp,
'Sitenames' =nn,
'FrequenceOF11'=FQ,
'FrequenceOF12'=FQ,
'FrequenceOF22'=FQ,
'FrequenceOF1'=FQ,
'FrequenceOF2'=FQ))

	# Jugement of linkage 
if(sites>1){
LDDM=matrix(0,nrow=sites,ncol=sites)
LDRM=matrix(0,nrow=sites,ncol=sites)
for(a in 2:ncol(x)){
for(b in 2:ncol(x)){
LDDM[a-1,b-1]=LDm(x[,c(1,a,b)])[1]
LDRM[a-1,b-1]=LDm(x[,c(1,a,b)])[2]
if(a==b){
LDDM[a-1,b-1]=1
LDRM[a-1,b-1]=1}}}
colnames(LDDM)=as.character(colnames(x)[-1])
rownames(LDDM)=as.character(colnames(x)[-1])
colnames(LDRM)=as.character(colnames(x)[-1])
rownames(LDRM)=as.character(colnames(x)[-1])
mapLDDM=LMAP(LDDM)
mapLDRM=LMAP(LDRM,LD="R2")
if (file.exists("./Rgenetics")==TRUE){cat("阁下目标文件夹 Rgenetics 已存在\n")}else{
dir.create("./Rgenetics", recursive=TRUE)
cat("目标文件夹 Rgenetics 已为阁下创建\n")}
setwd("./Rgenetics")
gNAMED=paste("阁下遗传连锁图D绘制已完成",gsub(":","_",Sys.time()),".png")
gNAMER=paste("阁下遗传连锁图R绘制已完成",gsub(":","_",Sys.time()),".png")
ggsave(filename=gNAMED,mapLDDM,dpi=600,width=8,height=8)
ggsave(filename=gNAMER,mapLDRM,dpi=600,width=8,height=8)
cat("阁下连锁遗传数据分析已完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)}
webtable}

# Search Haplotype    1    Generate group seeds
rephap=function(x){
t=x
for(i in 1:length(x)){
for(j in 1:2){
x=c(x,paste(as.character(x[i]),as.character(j),sep=""))}}
s=setdiff(x,t)
s}
# Search Haplotype    2    Generate groups seeds*
seedhap=function(n){
po=c()
if(n==1){
po=c("1","2")}else{
po=c("1","2")
for (a in 1:(n-1)){
po=rephap(po)}}
po}
# Search Haplotype    3     I am an genius
pxy=function(x,y){
if(x=="11"&y=="1"){pr=1}
if(x=="11"&y=="2"){pr=0}
if(x=="22"&y=="1"){pr=0}
if(x=="22"&y=="2"){pr=1}
if(x=="12"&y=="1"){pr=0.5}
if(x=="12"&y=="2"){pr=0.5}
pr}

# Main function 2
HapChen=function(gebd){
address=getwd()

	# Read and normalize file 
x=read.csv(deparse(substitute(gebd)))
LD=c()
if(("-" %in% as.character(x[,2]))|("-" %in% as.character(x[,3]))){
for(i in 2:3){LD=c(LD,which(as.character(x[,i])=="-"))}
x=x[-unique(LD),]}
x[x=="21"]="12"
x[x!="11"&x!="12"&x!="22"]="12"
sites=ncol(x)-1
	# Generate groups seeds
nhap=length(seedhap(sites))
pnxxs=c()
	# Calculate Haplotypes
for(j in 1:nhap){
pnxx=c()
for(a in 1:nrow(x)){
pn=c()
px=""
py=""
for(i in 1:sites){
px=substring(paste(as.character(x[a,2:ncol(x)]),collapse = ""),(i*2-1),i*2)
py=substring(as.character(seedhap(sites)[j]),i,i) 
pn=c(pn,pxy(px,py))}     
pnxx=c(pnxx,prod(pn))}  # I am an genius
cat("统计单倍型",seedhap(sites)[j],"完成\n")
pnxxs=c(pnxxs,(sum(pnxx)/nrow(x)))}

	# Haplotype Save
mhapy=data.frame(Haplotype=seedhap(sites),FreHaplotype=round(pnxxs,3))
if (file.exists("./Rgenetics")==TRUE){cat("阁下目标文件夹 Rgenetics 已存在\n")}else{
dir.create("./Rgenetics", recursive=TRUE)
cat("目标文件夹 Rgenetics 已为阁下创建\n")}
setwd("./Rgenetics")
NAME=paste("阁下单倍型-遗传统计已计算完成",gsub(":","_",Sys.time()),".csv")
write.csv(mhapy,NAME,row.names=FALSE)
cat("阁下遗传数据分析已完成,文件保存在",as.character(getwd()),"目录下\n")
setwd(address)

	# Haplotype Show
pp= formatter("span", 
style =  ~ style(
font.weight = "bold", 
color ="Cornislk"))
HYOtable=formattable(mhapy, align =c("l","c"),
list('FreHaplotype' =pp,'Haplotype' =pp))
HYOtable}






