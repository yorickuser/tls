## Test functions written by Hiroshi C. Ito (2024)
## Email: hiroshibeetle@gmail.com
##
## copyright (C) 2024 Hiroshi C. Ito
## This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License Version 2 as 
## published by the Free Software Foundation.
## http://www.r-project.org/Licenses/

library(treestats);
library(ape);
library(phytools);


gof <- function(){
    graphics.off();
}

flag_alive_only=0;
data_type=1;
runid=3;


odir0="tmp";
ofile_base="trial";

if(data_type==0)dtype="";
if(data_type==1)dtype="_discrete";
if(data_type==2)dtype="_sigmoid";

odir=paste0(odir0,dtype);
ofile=paste0(odir,"/",ofile_base,runid,".RData"); ## file name of simulation result
print(ofile);

load(ofile);

source("smvl_tls4e.R");
flag_check_tree2phylo=FALSE;


tree=adj_tree(a$res$tree_sort);



xran=max(abs(tree2traitu(tree,"x")));
nsamp=50;


tsamp=seq(a$timen*0.1,a$timen,,nsamp);


t_show=c(3,nsamp);


res=tree_sampling(tree,tsamp=tsamp,dp=0.0,out_analysis=TRUE,plot_each=FALSE);



tab=res$phe; 
res=res$res; 

gidu=sort(unique(tab$gid)); 
tuu=sort(unique(tab$t)); 


clade_size=matrix(rep(0,length(gidu)*length(tuu)),ncol=length(gidu))
clade_n=matrix(rep(0,length(gidu)*length(tuu)),ncol=length(gidu))
for(i in 1:length(tuu)){
    for(j in 1:length(gidu)){
        mask=((tab$t==tuu[i])*(tab$gid==gidu[j])>0);
        clade_size[i,j]=sum(mask);
        clade_n[i,j]=sum(tab$n[mask]);
    }
}



bgcol="black";
fgcol="white";
x11();par(bg=bgcol,fg=fgcol);

dev.hold();


plot(tuu,rowSums(clade_size),type="l",ylim=c(0,max(rowSums(clade_size))),col="gray",lty=2,fg=fgcol,col.axis=fgcol,col.lab=fgcol,col.main=fgcol,main="ex4: time changes of clade sizes",xlab="time",ylab="Clade size (number of phenotypes)");


cols=palfunc(list(gid=gidu));
for(i in 1:length(gidu))lines(tuu,clade_size[,i],col=cols[i]);

legend("topleft",legend=c("total",paste(gidu)),col=c("white",cols),lty=c(2,rep(1,length(gidu))),pch="",cex=1);
dev.flush();



x11();
par(bg=bgcol,fg=fgcol);
dev.hold();
plot(tuu,rowSums(clade_n),type="l",ylim=c(0,max(rowSums(clade_n))),col="gray",lty=2,fg=fgcol,col.axis=fgcol,col.lab=fgcol,col.main=fgcol,main="ex4: time changes of clade population",xlab="time",ylab="Clade population size (number of phenotypes)");


cols=palfunc(list(gid=gidu));
for(i in 1:length(gidu))lines(tuu,clade_n[,i],col=cols[i]);


legend("topleft",legend=c("total",paste(gidu)),col=c("white",cols),lty=c(2,rep(1,length(gidu))),pch="",cex=1);
dev.flush();



nphe=unlist(stree1(res,"nphe"));

nsum=unlist(stree1(res,"nsum"));

tspe=unlist(stree1(res,"tsamp"));


collspe=rep(0.0,length(res));
for(i in 1:length(collspe)){
        collspe[i]=colless_corr(res[[i]]$tree.phylo);
}



x11();
par0 = par();
mai = par()$mai;
mai[4]=mai[1];
par(mai=mai);

dev.hold();

plot(tspe,nphe,type="l",col="black",ylab="Number of lineages",xlab="Time",ylim=c(0,1.2*max(c(nphe,nsum))));
lines(tspe,nsum,col="black",lty=2);
par(new=T);

plot(tspe,collspe,type="l",ylim=c(0,1.0),col="blue",ylab="",xlab="",axes=FALSE,lwd=2);
axis(4);
mtext("Colless index (corrected)",side=4,line=2);

lines(c(1,1)*tsamp[t_show[1]],c(0,collspe[t_show[1]]),col="red",lty=2,lwd=2);
lines(c(1,1)*tsamp[t_show[2]],c(0,collspe[t_show[2]]),col="red",lty=1,lwd=1);

legend("topleft",legend=c("Number of lineages","Total population size","Colless index (corrected)"),col=c("black","black","blue"),lty=c(1,2,1),pch=-1,cex=1);
dev.flush();



plot_tsamp <- function(res,tagi){
    ted0=res[[tagi]]$tsamp;
    if(flag_alive_only==1){
        ca=res[[tagi]]$ca;
        tr2=res[[tagi]]$subtreea;
    }
    
    if(flag_alive_only==0){
        tr1=res[[tagi]]$subtree;
        stat1=tree2stat(tr1);
        maskr= as.numeric(stat1$ped-stat1$pst!=0);


        tr1=select_tree(tr1,which(maskr>0));


        tr2=adj_tree(trim_tree(tr1,tedge=a$pparam$tedge));
}



    res_phylo=tree2phylo(tr2,lab=NULL,pparam=a$pparam,timen=ted0,flag_check=flag_check_tree2phylo);

    trr1=res_phylo$tree.phylo;
    stat1=res_phylo$stat;
    stat=res_phylo$stat0;
    lab0=res_phylo$lab0;
    x=tree2trait_ed(tr2,"x");
    coll1=colless_corr(trr1);
    

    trr=res[[tagi]]$tree.phylo;
    coll=colless_corr(trr);


x11();
par(bg="#333333");
fgcol="#999999";
dev.hold();

plot(tab$x,tab$t,type="n",,cex.lab=1.2,xlim=c(-xran,xran),ylim=c(0,res[[tagi]]$tsamp*1.05),xlab="Trait x",ylab="Time",fg=fgcol,col.axis=fgcol,col.lab=fgcol,col.main=fgcol,main=sprintf("Time: %3.0f \nColless index (corrected): %2.3f",res[[tagi]]$tsamp,coll));

plot_tree(res[[tagi]]$subtree,col="black",lty=2,hold=1);
plot_tree(res[[tagi]]$subtreea,col=palfunc,hold=1);
ymax=res[[tagi]]$tsamp*1.05;

text(x,stat$ted+ymax*0.02,label=lab0,col=c("white","cyan")[(stat$alive>0)+1],cex=0.8,srt=90,pos=3);
dev.flush();


x11();
dev.hold();
plot.phylo(trr,direction="u",tip.color="blue",main=sprintf("phylo (alive lineages) \n Time: %3.0f \nColless index (corrected): %2.3f",res[[tagi]]$tsamp,coll));
axis(2, las = 1);
dev.flush();


x11();
dev.hold();
plot.phylo(trr1,direction="u",root.edge=T,tip.color=c("black","blue")[(res_phylo$stat$alive>0)+1],edge.color="black",main=sprintf("phylo (alive and dead lineages) \n Time: %3.0f \nColless index (corrected): %2.3f",res[[tagi]]$tsamp,coll1));

axis(2, las = 1);
dev.flush();
}


plot_tsamp(res,tagi=t_show[1]);
plot_tsamp(res,tagi=t_show[2]);

   


flag_out=0;
if(flag_out==1){
    
    ofile="tree_analysis1c";
    ofile1=ofile;
    if(data_type==1)ofile1=paste0(ofile,"_disc");
    if(data_type==-1)ofile1=paste0(ofile,"_sigm");
    
    for(i in dev.list()){
        geom=800;
        pngout(i,paste0(ofile1,"_",i),outeps=TRUE,geometry=geom);
    }

    
}
