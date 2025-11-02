## Test functions written by Hiroshi C. Ito (2024)
## Email: hiroshibeetle@gmail.com
##
## copyright (C) 2024 Hiroshi C. Ito
## This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License Version 2 as 
## published by the Free Software Foundation.
## http://www.r-project.org/Licenses/

source("smvl_tls4e.R");

show_sample_plot=FALSE; 

out_file_RData=TRUE; 
file_RData="ddflow_simulation.RData"; 

niche_dim=1;
sK=1.0; 
sa=0.2; 
b=0.1; 


flag_disc=TRUE;


m_rate_x=1.0; 
sm_x=0.01;
sm_y=2.0; 
m_rate_y=0.000025 

m_rate=m_rate_x
if(flag_disc)m_rate=m_rate_x+m_rate_y;
##m_rate_y=m_rate_x*sm_x^2/sm_y^2;




halt_edge_y=10.0; 
##halt_edge_y=3.0; 
oint=1000000000; 
##oint=100; 

pparam=pparam0; 

pparam$xid=1; 
pparam$yid=2; 
pparam$tedge=10e4; 
pparam$traj_line_color=T; 

pparam$bgcol="#FFFFFF"; 
pparam$fgcol="#000000"; 

if(pparam$traj_line_color){ 
    pparam$bgcol="#333333";
    pparam$fgcol="#999999";
}

pparam$winh=7;
pparam$winw=6;
pparam$winh_sub=5;
pparam$winw_sub=9;

pparam$omit_plot=TRUE;
##pparam$omit_plot=FALSE;

runid=1; 
seed=23;

##set.seed(19);


reset_win=T;


flag_show_divtime_only_alive=TRUE; 
flag_show_tedge=TRUE; 
file_data="out.txt";
mydir="work/"; 


if(niche_dim==1)phe_init=list(x=0.5,y=0.0,gid=0); 
if(niche_dim==2)phe_init=list(x=0.5,w=0.1,y=0.0,gid=0); 
n_init=0.9; 

q0 <- function(phe1,sigma){
    if(niche_dim==1)dis=((phe1$x)^2);
    if(niche_dim==2)dis=((phe1$x)^2+(phe1$w)^2);
    return(exp(-dis/(2*sigma^2)));
}

q <- function(phe1,phe,sigma){
    if(niche_dim==1)dis=((phe$x-phe1$x)^2);
    if(niche_dim==2)dis=((phe$x-phe1$x)^2+(phe$w-phe1$w)^2);
    return(exp(-dis/(2*sigma^2)));
}


qq <- function(phe,sigma){
    if(niche_dim==1)dis=(as.matrix(dist(phe$x)))^2;
    if(niche_dim==2)dis=(as.matrix(dist(phe$x)))^2+(as.matrix(dist(phe$w)))^2;
    return(exp(-dis/(2*sigma^2)));
}


K <- function(phe){
    return(q0(phe,sK));  
}

alpha <- function(phe1,phe){
    return(q(phe1,phe,sa));
}



growth_rate <- function(phe1){
    return(sqrt(b*phe1$y+1));
}





fitness0 <-function(phe1,phe,n){ 
    return( growth_rate(phe1) - sum(alpha(phe1,phe)*n)/K(phe1) );
}


fitness=fitness0;


set_parms <-function(phe,n){
    return(list(phe=phe,Kp=K(phe),ap=qq(phe,sa)));
}

pop_dynamics0 <- function(t,n,parms){
    list(n*(growth_rate(parms$phe)-rowSums(parms$ap%*%n)/parms$Kp));     
}


pop_dynamics=pop_dynamics0;


mutate1d <- function(phe){
    a$flag_new_group <<- 0;
        mgid=phe$gid;
        mx=phe$x+rnorm(1,mean=0.0,sd=sm_x);
        my=phe$y+rnorm(1,mean=0.0,sd=sm_x);
        return(list(x=mx,y=my,gid=mgid));
}

mutate_disc1d <- function(phe){
    a$flag_new_group <<- 0;
    if(runif(1)<m_rate_x/(m_rate_x+m_rate_y)){
        mx=phe$x+rnorm(1,mean=0.0,sd=sm_x);
        if(niche_dim==2)mw=phe$w+rnorm(1,mean=0.0,sd=sm_x);
        my=phe$y;
        mgid=phe$gid;
    }else{
        a$flag_new_group <<- 1;
        mgid=max(a$tree_phe$gid)+1;
        mx=phe$x;
        if(niche_dim==2)mw=phe$w;
        my=phe$y+rnorm(1,mean=0.0,sd=sm_y);    
    }
    return(list(x=mx,y=my,gid=mgid));
    
}


mutate2d <- function(phe){
    a$flag_new_group <<- 0;
        mgid=phe$gid;
        mx=phe$x+rnorm(1,mean=0.0,sd=sm_x);
        mw=phe$w+rnorm(1,mean=0.0,sd=sm_x);
        my=phe$y+rnorm(1,mean=0.0,sd=sm_x);
        return(list(x=mx,w=mw,y=my,gid=mgid));
}


mutate_disc2d <- function(phe){
    a$flag_new_group <<- 0;
    if(runif(1)<m_rate_x/(m_rate_x+m_rate_y)){
            mx=phe$x+rnorm(1,mean=0.0,sd=sm_x);
            mw=phe$w+rnorm(1,mean=0.0,sd=sm_x);
            my=phe$y;
            mgid=phe$gid;
        }else{
            a$flag_new_group <<- 1;
            mgid=max(a$tree_phe$gid)+1;
            mx=phe$x;
            mw=phe$w;
            my=phe$y+rnorm(1,mean=0.0,sd=sm_y);    
        }
    return(list(x=mx,w=mw,y=my,gid=mgid));
    
}


output <- function(timen,phe,n){
        
    options(scipen=100);
    
    divt=rep(0,length(phe$x));
    divy=rep(0,length(phe$x));
    divx=rep(0,length(phe$x));
    if(length(n)>1){
        if(length(a$res$divt)==length(n))divt=a$res$divt;
        if(length(a$res$divy)==length(n))divy=a$res$divy;
        if(length(a$res$divx)==length(n))divx=a$res$divx;
    }
    
    ym=sum(n*phe$y)/sum(n);
    ymax=max(phe$y);


    if(a$ninv==0)cat("ninv", 
                     "pid", 
                     "t",   
                     "n",  
                     "x",  
                     "y",  
                     "divt", 
                     "divx", 
                     "divy", 
                     "gid", 
                     "ym",  
                     "ymax",
                     "\n",file=a$file_data,append=TRUE);
    for(i in 1:length(phe$x)){
        cat(a$ninv,phe$pid[i],timen,a$n[i],phe$x[i],phe$y[i],divt[i],divx[i],divy[i],phe$gid[i],ym,ymax,"\n",file=a$file_data,append=TRUE);
        
    }
    options(scipen=0);
}


halt_func <- function(){
    yg=sum((a$phe$y)*(a$n))/sum(a$n);
    if(yg>halt_edge_y)a$sparam$flag_halt <<-TRUE;

    if(a$ninv%%500==0){
        cat("runid:",a$runid,"time:",a$timen, "ymax:",max(a$phe$y)," residents:",length(a$n), " invasion:", a$ninv,"amp:",a$sparam$amp_invf,"fit_over:",length(a$sparam$fit_over),"\n");
        
    }
}


gof <- function(){
    graphics.off();
}

show_comwin=0;
if(show_comwin==1)comwin(1,mydir=mydir);


if(niche_dim==1){
    mutate1=mutate1d;
    if(flag_disc)mutate1=mutate_disc1d;
}
if(niche_dim==2){
    mutate1=mutate2d;
    if(flag_disc)mutate1=mutate_disc2d;
}



arg = commandArgs(trailingOnly=TRUE);
print(arg);
if(length(arg)>0){
    cat("arguments:\n");
for(i in 1:length(arg)){
    cat(arg[i],"\n");
    eval(parse(text=arg[i]));
}
}


set.seed(seed);

time_st = proc.time()[3];



simevol(phe_init, 
        n_init, 
        fitness, 
        pop_dynamics=pop_dynamics, 
        set_parms=set_parms, 
        runid=runid,
        fitness_contour = F, 
        plot_func=plot_func00, 
        mutate=mutate1, 
        m_rate= m_rate, 
        bgid=3,
        pparam=pparam, 
        palfunc=palfunc,
        halt_func=halt_func,
        show_interval=oint,
        out_interval=oint,
        file_data=file_data,
        output=output, 
        reset_win=reset_win,
        mydir=mydir); 


if(pparam$omit_plot)analysis_and_plot(adj_tree(a$tree),pparam$xid,pparam$yid,pparam,analyze=TRUE);

time_ed = proc.time()[3];

cat("\n Time taken for simulation:",time_ed-time_st,"\n");


if(out_file_RData){
    print("saving the simulation data...");
    ##save(list="a",file=file_RData);
    save(list=ls(),file=file_RData);
    cat("Data was saved: \n",file_RData,"\n");
}



ex0 <- function(){
    p=pparam;
    x11(bg="black"); 

    x=tree2traitu(a$res$tree2,"x"); 
    ranx=max(abs(x)); 
    xx=seq(-ranx,ranx,,100); 
    KK=K(list(x=xx));

    
    plot(xx,KK*0.5*(a$timen),ylim=c(0,a$timen),type="l",col="green",lty=2, main="ex0",xlab="Trait x", ylab="Time",fg=p$fgcol,col.axis=p$fgcol,col.lab=p$fgcol,col.main=p$fgcol);

    
    plotmain(xid=1,yid=0,show_bran_alive=F,show_bran=F,show_tdiv=F,hold=TRUE);
    
}
    

ex1 <- function(){
    p=pparam;
    x11();
    dev.hold();
    
    res=a$res; 
    x=tree2traitu(res$tree2,"x"); 
    ranx=max(abs(x)); 

    
    
   
    plot(c(-ranx,ranx),c(0,a$timen),type="n",main="ex1",xlab="Trait x (niche trait)", ylab="Time");
    rect(-ranx*2,-1*a$timen,ranx*2,a$timen*1.1,col="#000033"); 
    
    plot_tree(res$tree2,xid=1,yid=0,col="red",lty=2,xlim=c(-1,1)*ranx,hold=TRUE); 
    plot_tree(res$tree2a,xid=1,yid=0,lwd=2,col=palfunc,lty=1,hold=TRUE); 
    points(res$ext$x,res$ext$t,pch=4,col="red"); 
    points(res$bran$x,res$bran$t,pch=16,col="orange",cex=0.7); 
    points(res$brana$x,res$brana$t,pch=16,col="white"); 

   
    phe=a$phe;
    n=a$n;
    divt=a$res$divt;
    for(i in 1:length(phe$x))lines(phe$x[i]*c(1,1),c(divt[i],a$timen),col="white",lty=2);


    points(phe$x,phe$x+a$timen,pch=16,col=palfunc(phe),cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
    dev.flush();
}


ex2 <-function(){
    x11();
    tab=read.table(file_data,header=T); 

    xran=max(abs(tab$x)); 
    divt1=tab$t-tab$divt; 
    xx=tab$x; 

   
    divtm=histm(xx,divt1,lev=seq(-xran,xran,,11)); 

    
    plot(xx,divt1,xlim=c(-1,1)*xran,ylim=c(0,1.5*max(divtm$val)),pch=3,col="#FFAAAA",cex=0.8,main="ex2: Divergence time",xlab="Trait x",ylab="Time");

   
    lines(divtm$lev,divtm$val,type="b",lwd=2,col="red");

    x1=seq(-xran,xran,,100);
    KK=K(list(x=x1));
    vy_max_est=0.5*sm_x^2*(0.5*b)/sqrt(2*pi);
    dh=-1*sign(x1)*sqrt(1/(KK+1e-10)-1);
    h=cumsum(dh)*(x1[2]-x1[1]);
    h=h-max(h);
    divt_est=abs(h)/vy_max_est;
    lines(x1,divt_est,lty=2,col="black");
    
    legend("topright",legend=c("divtime depth (mean)","divtime depth (mean)","prediction"),col=c("red","red","black"),lty=c(0,1,2),pch=c("+","o",""),cex=1);
    
}


ex3 <- function(){
    x11();
    cols=c("blue","red","purple");
    bran=a$res$bran; 
    ext=a$res$ext;  

    xran=max(abs(a$tree_phe$x))*1.2; 
    
    
    hb=hist(bran$x,breaks=seq(-xran,xran,,10),plot=F); 
    he=hist(ext$x,breaks=seq(-xran,xran,,10),plot=F);  
    yb=hb$counts;
    ye=he$counts;
    ybe=yb-ye;
    x=hb$mids;

    x1=seq(-xran,xran,,100);
    KK=K(list(x=x1));
    ybe_est=sqrt(KK/(1-KK+1e-10))*abs(x1/sK^2)*(KK-0.5);
    ybe_est=(max(ybe)/max(ybe_est))*ybe_est;

    ymax=max(c(yb,ye,ybe,ybe_est));
    ymin=min(c(yb,ye,ybe,ybe_est));


    plot(x,yb,xlim=c(-xran,xran),ylim=c(ymin,ymax),type="b",pch=16,col=cols[1],main="ex3: Branching and extinction",xlab="trait x",ylab="Frequency");
    lines(x,ye,col=cols[2],type="b",pch=16); 
    lines(x,ybe,col=cols[3],type="b",pch=16); 
    lines(x,ybe*0,col="black",lty=2); 

    lines(x1,ybe_est,col=cols[3],lty=2);
    legend("topright",legend=c("branching","extinction","net-branching","prediction(scaled)"),col=cols[c(1,2,3,3)],lty=c(1,1,1,2),pch="",cex=1);

}


ex4 <- function(){

    bgcol="black";
    fgcol="white";
    x11();par(bg=bgcol,fg=fgcol);

    tab=read.table(file_data,header=T); 
    
    gidu=sort(unique(tab$gid)); 
    
    ninvu=sort(unique(tab$ninv)); 
    tu=tab$t[match(ninvu,tab$ninv)]; 

    
    clade_size=matrix(rep(0,length(gidu)*length(ninvu)),ncol=length(gidu))
    for(i in 1:length(ninvu))for(j in 1:length(gidu))clade_size[i,j]=sum((tab$ninv==ninvu[i])*(tab$gid==gidu[j]));

    
    plot(tu,rowSums(clade_size),type="l",col="gray",lty=2,fg=fgcol,col.axis=fgcol,col.lab=fgcol,col.main=fgcol,main="ex4: time changes of clade sizes",xlab="time",ylab="Clade size (number of phenotypes)");

    
    cols=palfunc(list(gid=gidu));
    for(i in 1:length(gidu))lines(tu,clade_size[,i],col=cols[i]);
    
    legend("topleft",legend=c("total",paste(gidu)),col=c("white",cols),lty=c(2,rep(1,length(gidu))),pch="",cex=1);
    
}


ex5 <- function(){

    x11();
    tree2=a$res$tree2;
    
    vmx=tree2velocity(tree2,"x"); 
    vmy=tree2velocity(tree2,"y"); 
    
    xran=max(abs(vmx$trait));
    
    
    hdx=histsum(vmx$trait,vmx$dif,lev=seq(-xran,xran,,11));
    hdy=histsum(vmx$trait,vmy$dif,lev=seq(-xran,xran,,11));
    hdt=histsum(vmx$trait,vmx$tdif,lev=seq(-xran,xran,,11));

    vx=hdx$val/hdt$val;
    vy=hdy$val/hdt$val;
    
    xx=hdx$lev;

    
    plot(xx,vy,type="n",ylim=c(min(c(vx,-vy)),max(c(vx,vy))),xlim=c(-xran,xran),xlab="Trait x",ylab="Average velocity of directional evolution",main="ex5: Average velocity of directional evolution for each niche");
    
    lis=seq(1,length(vmx$trait),10);
    #points(vmx$trait[lis],vmy$vel[lis],pch=3,col="#FFAAAA",cex=0.4)
    #points(vmx$trait[lis],vmx$vel[lis],pch=3,col="#AAAAFF",cex=0.4)

    
    x1=seq(-xran,xran,,100);
    KK=K(list(x=x1));
    vy_max_est=0.5*sm_x^2*(0.5*b)/sqrt(2*pi); 
    
    
    lines(x1,vy_max_est*KK,col="red",lty=2);

    
    lines(x1,vy_max_est*sign(x1)*sqrt(KK*(1-KK)),col="blue",lty=2);

    
    lines(xx,vx,type="b",col="blue",lwd=2);
    lines(xx,vy,type="b",col="red",lwd=2);
    lines(c(-xran,xran),c(0,0),col="black",lty=2);    
    
    legend("topleft",legend=c("Trait x","Trait y","prediction","prediction"),col=c("blue","red","blue","red"),lty=c(1,1,2,2),pch=c("o","o","",""),cex=1);

    }



if(show_sample_plot==TRUE){
    ex0();
    ex1();
    ex2();
    ex3();
    ex4();
    ex5();
}




flag_out=0;
if(flag_out==1){
    
    ofile="ddflow_analysis3b1";
    
    for(i in dev.list()){
        geom=800;
        if(i==3)geom=1200;
        pngout(i,paste0(ofile,"_",i),outeps=TRUE,geometry=geom);
}


}








