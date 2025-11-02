# copyright (C) 2020 Hiroshi C. Ito
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License Version 2 as 
# published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available at
# http://www.r-project.org/Licenses/

#' Simulator of adaptive evolution.
#' @aliases simevol simevol-package
#' @keywords internal
"_PACKAGE"



library(deSolve);
##library(envstocker);

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/  Core function for plotting and analysis   _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##


analysis_and_plot <- function(tree1,xid,yid,p,analyze=TRUE){           
##.ee.append("analysis_and_plot",environment())
########### Analysis ######################
    if(analyze){
        ## a$tree is composed of many short branches, and which makes analysis
        ## on the tree structure difficult. Hence, we first reorganize it so that
        ## the tree is composed of small number of long branches and many short ones.
        ## Each long branch corresponds to each lineage's directional evolution.
        
        print("tree sorting...");
        cat("number of branches:",length(tree1)," ")
        tree1b=sort_tree(a$tree);

    cat("\n");print("tree reducing...");
    tree_sort=reduce_tree(tree1b);


    ## Since the reorganized tree has many short branches, and which can cause
    ## overestimates of branching and extinction, we remove the short branches
    ## by setting a threshold time length for them "tedge".

    ## triming the tree composed of both alive and dead lineages
    print("tree triming...");
    tree2=adj_tree(trim_tree(tree_sort,tedge=p$tedge));

    ## triming the tree composed of only alive lineages
    print("tree triming...");    
    tree2a=adj_tree(trim_tree(tree_sort,only_alive=TRUE));
    print("done...");
    

    ## set colors for branches based on their group ID
    if(p$traj_line_color){
        cols = .simevol_func$palfunc(tree2a);        
        for(kk in 1:length(tree2a)){            
            tree2a[[kk]]$col=cols[[kk]];
        }
    }else{
        for(kk in 1:length(tree2a))tree2a[[kk]]$col="blue";
        
    }

    ## store the reorganized trees in a list "red"
    res=list(tree_sort = tree_sort);
    res$tree2 = tree2;
        res$tree2a = tree2a;
        res$tree1 = tree1;

  
    ## calulate divergence times among the alive phenotypes
    ##if(flag_divtime_only_alive){
    
        ca2a=find_common_ans(tree2a,p);
        ##cpid=div_time(tree2a);
        xx2a=tree2trait_ed(tree2a,xid);
        yy2a=tree2trait_ed(tree2a,yid);
        al2a=tree2alive(tree2a);
        ppid=al2a$pid;
        
    ##}else{
    ca2=find_common_ans(tree2,p);
        ##cpid=div_time(tree2);
    xx2=tree2trait_ed(tree2,xid);
    yy2=tree2trait_ed(tree2,yid);
    al2=tree2alive(tree2);

   ## print("al2");print(al2);
    ##ppid=al$pid;
        
    ##}
    
##    ca$tdivp=ppid;
    
    tdiv=ca2a$tdiv;
    tdivy=ca2a$tdivy;
    tdivx=ca2a$tdivx;

       
    ## find branchings and extinctions by analyzing tree2 and tree2a

    bran=tree2phe_st(tree2);
        brana=tree2phe_st(tree2a);
        

   ## ext=tree2phe_last(a$tree2);
   ## lis=which(ext0$pid==-1);
   ## for(kk in 1:length(ext0))ext[[kk]]=ext0[[kk]][lis];

    
    al22=tree2stat(tree2);
    lis=which(al22$alive==-1);
    ext=tree2phe_ed(tree2);
    
    for(kk in 1:length(ext))ext[[kk]]=ext[[kk]][lis];
        ext$pid=al22$ped[lis];

        ##print(ppid);
    ##print(a$phe$pid);
    idd=match(a$phe$pid,ppid);
    
    ##print(idd);

    ## store these results in "res"
    res$divt = tdiv[idd];
    res$divy = tdivy[idd];
    res$divx = tdivx[idd];

    
    res$common_ans2a = ca2a;

    res$common_ans2 = ca2;
    
    res$bran = bran;
    res$ext = ext;
        res$brana = brana;
    
        res$p=p;
    a$res <<- res;
    a$res$ninv <<- a$ninv;
    }


########### Plotting ######################    
   if(!analyze){        
        res=a$res;
   }


if(!p$omit_plot){
plotmain(xid,yid,res,p,hold=TRUE);

dev.flush();
}

}

plotmain <- function(xid,yid,res=a$res,p=a$res$p,hold=FALSE,show_phe=TRUE,show_tree1=TRUE,show_tree2a=TRUE,show_tedge=flag_show_tedge, show_tdiv=TRUE,show_bran=TRUE,show_bran_alive=TRUE,show_ext=TRUE){
    if(!hold)plot_lim(xid,yid,p);
##.ee.append("plotmain",environment())    
    tree2a=res$tree2a;
    tree2=res$tree2;
    tree1=res$tree1;    
    ca2a=res$common_ans2a;
    ca2=res$common_ans2;
    bran=res$bran;
    brana=res$brana;
        ext=res$ext;

        
    xx2a=tree2trait_ed(tree2a,xid);
    yy2a=tree2trait_ed(tree2a,yid);
    al2a=tree2alive(tree2a);
    ppid=al2a$pid;
        

    xx2=tree2trait_ed(tree2,xid);
    yy2=tree2trait_ed(tree2,yid);
    al2=tree2alive(tree2);        


    xid1=xid;yid1=yid;
    if(xid==0)xid1=a$pdim+1;
    if(yid==0)yid1=a$pdim+1;

    branx=bran[[xid1]];
    brany=bran[[yid1]];
  
    
    branxa=brana[[xid1]];
    branya=brana[[yid1]];

    extx=ext[[xid1]];
    exty=ext[[yid1]];

    
    p$traj$col="orange";    
    if(show_tree1)plot_traj_line1(tree1,xid,yid,p);

    
    if(p$traj_line_color){
        lwd_tmp=p$traj$lwd;
        p$traj$lwd=2;
    }
    p$traj$col="blue"; 
    if(show_tree2a)plot_traj_line1(tree2a,xid,yid,p);
    p$traj$col="blue";
    
    
    
    if(p$traj_line_color)p$traj$lwd=lwd_tmp;
    
    col_ext="black";
    if(p$traj_line_color)col_ext="white";
    
    if(show_ext)points(extx,exty,col=col_ext,cex=0.7,pch=4);          
    if(show_bran)points(branx,brany,col="orange",cex=0.7,pch=16);
    if(show_bran_alive)points(branxa,branya,col="red",cex=1.0,pch=16);

    if(flag_show_divtime_only_alive){
        tdiv3=ca2a$tdiv;
        tdivy3=ca2a$tdivy;
        tdivx3=ca2a$tdivx;
        xx3=xx2a;
        yy3=yy2a;
    }else{
        tdiv3=ca2$tdiv;
        tdivy3=ca2$tdivy;
        tdivx3=ca2$tdivx;
        xx3=xx2;
        yy3=yy2;

    }
        
##    if(!analyze)print(ca);

    if(show_tdiv){
        if(yid==0){           
            ##print(tdiv3);
         
            for(kk in 1:length(tdiv3))lines(xx3[kk]*c(1,1),c(yy3[kk],tdiv3[kk]),col="red",lty=2);
            
            if(show_tedge)lines(xx3,xx3*0+a$timen-p$tedge,col="#999999",lty=2);
        }
        if(yid==2){           
            
            for(kk in 1:length(tdivy3))lines(xx3[kk]*c(1,1),c(max(yy3),tdivy3[kk]),col="red",lty=2);
        }
    }
    
    if(show_phe){
        if((xid>0)&&(yid>0)&&p$fitness_contour)plot_fitness_contour(land,p);
        ##plot_phe1(xid,yid,p);
        plot_phe1(xid,yid,p,phe=a$phe,n=a$n,timen=a$timen);
    }
}


## palette function used for coloring the tree branches and presently existing phenotypes
palfunc <- function(tree,uniq=NULL){
    npal=5;
    pal=(rev(rainbow(npal,start=0.1,s=0.5,end=0.72)));
    pal=pal[c(1,4,2,5,3)]

    if(class(tree[[1]])=="list"){
        
        gids=stree1(tree,"gid");
        ## print(gids);
        cids=gids;
        for(i in 1:length(gids)){
            if(length(uniq)==0)uniq=TRUE;
            if(uniq)gid1=unique(gids[[i]]);
            if(!uniq)gid1=(gids[[i]]);
            cid1=(gid1)%%(npal) +1;
            cids[[i]]=pal[cid1];
        }
        return(cids);
    }else{
            ##print(tree);
            if(length(uniq)==0)uniq=FALSE;
            
            if(uniq)gid=unique(tree$gid);
            if(!uniq)gid=(tree$gid);
            cid=(gid)%%(npal) +1;
            return(pal[cid]);
            
    }
    
        
}
    




##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/    Plotting function for ddflow analysis   _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
plotf <- function(out=FALSE,analyze=FALSE){
    plot_func00(out=out,analyze=analyze);
}

## the plotting function for the main and sub plotting windows
plot_func00 <- function(traj_line=TRUE,out=FALSE,analyze=FALSE){
    phe=a$phe;
    en=a$en;
    n=a$n;
    edim=a$edim;
    traj=a$traj;
    tree_phe=a$tree_phe;
    p=a$pparam;
    nspe=length(phe[[1]]);
    timen=a$timen;
    tedge=p$tedge;
    
    dev.set(a$winid[1]);
    par(bg=p$bgcol);
    par(fg=p$fgcol);
    dev.hold();
    
    xid=p$xid;yid=p$yid;
    trait_names=p$trait_names;
    env_names=p$env_names;
    plot_mask=p$plot_mask;

    
    
    if(class(.simevol_func$palfunc)=="function"){
        p$resi$bg=.simevol_func$palfunc(phe);
    }else{
        if(p$resi$bgid>-1){
            bgid=p$resi$bgid;
            npal=length(p$pal);
            if(bgid==0){
                bgval=n;
            }else{
                bgval=phe[[bgid]];
            }
            cmax=max(bgval)+1e-10;
            cmin=min(bgval)-1e-10;
            cid=as.integer((npal-1)*(bgval-cmin)/(cmax-cmin))+1;
            p$resi$bg=p$pal[cid];
        }
    }
        
    if(length(trait_names)==0)trait_names=names(phe);
    if((length(env_names)==0)&&(edim>0))env_names=paste0("Env",seq(edim));
        
    if((length(phe)-2)==1){
        plot_1dim(p);
        
    }else{        
       plot_lim(xid,yid,p);
       if((xid>0)&&(yid>0)&&p$fitness_contour){
           ranx=max(1.2*max(abs(tree_phe[[xid]])),0.001);        
           rany=max(1.2*max(abs(tree_phe[[yid]])),0.001);        
           
           land=calc_fitness_land(phe,en,xid,yid,xmin=-ranx,xmax=ranx,ymin=-rany,ymax=rany);
       }            
       tree1=adj_tree(a$tree);
       if(p$traj_line_color){
           
           if(class(.simevol_func$palfunc)=="function"){
               
           
               cols=.simevol_func$palfunc(tree1);
               for(kk in 1:length(tree1)){                       
                   tree1[[kk]]$col=cols[[kk]];
               }
                       
                   }else{
                       for(kk in 1:length(tree1)){                       
                           myx=tree1[[kk]][[bgid]];
                           colid=as.integer((npal-1)*(myx[length(myx)]-cmin)/(cmax-cmin))+1;
                           tree1[[kk]]$col=p$pal[colid];
                       }
                   }
       }
       
       analysis_and_plot(tree1,xid,yid,p,analyze=analyze);
    }
       
    if(a$show_subwin==TRUE){
        dev.set(a$winid[2]);
        par(bg=p$bgcol);
        par(fg=p$fgcol);
        dev.hold();
        
        for(i in 1:(length(phe)-2)){
            plot_lim(i,0,p);
            if(traj_line==TRUE){
                if((length(phe)-2)==1)tree1=adj_tree(a$tree);
                plot_traj_line1(tree1,i,0,p);
            }
            else{
                plot_traj(i,0,p);
            }
            
            
            plot_phe(i,0,p);
            
        }
        if(edim>0){
            for(i in 1:edim){
                plot(traj$e[,i],traj$te,col=p$env$col,lwd=p$env$lwd,type="l",xlab=env_names[i],ylab="Time",cex.lab=p$cex.lab);  
            }
        }
        
        dev.flush();        
    }
    
    if(out).simevol_func$output(a$timen,phe,n);
}



histm <- function(x,w,lev=NULL,nlev=10){
    if(length(lev)==0){
        lev=seq(min(x),max(x),nlev);
    }
    levd=diff(lev)*0.5;
    levd=c(levd[1],levd);
    lev1=c(lev-levd,lev[length(lev)]+levd[length(levd)]);        
    
    wm=rep(0,length(lev));
    mask0=rep(TRUE,length(lev));
    for(i in 1:length(lev)){
        mask= ((x>lev1[i])*(x<=lev1[i+1]))>0;
        if(sum(mask)>0){
            wm[i]=sum(w[mask])/sum(mask);
        }else{
            mask0[i]=FALSE
        }
    }
    return(list(val=wm[mask0],lev=lev[mask0]));
}

histsum <- function(x,w,lev=NULL,nlev=10){
    if(length(lev)==0){
        lev=seq(min(x),max(x),nlev);
    }
    levd=diff(lev)*0.5;
    levd=c(levd[1],levd);
    lev1=c(lev-levd,lev[length(lev)]+levd[length(levd)]);        
    
    wm=rep(0,length(lev));
    mask0=rep(TRUE,length(lev));
    for(i in 1:length(lev)){
        mask= ((x>lev1[i])*(x<=lev1[i+1]))>0;
        if(sum(mask)>0){
            wm[i]=sum(w[mask]);
        }else{
            mask0[i]=FALSE
        }
    }
    return(list(val=wm[mask0],lev=lev[mask0]));
}


##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/    Functions for tree structure adjustment _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##



find_ans <- function(mypid,treep,mask){
    found=NULL;
    fpid=NULL;
    ftid=NULL;

    lis=which(mask);
    for(j in lis){
        ans0=which(treep[[j]]==mypid);
        if(length(ans0)>0){
            ans0=ans0[1];
            if(ans0>1){
                ## cat("\nfind pid:", mypid, "found:", j, " ", ans0,"\n");
                found=c(j,ans0);
                ##  cat("\nfind:",mypid," ");
                ##  print(found);
                
                return(found);
            }
        }
    }
        return(NULL);
}


sort_tree2 <- function(tree){
    
    tree_adj=adj_tree(tree);
    tt=unlist(lapply(stree1(tree_adj,"t"),function(x){tail(x,n=1)}));
    
    if(length(tree)>1){
        
        al=tree2alive(tree);
        alt=tree2alive(tree_adj);
        
        pidlen=tree2lastpid(tree,len=TRUE);
        pidlast=pidlen$pid;
        nlen=pidlen$len;
        lis=order(-1*((tt+max(tt)*(al$alive>0))));
    
        tree0=tree;
        for(i in 1:length(tree)){
            tree[[i]]=tree0[[lis[i]]];
        }                                   
        pidlast=pidlast[lis];
        nlen=nlen[lis];
        
        treep=stree1(tree,"pid");
        pst=lapply(treep,function(x){x[1];});
        
        tree1=NULL;
        for(i in 1:length(tree))tree1 = c(tree1,list(list(pid=c(-100))));

        tree1b=pforeach(i=1:length(tree1), .combine='c',.errorhandling="stop")({    
           
            

        ##for(i in 1:length(tree1)){
        ##cat(i," ");
        tree1[[i]]=tree[[i]];

        mypid0=tree1[[i]]$pid[1];
        mypid=mypid0[1];
        mask= !((pidlast <= mypid)+(pst >= mypid0[length(mypid)])==2);
        ans=find_ans(mypid,treep,mask);
        
        while((length(ans)>0)){
            ##cat("i:", i, " pid:",tree1[[i]]$pid[1]);
            ##cat("find:",tree1[[i]]$pid[1]," ",ans,"\n");
            treej=getii(tree[[ans[1]]],st=1,ed=(ans[2]-1));
            
            tree1[[i]]=add_phenotype(treej,tree1[[i]]);
            
            mypid0=tree1[[i]]$pid[1];
            mypid=mypid0[1];
            mask= !((pidlast <= mypid)+(pst >= mypid0[length(mypid)])==2);
            
            ans=find_ans(mypid,treep,mask);
            
            pids=unlist(stree1(tree1,"pid"));
            dup=sum(pids==tree1[[i]]$pid[1]);
            
            if(dup>1)break;
            
        }
            
        
            list(tree1[[i]]);
        });
        tree1=tree1b;
        
    }else{
        tree1=tree;
    }
    
    return(tree1);
}


sort_tree <- function(tree){

    tree_adj=adj_tree(tree);
    tt=unlist(lapply(stree1(tree_adj,"t"),function(x){tail(x,n=1)}));
    
    if(length(tree)>1){
        
        al=tree2alive(tree);
        alt=tree2alive(tree_adj);
        
        pidlen=tree2lastpid(tree,len=TRUE);
        pidlast=pidlen$pid;
        nlen=pidlen$len;
        lis=order(-1*((tt+max(tt)*(al$alive>0))));
    
        tree0=tree;
        for(i in 1:length(tree)){
            tree[[i]]=tree0[[lis[i]]];
        }                                   
        pidlast=pidlast[lis];
        nlen=nlen[lis];
        
        treep=stree1(tree,"pid");
        pst=lapply(treep,function(x){x[1];});
        
        tree1=NULL;
        for(i in 1:length(tree))tree1 = c(tree1,list(list(pid=c(-100))));

        
        for(i in 1:length(tree1)){
            if(i%%100==0)cat(paste0(as.integer(100*i/length(tree1)),"% "));
        tree1[[i]]=tree[[i]];

        mypid0=tree1[[i]]$pid[1];
        mypid=mypid0[1];
        mask= !((pidlast <= mypid)+(pst >= mypid0[length(mypid)])==2);
            ans=find_ans(mypid,treep,mask);
            treepbuf=stree1(tree1,"pid");
            treepbuf[[i]]=c(-100);
            pids=unlist(treepbuf);
        while((length(ans)>0)){
            ##cat("i:", i, " pid:",tree1[[i]]$pid[1]);
            ##cat("find:",tree1[[i]]$pid[1]," ",ans,"\n");
            treej=getii(tree[[ans[1]]],st=1,ed=(ans[2]-1));

            bbb=add_phenotype(treej,tree1[[i]]);
            tree1[[i]]=bbb;
            
            mypid0=tree1[[i]]$pid[1];
            mypid=mypid0[1];
            mask= !((pidlast <= mypid)+(pst >= mypid0[length(mypid)])==2);
            ##mask[1:i]=FALSE;
            ans=find_ans(mypid,treep,mask);
            
            ##pids=unlist(stree1(tree1,"pid"));
            ##dup=sum(pids==tree1[[i]]$pid[1]);
            ##if(dup>1)break;
            dup=sum(pids==bbb$pid[1]);
            if(dup>0)break;

            
        }
        
        }
    }else{
        tree1=tree;
    }
    
    return(tree1);
}




reduce_tree <- function(tree1){
    rem=rep(0,length(tree1));
    nlen=length(tree1);
   if(length(tree1)>1){
       
        for(i in 1:(length(tree1)-1)){
            if(i%%100==0)cat(paste0(as.integer(100*i/length(tree1)),"% "));
         pidi=tree1[[i]]$pid;
            for(j in (i+1):length(tree1)){
                
                pidj=tree1[[j]]$pid;  
                
                if(!(pidj[1]>pidi[length(pidi)])||(pidj[length(pidj)]<pidi[1])){

                    if(sum(pidj==pidi[1])+sum(pidi==pidj[1])>0){
                    ##li=which(pidj==-1);
                    ##if(length(li)>0)pidj[li]=pidj[li]*0-100;
                    
                    if(pidj[length(pidj)]==-1)pidj[length(pidj)]=-100;
                    mask=rep(1,length(pidj));
                    mask[!(is.na(match(pidj,pidi)))]=0;
                    ##print(pidi);
                    ##print(pidj);
                    ##print(mask);
                
                lis=which(mask==1);
                
                
                    if(length(lis)>0){
                        tree1[[j]]=getii(tree1[[j]],st=min(lis)-1,ed=max(lis));
                    }else{
                        rem[j]=1;
                        cat("\n",j," th branch is removed!!\n");
                    }
                }
            }
            
            }
        }
   }

    tree1b=NULL;
    lis=which(rem==0);
    for(i in 1:length(lis)){
        tree1b = c(tree1b,list(tree1[[lis[i]]]));
    }
    
    
    return(tree1b);
    
}


reduce_tree_orig <- function(tree1){
    rem=rep(0,length(tree1));
    nlen=length(tree1);
   if(length(tree1)>1){
       
        for(i in 1:(length(tree1)-1)){
            if(i%%100==0)cat(paste0(as.integer(100*i/length(tree1)),"% "));
         pidi=tree1[[i]]$pid;
            for(j in (i+1):length(tree1)){
                
                pidj=tree1[[j]]$pid;  
                
                if(!(pidj[1]>pidi[length(pidi)])||(pidj[length(pidj)]<pidi[1])){

                    li=which(pidj==-1);
                    if(length(li)>0)pidj[li]=pidj[li]*0-100;
                    mask=rep(1,length(pidj));
                    mask[which(match(pidj,pidi)>0)]=0;
                    ##print(pidi);
                    ##print(pidj);
                    ##print(mask);
                
                lis=which(mask==1);
                
                
                    if(length(lis)>0){
                        tree1[[j]]=getii(tree1[[j]],st=min(lis)-1,ed=max(lis));
                    }else{
                        rem[j]=1;
                        cat("\n",j," th branch is removed!!\n");
                    }
                }
            }
            
        }
   }

    tree1b=NULL;
    lis=which(rem==0);
    for(i in 1:length(lis)){
        tree1b = c(tree1b,list(tree1[[lis[i]]]));
    }
    
    
    return(tree1b);
    
}


trim_tree <- function(tree2,edge=3,tedge=-1,omit_alive=FALSE,only_alive=FALSE){
    
    if(length(tree2)>1){
        al=tree2stat(tree2);
        nlens=al$len;
        alive=al$alive;
        tlens=(al$ted-al$tst)

        if(only_alive){
            mask= (alive>0);
        }else{

            if(omit_alive){
                ##alive=tree2trait_last(tree2,length(tree2[[1]]));
                ##print(alive);
                if(tedge<0){
                    mask= ((nlens>edge)+(alive>0))>0;
                }else{
                    mask= ((tlens>tedge)+(alive>0))>0;
                }
            }else{
                
                if(tedge<0){
                    mask=(nlens>edge);
                }else{

                    mask=(tlens>tedge);
                    ##print(tedge);
                    ##print(tlens);
                    ##print(as.numeric(mask));
                }
            }
        }
        

        
              
        tree2b=NULL;
        lis=which(mask);

        if(length(lis)>0){
            ##print(tree2);
            for(i in 1:length(lis)){
                tree2b = c(tree2b,list(tree2[[lis[i]]]));
            }
        }else{
            tree2b=tree2;
        }
        
    }else{
        tree2b=tree2;
    }
    
    return(tree2b);
    
}


cut_tree <- function(tree,tst=-1,ted=NULL){
    stat=tree2stat(tree);
    
    mask= (!((stat$ted<tst)+(stat$tst>ted)>0));
    mask_cut= (stat$tst<tst)+(stat$ted>ted)>0;
    
    treeb=NULL;
    for(i in 1:length(tree)){
        if(mask[i]){
            if(mask_cut[i]){
                phe0=tree[[i]];
                myt=phe0$t;
                lis=which((myt>=tst)*(myt<=ted)>0);
                
                mytlis=myt[lis];
                phe=geti(phe0,lis);
                mypid=phe$pid;
                
            if((myt[1]<tst)&&(mytlis[1]>tst)){
                phe$t[1]=tst;
            }

            
            if((mypid[length(mypid)]>0)&&(myt[length(myt)]>ted)&&(mytlis[length(mytlis)]<ted)){
                ##lis=c(lis,lis[length(lis)]);
                ##phe=geti(phe0,lis);
                phe=add_phenotype(phe,geti(phe,length(phe$pid)));
                phe$t[length(phe$t)]=ted;
                ##print(phe$t[length(phe$t)]);
            }

            
        }else{
            phe=tree[[i]];

        }

            if(phe$pid[1]>-1)treeb = c(treeb,list(phe));               

    }
    }


    ##stat=tree2stat(treeb);
    ##treeb1=select_tree(treeb,which(stat$ped!=stat$ped[i]));

    if(0){
        stat=tree2stat(treeb);
        pedu=sort(unique(stat$ped))
        treeb1=NULL;
        for(i in 1:length(pedu)){
            lis=which(stat$ped==pedu[i]);
            id=which.max(stat$len[lis]);
            mybr=treeb[[lis[id]]];
            treeb1 = c(treeb1,list(mybr));               
        }
    }

    if(1){
    stat=tree2stat(treeb);
    maskr= as.numeric(stat$ped-stat$pst!=0)
    maske=maskr*0;
    for(i in 1:length(maske)){
        if(length(which(stat$ped==stat$ped[i]))==1)maske[i]=1;
    }

    ##print(maske);
    ##print(maskr);
    treeb1=select_tree(treeb,which(maskr+maske>0));
    ##treeb1=select_tree(treeb,which(maskr>0));
    }

    
    return(treeb1);
}



div_time <- function(tree1){
    if(length(tree1)==1)cpid=tree1[[1]]$pid[1];
    if(length(tree1)>1){

        cpid=rep(0,length(tree1));
       
        for(i in 1:length(tree1)){
           
            pidi=tree1[[i]]$pid;
            for(j in 1:length(tree1)){
                if(i!=j){
                pidj=tree1[[j]]$pid;
                cpij=max(pidi[which(match(pidi,pidj)>0)]);
                if(cpij>cpid[i])cpid[i]=cpij;
             }
            }
        }
    }

    
   return(cpid);

}



find_common_ans <- function(tree2ba, p){
    stat=tree2stat(tree2ba);
    

    psts=stat$pst[stat$tst==min(stat$tst)];
        
    trp20=stree1(tree2ba,"pid");
    trp2=trp20;

    
    laspid=tree2lastpid(tree2ba);
    ##    cpi=c(unlist(lapply(trp20,function(x){x[1]})),unlist(lapply(trp20,function(x){tail(x,n=1)})))
        cpi=c(unlist(lapply(trp20,function(x){x[1]})),laspid)
    cpi=sort(cpi);
    
    for(i in 1:length(trp2)){
        lis=which(match(trp20[[i]],c(0,cpi))>0);
        trp2[[i]] = trp20[[i]][lis];
    }
        
    tp2=trp2;

    ##print(psts);
    
    for(i in 1:length(tp2)){
        
##        cat("i:",i," ");
        
        ##        while(sum(tp2[[i]]==0)==0){
        while(sum(psts==tp2[[i]][1])==0){
            for(j in 1:length(tp2)){
                if((i!=j)&&(tp2[[i]][1]>tp2[[j]][1])){
                
                li=which(match(tp2[[j]],tp2[[i]])>0);
                if(length(li)>0){
                    tj=tp2[[j]][1:(li[1]-1)];
                    tp2[[i]]=c(tj,tp2[[i]]);
                }
                    
                }
                
                
            }
        }
    }
    
    
 
    trp=tp2;
    n2=length(trp);
    ca=diag(n2)*0;
    ga=ca;
    for(i in 1:n2){
        
        for(j in 1:n2){
            if(ca[i,j]==0){
                if(i==j){
                    ca[i,j]=trp[[i]][length(trp[[i]])];
                    ga[i,j]=i;
                    
                }else{
                    lis=which(match(trp[[i]],trp[[j]])>0);
                    if(length(lis)>0){
                        ca[i,j]=trp[[i]][lis[length(lis)]];
                        ga[i,j]=j;
                    }
                }
            }
        }
    }


    cpid=rep(0,ncol(ca));

    if(length(cpid)>1){
        for(i in 1:length(cpid)){
            cpid[i]=max(ca[i,-i]);
            ##if(max(ca[i,-i])==-Inf){print(a$ninv);print(ca);print(tree2ba);}
        }
    }else{
        cpid=ca[1,1];
    }
        tdiv0=a$tree_phe$t[cpid+1];
        tdivy=a$tree_phe$y[cpid+1];
        tdivx=a$tree_phe$x[cpid+1];
    
        dismt= matrix(a$tree_phe$t[ca+1],ncol=ncol(ca));

  
    
    tdiv=tdiv0*0+a$timen-p$tedge;
    ##tdiv=ca$tdiv*0;

    
    for(i in 1:length(tdiv)){
        buf=dismt[i,-i];
        buft=a$timen-dismt[i,-i];
        li=which(buft>p$tedge);
        buf=buf[li];
        buft=buft[li];
        if(length(li)>0)tdiv[i]=buf[which.min(buft)];
    }
    
    ##    phe_ed=geti(a$tree_phe,(laspid+1));
    
    phe_ed=tree2phe_ed(tree2ba);
    
    return(list(ca=ca,cpid=cpid,tdiv=tdiv,tdiv0=tdiv0,tdivy=tdivy,tdivx=tdivx,dismt=dismt,treep=trp,treep0=trp2,phe_ed=phe_ed));
}


tree_sampling <- function(tree,st=NULL,ed=NULL,nsamp=10,tsamp=NULL,dp=0.1,tree_phe=a$tree_phe,n_all=a$n_all,pparam=a$pparam,out_analysis=FALSE,plot_each=TRUE){
    
    tu=sort(unique(tree_phe$t));
    ##n_all=a$n_all;
    ninvu=sort(unique(n_all$ninv));

    if(length(tsamp)>0){
        st=min(tsamp);
        ed=max(tsamp);
    }else{
        if(length(ed)==0)ed=tu[length(tu)];
        ##if(length(st)==0)st=0.5*(tu[1]+tu[2]);
        if(length(st)==0)st=tu[1];
        if(length(tsamp)==0)tsamp=seq(st,ed,,nsamp);
    }
    print(tsamp);
ph=NULL;
res=NULL;

if(plot_each)x11();
for(i in 1:length(tsamp)){
##i=1;
  
    cat("\n samp:",i,"   time:",tsamp[i]);
    ted=tsamp[i];

    treeb0=cut_tree(tree,tst=-1,ted=ted);
    
    treeb=NULL;
    if(length(treeb0)>1){
        treeb=select_tree(treeb0,which(tree2stat(treeb0)$alive>0));
        statb=tree2stat(treeb,timen=ted);
        caa=find_common_ans(treeb,pparam);
    }else{

        if(treeb0[[1]]$pid[length(treeb0[[1]]$pid)]>0)treeb=treeb0;
    }
    
    if(plot_each){
        dev.hold();
        plot_tree(tree,lty=2);
        if(length(treeb)>0)plot_tree(treeb,col="red",hold=1);
    

        lines(xran*c(-1,1),ted*c(1,1),col="black",lty=2);
        ##lines(xran*c(-1,1),tst*c(1,1),col="black",lty=2);
        lines(xran*c(-1,1),st*c(1,1),col="black",lty=2);
        
        if(length(treeb0)>1){
            points(caa$tdivx,caa$tdiv);
            points(caa$phe_ed$x,caa$phe_ed$t,pch=1,col="blue");
            
            for(kk in 1:length(caa$tdiv))lines(caa$phe_ed$x[kk]*c(1,1),c(caa$phe_ed$t[kk],caa$tdiv[kk]),col="black",lty=2);
        }
    
        dev.flush();
        Sys.sleep(dp);
    }
        
    if(length(treeb0)>1){
    phe=caa$phe_ed;
    phe$tdiv=caa$tdiv;
    phe$tsamp=rep(ted,length(phe$tdiv));
    }else{
        phe=tree2phe_ed(treeb0);
        phe$tdiv=0;
        phe$tsamp=ted;
        

    }
    tid = max(which(tu<=ted))-1;
    ##tid = min(which(tu>=ted));
    lis=which(n_all$ninv==tid)
    pid0=n_all$pid[lis];
    n0=n_all$n[lis];
    lism=match(phe$pid,pid0);
    print(paste(length((phe$pid)),length(unique(phe$pid)),length(pid0)));
    print(sort(unique(phe$pid)));
    print(sort(pid0));
    print(" ");
    
    if(sum(is.na(lism))){
        print("pid not found!!!");
        print(phe$pid);
        print(pid0);
    }
    
    phe$n=n0[lism];

    lismm=match(unique(phe$pid),phe$pid);
    phe=geti(phe,lismm);
    
    ph=c(ph,list(phe));

    if(out_analysis){
        dis=ted-caa$dismt;
        for(kk in 1:ncol(dis))dis[kk,kk]=0.0;

        lab=paste0(caa$phe_ed$gid,":",1:length(caa$phe_ed$gid));
        colnames(dis)=lab;
        rownames(dis)=lab;
        ##tree.phylo = as.phylo(hclust(as.dist(dis),method="average"));
        tree.phylo = as.phylo(hclust(as.dist(2*dis),method="average"));

       
        
        resb=list(tsamp=ted,nphe=length(phe[[1]]),nsum=sum(n0),lism=lism,phe=phe,n=n0,ca=caa,subtreea=treeb,subtree=treeb0,tree.phylo=tree.phylo,dis=dis);
        res=c(res,list(resb));
    }
}

print("");

    ph=data.frame(merge_phenotype(ph));
    if(out_analysis){
        return(list(phe=ph,res=res));
    }else{
        return(ph);
    }
    
}

tree2phylo <- function(tr2,lab=NULL,pparam=a$pparam,timen=a$timen,flag_check=TRUE){
    ted0=timen;
    ca=find_common_ans(tr2,pparam);
    tp2=ca$treep0;

    gid=tree2trait_ed(tr2,"gid");
    stat=tree2stat(tr2,timen=ted0);
    tst=stat$tst;
    ted=stat$ted;
    
    ##    lab=paste0(gid,"_",1:length(gid));
    if(length(lab)==0)lab=paste0(gid,":",1:length(gid));
    
    lab1=gsub(":","_",lab);
    
    pidt=a$tree_phe$t;
    ##nwk0=as.character(1:length(tp2));
    nwk0=lab1;
    
    nwk=nwk0;
    nwk1=nwk;
    ndt=rep(0.0,length(nwk0));
    
    ##li=which(stat$alive>0);

    ##ndt[li]=ted0-pidt[stat$ped[li]+1];

    bpid0=unlist(lapply(tr2,function(x){x$pid[2]}));
    bpid0[which(bpid0==1)]=0;
    bst0=unlist(lapply(tr2,function(x){x$t[2]}));
    bst0[which(bpid0==0)]=0.0;
    bed0=stat$ted;

    lis=order(-bpid0);
   
    
    bpid=bpid0[lis];
    bst=bst0[lis];
    bed=bed0[lis];
    bed1=bed;
    ##bx=x[lis];
    tr2b=select_tree(tr2,lis);
    tp2b=select_tree(tp2,lis);
    labb=lab1[lis];
    gidb=gid[lis];
    nwk0=labb;
    nwk=nwk0;
    for(k in 1:(length(bpid)-1)){
        myp=tp2b[[k]][1];
        cba0=unlist(lapply(tp2b,function(x){length(which(x[2:length(x)]==myp))}));
        cba0[k]=0;
        bid=which(cba0>0);
        
        if(flag_check)print(paste(k,bid));
        
        dt_k=bed1[k]-bst[k];
        dt_bid=bed1[bid]-bst[k];
        bed1[bid]=bst[k];
        
        nwk[bid]=paste0("(",nwk[bid],":",dt_bid,",",nwk[k],":",dt_k,")");
        nwk[k]="";
    }
    
    nwk1=nwk[nwk!=""];
    nwk=nwk1;
    
    if(flag_check)print(nwk);


    trn=ape::read.tree(text=paste0(nwk,";"));
    trnp=trn;

    trnp$tip.label=gsub("_",":",trnp$tip.label);
    ##trnp$tip.label= as.character(1:length(nwk));

    trnp$root.edge=bed1[length(bed1)];

    if(flag_check){

        times_tree=sort(unique(c(bst,bed)));
        times_trnp=c(0,sort(unique(as.numeric(nodeHeights(trnp))))+trnp$root.edge);
 
        print("comparison of times for tip and nodes")

    
        
        print("tree: tree");
        print(times_tree);
        
        print("phylo: trpn");
        print(times_trnp);
        
        tdiff=times_tree*0;
        for(i in 1:length(tdiff)){
            tdiff[i]=min(abs(times_trnp-times_tree[i]));
        }
        
        ##tdiff=times_tp2-times_tree[1:length(times_tree)];
        print(paste("Time difference max:",max(tdiff)))
        
        print(tdiff);
        print(paste("Time difference max:",max(tdiff)))


    }

    lis_order=match(trnp$tip.label,lab);
    stat1=stat;
    for(i in 1:length(stat1))stat1[[i]]=stat1[[i]][lis_order];

    
    return(list(tree.phylo=trnp,stat=stat1,stat0=stat,lab0=lab,ord=lis_order));
}


##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/ Basic functions for tree manipulation  _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

## convert a$tree into a tree containing only the element "ename" (e.g., "x", "y", "t", or "pid")
stree <- function(ename){
    lapply(a$tree,eval(parse(text=sprintf('function(x)x$%s',ename))));
}

## convert tree into a tree containing only the element "ename"
stree1 <- function(tree,ename){
    lapply(tree,eval(parse(text=sprintf('function(x)x$%s',ename))));
}
## convert tree into a tree containing only the element "ename"
stree1u <- function(tree,ename){
    unlist(stree1(tree,ename));
}

## convert tree into a tree containing only the element specified by
## "ename" that is the name or order of the element.
tree2trait <- function(tree,trait){
   
    if(is.numeric(trait)){
        if(trait==0){
            (lapply(tree, function(phe){phe$t}))
        }else{
            (lapply(tree, function(phe){phe[[trait]]}))
        }
    }else{
         stree1(tree,trait);
    }

}

## unlist version of tree2trait
tree2traitu <- function(tree,trait){
  ##  print(is.character(trait));
    
    if(is.numeric(trait)){
        if(trait==0){
            unlist(lapply(tree, function(phe){phe$t}))
        }else{
            unlist(lapply(tree, function(phe){phe[[trait]]}))
        }
    }else{
        unlist(stree1(tree,trait));
    }
}


## extract information from each branche of tree, including whether each branch
## is still alive or not 
tree2alive <- function(tree){
    treep=stree1(tree,"pid");
    treet=stree1(tree,"t");
    pid=unlist(lapply(treep,getlastpid));
    len=unlist(lapply(treep,length));
    alive=unlist(lapply(treep,function(x){tail(x,n=1)}));
    t=unlist(lapply(treet,function(x){tail(x,n=1)}));
        return(list(pid=pid,len=len,alive=alive,t=t));
}

## extended one from tree2alive

tree2stat <- function(tree,timen=a$timen){
    treep=stree1(tree,"pid");
    treet=stree1(tree,"t");
    pst=unlist(lapply(treep,function(x){x[1]}));
    ped=unlist(lapply(treep,getlastpid));
    len=unlist(lapply(treep,length));
    alive=unlist(lapply(treep,function(x){tail(x,n=1)}));

    tst=unlist(lapply(treet,function(x){x[1]}));
    ted=unlist(lapply(treet,function(x){tail(x,n=1)}));
    ted[alive>0]=ted[alive>0]*0+timen;
    tlen=ted-tst;
        return(list(alive=alive,len=len,pst=pst,ped=ped,tst=tst,ted=ted,tlen=tlen));
}


## get a part of phe (resident phenotypes) specified by st and ed
getii <- function(phe,st=1,ed=length(phe[[1]])){
    phe1=phe;
    for(j in 1:length(phe1))phe1[[j]]=(phe[[j]])[st:ed];
    return(phe1);
    ##lapply(z,function(z1)(z1[i]));
}

## get the first and last element from each branch
tree2trait_range <- function(tree,trait){
    
    if(is.numeric(trait)){
        
        if(trait==0){
            v=(lapply(tree, function(phe){phe$t}))
        }else{
            v=(lapply(tree, function(phe){phe[[trait]]}))
        }
        name=names(tree[[1]])[trait];
    }else{
        
        v=stree1(tree,trait);
        name=trait;
        
        
    }
    
    
    vst=unlist(lapply(v,function(x){x[1]}));
    ved=unlist(lapply(v,function(x){tail(x,n=1)}));
    return(list(st=vst,ed=ved,name=name));
}

## get the first element of the trait from each branch
tree2trait_st <- function(tree,trait){
   tree2trait_range(tree,trait)$st;
}
## get the last element of the trait from each branch
tree2trait_ed <- function(tree,trait){
   tree2trait_range(tree,trait)$ed;
}

## get the last pid from an array of pid
getlastpid <- function(pid){
    nn=length(pid);
       if(pid[nn]==-1){
           pidl=pid[nn-1];
       }else{
           pidl=pid[nn];
           
       }
    return(pidl);
}

## extract last pid from each branch
tree2lastpid <- function(tree,len=FALSE){
    treep=stree1(tree,"pid");
    pid=unlist(lapply(treep,getlastpid));
    if(len==FALSE){
        return(pid);
    }else{
        len=unlist(lapply(treep,length));
        alive=unlist(lapply(treep,function(x){tail(x,n=1)}));
        return(list(pid=pid,len=len,alive=alive));
    }
}


## extract first phenotype from each branch
tree2phe_st <-function(tree){
    if(length(tree)>1){
        tree1=lapply(tree, function(phe){geti(phe,1)});
        merge_phenotype(tree1);
    }else{
        geti(tree[[1]],1);
    } 
}

## extract last phenotype from each branch
tree2phe_ed <-function(tree){
    if(length(tree)>1){
        tree1=lapply(tree, function(phe){geti(phe,length(phe[[1]]))});
        merge_phenotype(tree1);
    }else{
        geti(tree[[1]],length(tree[[1]]$pid));
    }
}

## merge branches of a tree into a single branch
merge_phenotype <- function(tree){
    phe0=tree[[1]];
    for(j in 2:length(tree)){
        phe0=add_phenotype(phe0,tree[[j]]);
    }
    return(phe0);
}

## get i-th phenotypes from phe
#' @export
geti <- function(phe,i){
    phe1=phe;
    for(j in 1:length(phe1))phe1[[j]]=(phe[[j]])[i];
    return(phe1);
    ##lapply(z,function(z1)(z1[i]));
}

select_tree <- function(tree,lis){
        treeb=NULL;
        for(i in 1:length(lis)){
            treeb = c(treeb,list(tree[[lis[i]]]));
        }
        return(treeb);


}

## add a phenotype to phe
#' @export
add_phenotype <- function(phe,phe1){
    for(j in 1:length(phe)){
        phe[[j]]=c(phe[[j]],phe1[[j]]);
    }
    return(phe);
}

## calculate velocity for directionary evolution of trait
tree2velocity <- function(tree,trait,merg=TRUE){
    ##trait="x";
    treex=tree2trait(tree,trait);
    ##trait="y";
    ##treey=tree2trait(tree2,trait);
    
    treep=tree2trait(tree,"pid");
    treet=tree2trait(tree,"t");
    
    tagid=c(1,2);
    
    treev=NULL;
    for(i in 1:length(treex)){
        tt=treet[[i]];
        xx=treex[[i]];
        ##  yy=treey[[i]];
        pp=treep[[i]];
        if(pp[length(pp)]==-1)pp[length(pp)]=pp[length(pp)-1];
        lis= c(which(diff(tt)>0) ,length(tt));
        xx=xx[lis];
        ##yy=yy[lis];
        tt=tt[lis];
        pp=pp[lis];
        
        if(length(tt)>1){
            xdif=diff(xx);
            tdif=diff(tt);
            vx=xdif/tdif;
            ##vx=(diff(xx)/diff(tt));
            ##  vy=(diff(yy)/diff(tt));
            xx=xx[1:length(vx)];
            ##yy=yy[1:length(vx)];
            tt=tt[1:length(vx)];
            pp=pp[1:length(vx)];
            
            treev=c(treev,list(list(vel=vx,trait=xx,t=tt,pid=pp,dif=xdif,tdif=tdif)));
        }
    }
    if(merg)treev=merge_phenotype(treev);
    return(treev);   
}



##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/      Functions for simulation      _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##






#' @export
fitness_all <- function(phe,en){
    buf=phe[[1]]*0.0;
    for(i in 1:length(phe[[1]]))buf[i]=.simevol_func$fitness(geti(phe,i),phe,en);
    return(buf);
}


max_fitness <- function(phe,en){
    buf=0.0;
    for(i in 1:length(phe[[1]])){
        buf1=abs(.simevol_func$fitness(geti(phe,i),phe,en));
        if(buf<buf1)buf=buf1;
    }
    return(buf);
}


#' @export
pop_dynamics0 <- function(t,n,parms){    
    dn=n*0;
    for(i in 1:length(dn)){
          dn[i]=n[i]*.simevol_func$fitness(geti(parms$phe,i),parms$phe,n);
    }
    list(dn);
}

#' @export
mutate0 <- function(phe){
    ##mut=phe;
    mut=phe[1:(length(phe)-2)];
    for(i in 1:(length(phe)-2))mut[[i]]=phe[[i]]+rnorm(1,mean=0.0,sd=a$sparam$m_sd);
    return(mut);
}


#' @export
invader <- function(phe=a$phe,en=a$en,mutate=.simevol_func$mutate,timen=0.0,amp_invf=a$sparam$amp_invf){
    n=en[1:length(phe[[1]])];
    sum_n=sum(n);
    flag_invade=0;
    x_mut=0.0;
    nsp=length(n);
    while(flag_invade==0){
        timen = timen-(1.0/a$sparam$m_rate)*amp_invf*(1.0/sum_n)*log(runif(1)+1e-20);
        buf=0.0;
        m_spe=0;
        rv=runif(1);
        while((rv >= buf)&&(m_spe<nsp)){
            m_spe=m_spe+1;
            buf=buf+n[m_spe]/sum_n;
        }
        phe_par=geti(phe,m_spe);
        mutant=mutate(phe_par);
        
        ##fb=fitness(mutant,phe,en)-fitness(phe_par,phe,en);
        fb=fitness(mutant,phe,en);
        
        if(fb*amp_invf>1.0){
            cat("Invasion fitness larger than 1 !! ",fb,fb*amp_invf,"\n");
            a$sparam$fit_over <<- c(a$sparam$fit_over,fb);
        }
        if( fb*amp_invf> runif(1)){
            flag_invade=1;
        }
    }
    mutant=c(mutant,list(t=timen,pid=a$ninv+1));
    return(list(phe=mutant,n=a$sparam$n_mutant_init,f=fb,timen=timen,pid_par=phe_par$pid));
}



#' @export
add_invader <- function(phe,en,inv){    
    edim=length(en)-length(phe[[1]]);
    nspe=length(phe[[1]]);
    ebuf=c(NULL);
    if(edim>0)ebuf=en[(nspe+1):length(en)];
    n=en[1:nspe];

    phe=add_phenotype(phe,inv$phe);        
    n=c(n,inv$n);
    en=c(n,ebuf);
    nspe=length(n);

    return(list(phe=phe,en=en,n=n,nspe=nspe,inv=inv));
}



#' @export
remove_extinct <- function(en,phe,edge_extinct=NULL){
    die=c(NULL);
    nspe=length(phe[[1]]);
    
    for(i in 1:nspe){
        if(length(which(en[i]<edge_extinct))>0){
            die=c(die,i);
        }
    }

    
    if(length(die)>0){
                
        en=en[-die];
        for(i in 1:length(phe)){
            val=phe[[i]];
            val=val[-die];
            phe[[i]]=val;
        }
        
    }

    return(list(phe=phe,en=en,die=die));
}


#' @export
rootfunc <- function(t,n,parms){
    if(parms$edim==0){
        return(n-0.1*parms$edge_extinct);
    }
    else{
        n=n-0.1*parms$edge_extinct;
        n[(parms$nspe+1):length(n)]=1.0;
        return(n);
    }
}

#' @export
eventfunc <- function(t,n,parms){
    lis=which(n<parms$edge_extinct);
     if(length(lis)>0)n[lis]=n[lis]*0.0;
     
    return(n);
}


#' @export
simpop_check <- function(phe=a$phe,en=a$en,inv=invader(),pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,check=TRUE,nrad=7,drad=1.0,divrad=64,xlim=c(NULL),out=FALSE,reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
    res=simpop_invade(phe=phe,en=en,inv=inv,pop_dynamics=pop_dynamics,set_parms=set_parms,edge_extinct=edge_extinct,edge_fit=edge_fit,check=check,nrad=nrad,drad=drad,divrad=divrad,xlim=xlim,reset_win=reset_win,logt=logt,param_desolve=param_desolve);
    if(out==TRUE)return(res);
}
    
#' @export
simpop_invade <- function(phe=a$phe,en=a$en,inv=a$inv,pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,check=FALSE,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad,xlim=c(NULL),reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
##.ee.append("simpop_invade",environment())

    state=add_invader(phe,en,inv);       
    state_new=simpop(state$phe,state$en,pop_dynamics=pop_dynamics,set_parms=set_parms,edge_extinct=edge_extinct,edge_fit=edge_fit,check=check,nrad=nrad,drad=drad,divrad=divrad,xlim=xlim,reset_win=reset_win,logt=logt,param_desolve=param_desolve);
    return(state_new);
    
}


#' @export
simpop_set_parms <- function(phe1,en1,set_parms=NULL){
            if(class(set_parms)=="function"){
                parm=set_parms(phe1,en1);
            }else{
                parm=list(phe=phe1);               
            }
return(parm);
}

#' @export
get_last_en <- function(en_next){
    lent=length(en_next[,1]);
    en=as.numeric(en_next[lent,]);
    return(en[2:length(en)]);    
}


#' @export
simpop <- function(phe1=a$phe,en1=a$en,pop_dynamics=.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extince,edge_fit=a$sparam$edge_fit,check=FALSE,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad,xlim=c(NULL),reset_win=FALSE,logt=TRUE,param_desolve=a$sparam$param_desolve){
    hini=param_desolve$hini;
    hmax=param_desolve$hmax;
    rtol=param_desolve$rtol;
    atol=param_desolve$atol;
    method=param_desolve$method;
    
##.ee.append("simpop",environment())
    nspe0=length(phe1[[1]]);
    nspe1=nspe0;
    edim1=length(en1)-nspe1;
    if(check==FALSE){
        for(irad in 1:nrad){            
            parm=c(simpop_set_parms(phe1,en1,set_parms),list(edge_extinct=edge_extinct,nspe=nspe1,edim=edim1));
             mytime=10^seq(drad*irad,drad*(irad+1),,divrad);
            n_next=deSolve::ode(y=en1,times=mytime,func=pop_dynamics,parms=parm,method=method,hini=hini,hmax=hmax,rootfun=rootfunc,rtol=rtol,atol=atol);
            
            res=remove_extinct(get_last_en(n_next),phe1,edge_extinct=edge_extinct);
            phe1=res$phe;
            en1=res$en;
            nspe1=length(phe1[[1]]);
            edim1=length(en1)-nspe1;
            max_fit=max(abs(fitness_all(phe1,en1)));
            max_t=max(n_next[,1]);
            ##if((max_t<times[irad+1])&&(length(phe1[[1]])==nspe0)){
            ##    cat("population dynamis failed!! time:", max_t,"max_fit:",max_fit,"\n");
            ##    a$sparam$flag_halt<<-T;
            ##}

            if(length(which(res$die==nspe0))>0){
                cat("inveder extinct!!\n");
                ##if(a$sparam$error_stop==TRUE);                
                ##a$sparam$flag_halt<<-T;
                }
                        

            if(max_fit<edge_fit)break;
            ##hini=min(10^(drad*irad),1e-3/max_fit);
            hini=min(10^(drad*irad-1),1e-4/max_fit); ## this could be improved
            ##cat("hini:",hini,"\n")
        }
        
        return(list(phe=phe1,en=en1,irad=irad));
        
    }
    
    if(check==TRUE){
        phe0=phe1;
        en0=en1;
        
        parm=c(simpop_set_parms(phe1,en1,set_parms),list(edge_extinct=edge_extinct,nspe=nspe1,edim=edim1));        
        for(irad in 1:nrad){    
            mytime=10^seq(drad*irad,drad*(irad+1),,divrad);
            n_next=deSolve::ode(y=en1,times=mytime,func=pop_dynamics,parms=parm,method=method,hini=hini,hmax=hmax,rootfun=rootfunc,atol=atol,rtol=rtol);
            en1=(n_next[nrow(n_next),])[2:ncol(n_next)];
            if(irad==1)n_nex=n_next;
            if(irad>1)n_nex=rbind(n_nex,n_next);
            max_fit=max(abs(fitness_all(phe1,en1)));
            if(max_fit<edge_fit)break;
            hini=min(10^(drad*irad-1),1e-4/max_fit);
        }
        n_next=n_nex;
      
        res=remove_extinct(get_last_en(n_next),phe1,edge_extinct=edge_extinct);
        phe1=res$phe;
        en1=res$en;
        if(length(which(res$die==nspe0))>0){
                cat("inveder extinct!!\n");
        }
        plot_simpop(n_next,phe0,edge_extinct,reset_win=reset_win,logt=logt,xlim=xlim);
        print("fitness");
        print(fitness_all(phe1,en1));
        
        return(list(phe=phe1,en=en1,n_next=n_next));
        
        
    }
}

#' @export
plot_simpop <- function(n_next,phe0,edge_extinct,reset_win=FALSE,logt=FALSE,xlim=NULL){
    nspe=length(phe0[[1]]);
    edim=ncol(n_next)-nspe-1;
    
            nam=colnames(n_next);
            name_n=nam[2:(nspe+1)];
            colnames(n_next)[2:(nspe+1)]=paste0("n",name_n);
            
            en_last=n_next[nrow(n_next),];
            n_last=en_last[2:(nspe+1)];
            t_last=en_last[1];
            if(edim>0){
                name_e=nam[(nspe+2):length(nam)];
                colnames(n_next)[(nspe+2):length(nam)]=paste0("e",seq(1,length(name_e)));
                e_last=en_last[(nspe+2):length(en_last)];
            }
            lis=which(n_last<edge_extinct);
            color=rep("blue",ncol(n_next)-1);
            color[nspe]="green";
            if(length(lis)>0)color[lis]="red";
            
            npanel=nspe+edim;
            ncols=3;
            nrows=as.integer(npanel/ncols)+as.integer(npanel%%ncols>0);
            
        if((reset_win==FALSE)&&(a$winid[3]>0)){
            dev.set(a$winid[3]);
            plot.new();
        }
        else{
            X11(width=ncols*3.0,height=nrows*1.5);
            a$winid[3]<<-cur.dev();
        }
        om_left=1;
        om_right=0;
        om_bottom=1;
        om_top=0;

        par(oma = c(om_bottom, om_left, om_top, om_right),mar=c(3,3,2,3),mgp=c(2,0.7,0)); ##bottom,left,top,right
       
        par(mfrow=c(nrows,ncols));
        name=colnames(n_next);
        for(i in 1:(ncol(n_next)-1)){
            if(i<=nspe){
                ylim=c(0,max(n_next[,(i+1)])*1.2);
                ##ylim=c(NULL);
                tit="";
               
                for(j in 1:(length(phe0)-2))tit=paste(tit,sprintf("%s:%f",names(phe0)[j],phe0[[j]][i]));
                if(logt){plot(n_next[,1],n_next[,(i+1)],log="x",col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,xlim=xlim,main=tit);}
                else{plot(n_next[,1],n_next[,(i+1)],col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,xlim=xlim,main=tit);}
            }
            else{
                if(logt){plot(n_next[,1],n_next[,(i+1)],log="x",col="orange",type="l",xlab=name[i+1],ylab="value",xlim=xlim);}
                else{plot(n_next[,1],n_next[,(i+1)],col=color[i],type="l",xlab=name[i+1],ylab="Density",ylim=ylim,main=tit,xlim=xlim);}
            }
        }
        print(phe0);
        print(t_last);
        print(n_last);
        if(edim>0)print(e_last);
}


##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/        Plotting functions          _/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

#' @export
calc_fitness_land <- function(phe,en,xid=1,yid=2,xmin=-1.0,xmax=1.0,ymin=-1.0,ymax=1.0,ndiv=64){
    n=en[1:length(phe[[1]])];
   
    xx=seq(xmin,xmax,,ndiv);
    yy=seq(ymin,ymax,,ndiv); 
    
    land=matrix(ndiv*ndiv,nrow=ndiv,ncol=ndiv)*0.0;

    if(length(a$pparam$fitness_contour_phe)>0){
        phe1=a$pparam$fitness_contour_phe;
    }
    else{
        phe1=phe;
        pdim=length(phe)-2;
        if(pdim>1){
            for(k in 1:pdim)phe1[[k]]=sum(phe1[[k]]*n)/sum(n);
        }
    }
    for(i in 1:ndiv){
        for(j in 1:ndiv){
            phe1[[xid]]=xx[i];
            phe1[[yid]]=yy[j];
            land[i,j]=.simevol_func$fitness(phe1,phe,en);
        }
    }
    return(list(x=xx,y=yy,fit=land));
}

#' @export
plot_fitness_contour <- function(land,p){
    contour(land$x,land$y,land$fit,add=1,levels=p$fcon$levels,col=p$fcon$col,lty=p$fcon$lty);
}


#' @export
plot_init <- function(p,reset_win=TRUE){
    if(reset_win){
    X11(width=p$winw,height=p$winh,title=paste("runid:",a$runid,"  ",a$runname,"    dev:",(cur.dev()+1)));
    a$winid[1] <<-cur.dev();
    }
    else{
        dev.set(a$winid[1]);
    }
    par(pty=a$pparam$win_style);
        
    npanel=length(a$phe)-2+a$edim;
    ncols=npanel;
    ##nrows=as.integer(npanel/ncols)+as.integer(npanel%%ncols>0);
    nrows=1;
    om_left=1;
    om_right=0;
    om_bottom=1;
    om_top=0;

   
    if(a$show_subwin==TRUE){
        if(reset_win){
        ##X11(width=ncols*3.0,height=nrows*4);
            X11(width=p$winw_sub,height=nrows*p$winh_sub,title=paste("runid:",a$runid,"  ",a$runname,"    dev:",(cur.dev()+1)));

            a$winid[2]<<-cur.dev();
        }
        else{
            dev.set(a$winid[2]);        
        }
        
        par(oma = c(om_bottom, om_left, om_top, om_right),mar=c(3,3,2,3),mgp=c(2,0.7,0)); ##bottom,left,top,right
        ##par(family=font_family) ;
        par(mfrow=c(nrows,ncols));
    }
        ##return(c(win0,win1));
    }

#' @export
pparam0=list(
    fitness_contour=TRUE,
    fitness_contour_phe=NULL,
    time_lab="Time",
    cex.lab=1.1,
    xid=1,
    yid=2,
    xlim=c(NULL),
    ylim=c(NULL),
    win_style=NULL,
    fcon=list(levels=c(0.0,0.1),col=c("red","gray"),lwd=1,lty=c(1,1)),
    resi=list(col="black",bg="green",cex=0.4,pch=21,amp=1.0,ampn=10,bg2="red",bgid=-1),
    traj=list(col="blue",cex=1.0,pch=17,every=10,lwd=0.5),
    env=list(col="orange",lwd=1),
    plot_mask=NULL,
    trait_names=NULL,
    nv_names=NULL,
    fitness_contour=TRUE,
    fitness_contour_phe=NULL,
    pal=(rev(rainbow(100,end=0.7))),
    palfunc=NULL,
    winh=5,
    winw=5,
    winh_sub=4,
    winw_sub=7,
    omit_plot=FALSE
);

#' @export
adjust_n <- function(n,amp,ampn,offset){
    return (offset+amp*(log(n*ampn+1)/log(ampn+1)));
}

#' @export
plot_lim <- function(xid,yid,p,main=NULL){
    if((xid>0) && (yid>0)){
        xp=c(min(a$tree_phe[[xid]]),max(a$tree_phe[[xid]]));
        yp=c(min(a$tree_phe[[yid]]),max(a$tree_phe[[yid]]));
        xlab=p$trait_names[xid];
        ylab=p$trait_names[yid];
    }
    else{
        if(xid>0){
            xp=c(min(a$tree_phe[[xid]]),max(a$tree_phe[[xid]]));
            yp=c(min(a$tree_phe$t),max(a$tree_phe$t));
            xlab=p$trait_names[xid];
            ylab=p$time_lab;
            
        }
        if(yid>0){
            xp=c(min(a$tree_phe$t),max(a$tree_phe$t));
            yp=c(min(a$tree_phe[[yid]]),max(a$tree_phe[[yid]]));
            xlab=p$time_lab;
            ylab=p$trait_names[yid];
            
        }
        
    }
    

        plot(xp,yp,type="n",,cex.lab=p$cex.lab,xlim=p$xlim,ylim=p$ylim,xlab=xlab,ylab=ylab,fg=p$fgcol,col.axis=p$fgcol,col.lab=p$fgcol,main=main);
}

#' @export
plot_lim_sub <- function(xid,p){
    xp=c(min(a$tree_phe[[xid]]),max(a$tree_phe[[xid]]));
    yp=c(min(a$tree_phe$t),max(a$tree_phe$t));
    xlab=p$trait_names[xid];
    ylab=p$time_lab;
    
    plot(xp,yp,type="n",,cex.lab=p$cex.lab,xlim=p$sub_xlim[[xid]],ylim=p$sub_ylim[[xid]],xlab=xlab,ylab=ylab,fg=p$fgcol,col.axis=p$fgcol,col.lab=p$fgcol);
}


adj_tree <- function(tree){
    for(i in 1:length(tree)){
        b=tree[[i]];
        lis=which(a$phe$pid==b$pid[length(b$pid)]);
        if(length(lis)>0){
            b=geti(b,length(b$pid));
            b$t=a$timen;
            tree[[i]]=add_phenotype(tree[[i]],b);
        }
    }
    return(tree);
}

        
pline <- function(b,xid,yid,...){
        lines(b[[xid]],b[[yid]],...);
}



pline1 <- function(b,xid,yid,...){
    if(length(b$col)==1){   
        lines(b[[xid]],b[[yid]],col=b$col,...);
    }else{
        lis=unique(b$gid);
        li=which(b$gid==lis[1]);
        lines(b[[xid]][li],b[[yid]][li],col=b$col[1],...);
        
        for(i in 2:length(lis)){
            li=which(b$gid==lis[i]);
            li=c((li[1]-1),li);
            lines(b[[xid]][li],b[[yid]][li],col=b$col[i],...);
        }
    }
}

plot_traj_line1 <-function(tree,xid,yid,p,mask=NULL){
  
    if(xid==0)xid=a$pdim+1;
    if(yid==0)yid=a$pdim+1;
    ##print(p$traj_line_color);
    if(p$traj_line_color){
        lapply(tree,pline1,xid,yid,pch=p$traj$pch,cex=p$traj$cex,lwd=p$traj$lwd);
    }else{
        lapply(tree,pline,xid,yid,col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex,lwd=p$traj$lwd);    
        
        
    }
   
    

}

plot_tree <-function(tree,xid=1,yid=0,col="blue",xlim=NULL, ylim=NULL,cex=1,lwd=1,lty=1,p=NULL,hold=FALSE){
    
    if(xid==0)xid=a$pdim+1;
    if(yid==0)yid=a$pdim+1;

    if(length(p)>0){
        col=p$traj$col;
        cex=p$traj$cex;
        lwd=p$traj$lwd
    }


    
    if(hold==0){
        x=tree2traitu(tree,xid);
        y=tree2traitu(tree,yid);
        
        ##plot(x,y,type="n",xlim=c(min(x),max(x)),ylim=c(min(y),max(y)));
        if(length(p)>0){
            plot(x,y,type="n",xlim=xlim,ylim=ylim,fg=p$fgcol,col.axis=p$fgcol,col.lab=p$fgcol);
            }else{
        plot(x,y,type="n",xlim=xlim,ylim=ylim);
}
##              plot(x,y,type="n",xlim=xlim,ylim=ylim,fg=a$pparam$fgcol,col.axis=a$pparam$fgcol);
    }

    if(class(col)=="function"){
        cols = col(tree);
        for(kk in 1:length(tree))tree[[kk]]$col=cols[[kk]];
        
        v=lapply(tree,pline1,xid,yid,cex=cex,lwd=lwd,lty=lty);    

    }else{
        v=lapply(tree,pline,xid,yid,col=col,cex=cex,lwd=lwd,lty=lty);    
        
        
    }
   
    

}



#' @export
plot_traj_line <-function(tree,xid,yid,p,mask=NULL){
    if(xid==0)xid=a$pdim+1;
    if(yid==0)yid=a$pdim+1;
        lapply(tree,pline,xid,yid,col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex,lwd=p$traj$lwd);    
}


#' @export
plot_traj <-function(xid,yid,p,mask=NULL){
    if(length(mask)==0){
        if((xid>0)&&(yid>0)){
            points(a$traj$phe[[xid]],a$traj$phe[[yid]],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
        else{        
            if(xid>0)points(a$traj$phe[[xid]],a$traj$t,col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
            if(yid>0)points(a$traj$t,a$traj$phe[[yid]],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
    }
    else{
        if((xid>0)&&(yid>0)){
            points((a$traj$phe[[xid]])[mask],(a$traj$phe[[yid]])[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
        else{        
            if(xid>0)points((a$traj$phe[[xid]])[mask],(a$traj$t)[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
            if(yid>0)points((a$traj$t)[mask],(a$traj$phe[[yid]])[mask],col=p$traj$col,pch=p$traj$pch,cex=p$traj$cex);
        }
    }
}

#' @export
plot_phe <-function(xid,yid,p){
    timen=a$timen;
    nspe=length(a$phe[[1]]);
    if((xid>0)&&(yid>0)){
        points(a$phe[[xid]],a$phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    else{
        if(xid>0)points(a$phe[[xid]],rep(timen,nspe),col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
        if(yid>0)points(rep(timen,nspe),a$phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(a$n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    
}

plot_phe1 <-function(xid,yid,p,phe=a$phe,n=a$n,timen=a$timen){
    nspe=length(phe[[1]]);
    if((xid>0)&&(yid>0)){
        points(phe[[xid]],phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    else{
        if(xid>0)points(phe[[xid]],rep(timen,nspe),col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
        if(yid>0)points(rep(timen,nspe),phe[[yid]],col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
    }
    
}

#' @export
plot_1dim <- function(p){
    phe=a$phe;
    en=a$en;
    n=a$n;
    edim=a$edim;
    nspe=length(phe[[1]]);
    
    ndiv=128;
    ranx=1.2*max(abs(a$tree_phe[[1]]));
    xx1=seq(-ranx,ranx,,ndiv);
    land=xx1*0.0;
    fit1=phe[[1]]*0.0;
    x1=phe[[1]];
    phe1=geti(phe,1);
        for(i in 1:length(xx1)){
            phe1[[1]]=xx1[i];
            land[i]=fitness(phe1,phe,en);
        }
        for(i in 1:length(fit1)){
            phe1[[1]]=x1[i];
            fit1[i]=fitness(phe1,phe,en);
        }
        plot(xx1,land,col=p$fcon$col[1],lty=p$fcon$lty[1],type="l",xlab=p$trait_names[1],ylab="Fitness",ylim=c(-max(land)*0.2,max(land)*1.0),cex.lab=p$cex.lab);
        lines(xx1,land*0,col=p$fcon$col[2],lty=p$fcon$lty[2]);
        
        points(x1,fit1,col=p$resi$col,bg=p$resi$bg,pch=p$resi$pch,cex=adjust_n(n,p$resi$amp,p$resi$ampn,p$resi$cex));
}


#' @export
plot_func0 <- function(traj_line=TRUE){
    phe=a$phe;
    en=a$en;
    n=a$n;
    edim=a$edim;
    traj=a$traj;
    p=a$pparam;
    nspe=length(phe[[1]]);
    timen=a$timen;
    
    dev.set(a$winid[1]);
    
    xid=p$xid;yid=p$yid;
    trait_names=p$trait_names;
    env_names=p$env_names;
    plot_mask=p$plot_mask;

    dev.hold();
    if(class(.simevol_func$palfunc)=="function"){
            npal=length(p$pal);
            bgval=.simevol_func$palfunc(phe,en);
        
            ##cmax=max(bgval)+1e-10;
            ##cmin=min(bgval)-1e-10;
            ##cid=as.integer((npal-1)*(bgval-cmin)/(cmax-cmin))+1;
            cid=as.integer((npal-1)*bgval)+1;
            p$resi$bg=p$pal[cid];
    }else{
        if(p$resi$bgid>-1){
            bgid=p$resi$bgid;
            npal=length(p$pal);
            if(bgid==0){
                bgval=n;
            }else{
                bgval=phe[[bgid]];
            }
            cmax=max(bgval)+1e-10;
            cmin=min(bgval)-1e-10;
            cid=as.integer((npal-1)*(bgval-cmin)/(cmax-cmin))+1;
            p$resi$bg=p$pal[cid];
        }
    }
        
    
    if(length(trait_names)==0)trait_names=names(phe);
    if((length(env_names)==0)&&(edim>0))env_names=paste0("Env",seq(edim));
        
    if((length(phe)-2)==1){
        plot_1dim(p);
    }else{        
       plot_lim(xid,yid,p);
       if((xid>0)&&(yid>0)&&p$fitness_contour){
           ranx=max(1.2*max(abs(traj$phe[[xid]])),0.001);        
           rany=max(1.2*max(abs(traj$phe[[yid]])),0.001);        
           
           land=calc_fitness_land(phe,en,xid,yid,xmin=-ranx,xmax=ranx,ymin=-rany,ymax=rany);
       }            
       if(traj_line==TRUE){
           tree1=adj_tree(a$tree);
           plot_traj_line(tree1,xid,yid,p);
       }
       else{
           plot_traj(xid,yid,p);
       }
        if((xid>0)&&(yid>0)&&p$fitness_contour)plot_fitness_contour(land,p);
        plot_phe(xid,yid,p);
    }

        
    dev.flush();
    
    if(a$show_subwin==TRUE){
        dev.set(a$winid[2]);
        dev.hold();
    for(i in 1:(length(phe)-2)){
        ##plot_lim(i,0,p);
        plot_lim_sub(i,p);
        if(traj_line==TRUE){
            if((length(phe)-2)==1)tree1=adj_tree(a$tree);
                       plot_traj_line(tree1,i,0,p);
       }
       else{
                       plot_traj(i,0,p);
       }


            plot_phe(i,0,p);
        
    }
    if(edim>0){
        for(i in 1:edim){
            plot(traj$e[,i],traj$te,col=p$env$col,lwd=p$env$lwd,type="l",xlab=env_names[i],ylab="Time",cex.lab=p$cex.lab);  
        }
    }
        dev.flush();
    }

    
}

#' @export
ptraj <- function(li=1){
    if(li==1)plot_func0(traj_line=TRUE);
    if(li==0)plot_func0(traj_line=FALSE);

}

##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/      Default functions for data output     _/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##


#' @export
output0 <- function(timen,z,n){
    options(scipen=100); 
    nspe=length(z[[1]]);
    cat(timen,nspe,"\n",file=a$file_data,append=TRUE);
    cat(a$phe$pid+1,"\n",file=a$file_data,append=TRUE);
    for(i in 1:nspe)cat(z$x[i],z$y[i],n[i],"\n",file=a$file_data,append=TRUE);
    cat("\n",file=a$file_data,append=TRUE);
    options(scipen=0);
}


#' @export
output_tree0 <- function(fname){
    options(scipen=100);
    q=unlist(lapply(a$tree,function(z)(length(z$pid))));
    cat(length(a$tree),max(q),"\n\n",file=fname,append=FALSE);
    for(i in 1:length(a$tree)){
        time_ed=-1.0;
        blen=length(a$tree[[i]]$pid);
        if(a$tree[[i]]$pid[blen]==-1)time_ed=a$tree[[i]]$t[blen];
        time_st=a$tree[[i]]$t[1];
        
        cat(i,blen,time_st,time_ed,"\n",file=fname,append=TRUE);
        write((a$tree[[i]]$pid+1),ncolumns=50,file=fname,append=TRUE);
        cat("\n",file=fname,append=TRUE);
    }

    cat("\n\n",file=fname,append=TRUE);

    cat(a$pdim+2,length(a$tree_phe$pid),"\n",file=fname,append=TRUE);

    for(i in 1:a$pdim){
        buf=unlist(a$tree_phe[[i]]);
        write(buf,file=fname,append=TRUE,ncolumns=50);
       cat("\n",file=fname,append=TRUE);    
    }

    
    write(a$tree_phe$t,file=fname,append=TRUE,ncolumns=50);
    cat("\n",file=fname,append=TRUE);    
    write((a$tree_phe$pid+1),file=fname,append=TRUE,ncolumns=50);
    cat("\n",file=fname,append=TRUE);
    options(scipen=0);
}



##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##
##_/_/_/_/_/_/     Functions for simulation controll      _/_/_/_/_/_/##
##_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/##

#' @export
comwin <- function(make_win=0,mydir=a$sparam$mydir){
    ##cat("library(simevol);\n",file=".Rprofile");
    ##cat(sprintf('source("simevol1g.R");file_command<<-"%scommand";\n',mydir),file=".Rprofile");
    cat(sprintf('library(simevol);file_command<<-"%scommand";\n',mydir),file=".Rprofile");
    if(make_win!=0)system("gnome-terminal --geometry=40x10 -- R");
    
}


##' @title send commands to simulators
#' @export
scom <- function(scom,runid=1){
    for(i in 1:length(runid))cat(sprintf("try(eval(parse(text=\'%s\')))\n",scom),file=paste0(file_command,runid[i],".R"));
}


#' @export
plot_var<-function(xid=1,yid=2){
    a$pparam$xid<<-xid;
    a$pparam$yid<<-yid;
    .simevol_func$plot_func();
}


##' @title change valuables for plotting
#' @export
pvar <- function(xid=1,yid=2,runid=1){
    for(i in 1:length(runid))cat(sprintf("plot_var(%d,%d)\n",xid,yid),file=paste0(file_command,runid[i],".R"));
}


#' @export
entercom <- function(){
    ##cat("enter command:");
    ##coms=scan("stdin",character(),n=1);
    coms=readline("enter command: ");
    if(coms!=""){
        print(try(eval(parse(text=coms),envir=.GlobalEnv)));
        entercom();
        }
    else{
        cat("command-mode ended\n");
    }
    
}


##' @title command prompt
#' @export
com <- function(runid=1){
    cat("entercom()\n",file=paste0(file_command,runid,".R"));
}


#' @export
halt <- function(void){
    print("simulation ended.");
    a$sparam$flag_halt<<-TRUE;
}


##' @title halt simulation
#' @export
hal <- function(runid=1){
    for(i in 1:length(runid))cat("halt()\n",file=paste0(file_command,runid[i],".R"));
}

#' @export
resetrange<-function(){
    a$pparam$xlim<<-c(NULL);
    a$pparam$ylim<<-c(NULL);
    .simevol_func$plot_func();
}

#' @export
setrr<-function(runid=1){
   for(i in 1:length(runid))cat("resetrange()",file=paste0(file_command,runid[i],".R"));
}



#' @export
setrange <-function(x0=NULL,x1=NULL,y0=NULL,y1=NULL){
    a$pparam$xlim<<-c(x0,x1);
    a$pparam$ylim<<-c(y0,y1);
    .simevol_func$plot_func();
}

#' @export
setrangex <-function(x0=NULL,x1=NULL){
    a$pparam$xlim<<-c(x0,x1);
    .simevol_func$plot_func();
}

#' @export
setrangey <-function(y0=NULL,y1=NULL){
    a$pparam$ylim<<-c(y0,y1);
    .simevol_func$plot_func();
}

#' @export
setr <- function(x0=NULL,x1=NULL,y0=NULL,y1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrange(%f,%f,%f,%f)",x0,x1,y0,y1),file=paste0(file_command,runid[i],".R"));
}

#' @export
setrx <- function(x0=NULL,x1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrangex(%f,%f)",x0,x1),file=paste0(file_command,runid[i],".R"));
}

#' @export
setry <- function(y0=NULL,y1=NULL,runid=1){
   for(i in 1:length(runid))cat(sprintf("setrangey(%f,%f)",y0,y1),file=paste0(file_command,runid[i],".R"));
}


#' @export
resetrange_sub<-function(id){
    a$pparam$sub_xlim[id]<<- list(NULL);
    a$pparam$sub_ylim[id]<<- list(NULL);
    .simevol_func$plot_func();
}

#' @export
setrange_sub <-function(x0=NULL,x1=NULL,id=1){
    a$pparam$sub_xlim[[id]]<<-c(x0,x1);
    .simevol_func$plot_func();
}


#' @export
cur.dev <- function(){
    v=as.numeric(dev.list()[length(dev.list())]);
    if(length(v)==0)v=1;
return(v);
}


#' @export
cpal <- function(palid=0,bgid=-2){
        if(class(palid)=="character"){
            a$pparam$resi$bgid <<- -1;
            a$pparam$resi$bg <<- palid;
            cat("\n bg color: ",a$pparam$resi$bg,"\n");
            
        }else{
            palname=c("1: ranbow", "2: heat.colors", "3: terrain.colors", "4: topo.colors", "5: cm.colors");
##            palname.squash=c('6: rainbow2', '7: jet', '8: grayscale', '9: heat', '10: coolheat', '11: blueorange', '12: bluered', '13: darkbluered');
##            palname=c(palname,palname.squash);
            
            
            if((palid==0)&&(bgid==-2)){
            cat("\nPalettes \n", paste(palname,collapse="\n "), "\n");
            }else{
                if(bgid==-2){
                    if(a$pparam$resi$bgid==-1)a$pparam$resi$bgid<<-1;
                }else{
                    a$pparam$resi$bgid<<-bgid;
                }
               
                if(palid==1)a$pparam$pal <<- rev(rainbow(100,end=0.7));
                if(palid==2)a$pparam$pal <<- heat.colors(100);
                if(palid==3)a$pparam$pal <<- terrain.colors(100);
                if(palid==4)a$pparam$pal <<- topo.colors(100);
                if(palid==5)a$pparam$pal <<- cm.colors(100);
                ##if(palid==6)a$pparam$pal <<- rainbow2(100);
                ##if(palid==7)a$pparam$pal <<- jet(100);
                ##if(palid==8)a$pparam$pal <<- grayscale(100);
                ##if(palid==9)a$pparam$pal <<- heat(100);
                ##if(palid==10)a$pparam$pal <<- coolheat(100);
                ##if(palid==11)a$pparam$pal <<- blueorange(100);
                ##if(palid==12)a$pparam$pal <<- bluered(100);
                ##if(palid==13)a$pparam$pal <<- darkbluered(100);
                cat("\n Pallete: ",palname[palid], " bgid:",a$pparam$resi$bgid,"\n");
            }    
        }
}

        
#' @export
pngout<-function(dev_id=as.numeric(dev.list()[length(dev.list())]),plotfile="testout",density=150,geometry=600,outeps=FALSE,prefix="./",show=TRUE,outpng=TRUE){
    
    dens=as.character(density);
    geom=as.character(as.integer(geometry));
    dev.set(dev_id);

    tmpf=paste0(a$mydir,"R_pngout_temp");
    
    dev.copy2eps(file=sprintf("%s.eps",tmpf));
   
     if(outeps){
         ##system(paste("cp .R_pngout_temp.eps ",prefix,plotfile,".eps",sep=""));
         system(sprintf("cp %s.eps %s%s.eps",tmpf,prefix,plotfile));
         cat("eps output:",paste(prefix,plotfile,".eps",sep=""),"\n")
     }
    
    if(outpng){
        system(sprintf("convert -density %sx%s -geometry %s  -background white -alpha remove %s.eps %s.png",dens,dens,geom,tmpf,tmpf));
        system(sprintf("mv %s.png %s%s.png",tmpf,prefix,plotfile));
        cat(sprintf("png output: %s%s.png \n",prefix,plotfile));
        if(show)system(sprintf("display %s%s.png&",prefix,plotfile));
    
    }
    system(sprintf("rm %s.eps",tmpf));
    
}




#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
#_/_/_/_/ MAIN FUNCTION FOR SIMULATION OF ADAPTIVE EVOLUTION _/_/_/_/#
#_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/#
##' This function simulates adaptive evolution by means of the oligomorphic stochastic model (OSM). The OSM assumes very rare mutations in comparison with the time scale of population dynamics, so that evolutionary dynamics can be described as a trait substitution sequence, engendered by repeated mutant invasions.
##'
##' The population dynamics triggerred by each mutant invasion is calculated with R-package "deSolve".
##' @title simulates adaptive evolution under given fitness function and mutation function.
##' @param phe list: phenotypes of coexisting residents.
##' @param en array: population densities of coexisting residents and environmental variables.
##' @param fitness function: fitness function.
##' @param mutate function: function for mutation.
##' @param pop_dynamics function: function for population dynamics. When this function is not given, the "fitness" function is used for calculation of population dynamics.
##' @return no direct output (the variable "a" contains all simulation data as well as parameters)
##' @author Hiroshi C. Ito
#' @export
simevol <- function(phe=a$phe,en=a$en,## state values
                    fitness=NULL,## functions
                    mutate=mutate0,
                    pop_dynamics=pop_dynamics0,
                    set_parms=NULL,
                    plot_func=plot_func0,
                    output=output0,
                    output_tree=output_tree0,
                    halt_func=NULL,
                    yorick_plot=NULL,
                    tmax=100000000, ## parameters for simulation and output
                    out_interval=10,
                    show_interval=10,
                    file_data="test.dat",
                    file_data_tree="test_tree.dat",
                    file_data_pid="test_pid.dat",
                    continue=FALSE,
                    runid=1,
                    runname="",
                    mydir=".simevol/",
                    fitness_contour=TRUE,## parameters for plotting
                    fitness_contour_phe=NULL,
                    plot_mask=NULL,
                    reset_win=TRUE,
                    show_subwin=TRUE,
                    trait_names=NULL,
                    env_names=NULL,
                    bgid=-1,
                    palid=1,
                    palfunc=NULL,
                    pparam=pparam0,
                    amp_invf=0.1,## parameters for mutant invasion
                    level_invf=0.02,
                    amp_invf_fix=FALSE,
                    m_rate=1.0,
                    m_sd=0.01,
                    n_mutant_init=1e-6,## parameters for population dynamics
                    edge_extinct=1e-8,
                    edge_fit=1e-13,
                    nrad=12,
                    drad=1.0,
                    divrad=2,
                    param_desolve=list(method="radau",hini=1e-4,hmax=1e9,rtol=1e-4,atol=1e-20)                      
                  ){


       ##.ee.append("simevol",environment())
 
    if(continue==FALSE){ 

        if(!file.exists(mydir))system(sprintf("mkdir %s",mydir));
        file_outcount=paste0(mydir,"outcount.dat");
        file_command=paste0(mydir,"command",runid,".R");
        
        timen=0.0;        
        outcount=1;
        
        nspe=length(phe[[1]]);
        n=en[1:nspe];
        pdim=length(phe);
        edim=length(en)-nspe;
        if(length(trait_names)==0)trait_names=names(phe);
        
        .simevol_func <<- list(fitness=fitness,pop_dynamics=pop_dynamics,mutate=mutate,set_parms=set_parms,plot_func=plot_func,output=output,output_tree=output_tree,halt_func=halt_func,palfunc=palfunc);
        
     
        pparam$plot_mask=plot_mask;
        pparam$trait_names=trait_names;
        pparam$env_names=env_names;
        pparam$fitness_contour=fitness_contour;
        pparam$fitness_contour_phe=fitness_contour_phe;
        pparam$resi$bgid=bgid;

        pparam=c(pparam,list(sub_xlim=vector("list",length=pdim),sub_ylim=vector("list",length=pdim)));

        
        traj=list(phe=phe,n=c(n),t=c(rep(0.0,nspe)),e=c(NULL),te=c(0.0));
        ##phe=c(phe,list(t=0.0,pid=0));
        ##tree=list(phe);
        ##tree_phe=c(phe,list(pid_par=-1));

        phe=c(phe,list(t=rep(0.0,nspe),pid=(seq(nspe)-1)));
        tree=list(phe);
        tree_phe=c(phe,list(pid_par=rep(-1,nspe)));
        
       
        
        sparam=list(m_rate=m_rate,
                    m_sd=m_sd,
                invf=c(0.1),
                fit_over=c(NULL),
                flag_halt=FALSE,
                edge_extinct=edge_extinct,
                edge_fit=edge_fit,
                amp_invf=amp_invf,
                level_invf=level_invf,
                amp_invf_fix=amp_invf_fix,
                n_mutant_init=n_mutant_init,
                mydir=mydir,
                file_command=file_command,
                file_outcount=file_outcount,
                outcount=outcount,
                out_interval=out_interval,
                show_interval=show_interval,
                nrad=nrad,
                drad=drad,
                divrad=divrad,
                param_desolve=param_desolve
                );

   

    a<<-list(
        phe=phe,
        en=en,
        n=n,
        e=en[(nspe+1):length(en)],
        timen=timen,
        ninv=0,
        inv=c(NULL),
        edim=edim,
        pdim=pdim,
        traj=traj,
        tree=tree,
        tree_phe=tree_phe,
        winid=c(2,3,0),
        sparam=sparam,
        pparam=pparam,
        file_data=file_data,
        file_data_tree=file_data_tree,
        file_data_pid=file_data_pid,

        runid=runid,
        runname=runname,
        show_subwin=show_subwin
        );
        options(scipen=100);
        write(c(1,0),file=a$file_data_pid,append=FALSE,ncolumns=2);
        options(scipen=0);
        ##cat("\n",file=a$file_data_pid,append=TRUE);


        cpal(palid);
        
        ##comwin();   

        
        if(a$edim>0)a$traj$e<<-c(a$traj$e,en[(nspe+1):length(en)]);
        
        if(file.exists(file_data))file.remove(file_data);
        
        .simevol_func$output(a$timen,phe,n);
        cat(a$sparam$outcount,file=a$sparam$file_outcount);
        a$sparam$outcount<<-a$sparam$outcount+1;
        
        if(class(yorick_plot)=="function"){
            yorick_plot();
            system("xterm -fn 7x14 -bg navy -fg white -e 'rlwrap -c ./yorick_idl_follow.sh'&");
        }

if(!pparam$omit_plot){
        if(reset_win || (length(dev.list())<2)){
            graphics.off();
            plot_init(pparam);
        }
}

        res=simpop(phe,en,.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad);
        

        phe=res$phe;
        nspe=length(phe[[1]]);
        en=res$en;
        n=en[1:nspe];

        a$n_all <<- data.frame(n=n,t=rep(a$timen,length(n)),ninv=rep(a$ninv,length(n)),pid=phe$pid);


        
        if(file.exists(file_command))file.remove(file_command);
    }
    else{
        a$sparam$flag_halt<<-FALSE;
    }

    
   
    for(t in 1:tmax){
        if(class(.simevol_func$halt_func)=="function").simevol_func$halt_func();
        if(a$sparam$flag_halt==T)break;

        if(file.exists(file_command)){
            ##source(file_command);
            sys.source(file_command,envir=.GlobalEnv);
            file.remove(file_command);
        }
        
        
        inv=invader(phe,en,mutate=.simevol_func$mutate,a$timen,amp_invf=a$sparam$amp_invf);
        
        a$timen<<-inv$timen;
        a$sparam$invf<<-c(a$sparam$invf,inv$f);
        a$inv <<- inv;
        a$ninv <<- a$ninv+1;
        
        state_new=simpop_invade(phe,en,inv,.simevol_func$pop_dynamics,set_parms=.simevol_func$set_parms,edge_extinct=a$sparam$edge_extinct,edge_fit=a$sparam$edge_fit,nrad=a$sparam$nrad,drad=a$sparam$drad,divrad=a$sparam$divrad);

        
        par_alive= which(state_new$phe$pid==inv$pid_par);
        a$tree_phe <<- add_phenotype(a$tree_phe,c(inv$phe,list(pid_par=inv$pid_par)));

        if(0){
        options(scipen=100);
        write(c(inv$phe$pid+1,inv$pid_par+1),file=a$file_data_pid,append=TRUE,ncolumns=2); ## to be improved so that pids in simevol and file_data_pid are the same.
        options(scipen=0);
              }  
        
        ##if(length(par_alive)>0){
        ##phe_par=geti(state_new$phe,par_alive);
        if((length(par_alive)>0)||(a$flag_new_group==1)){
            if(length(par_alive)>0){
            phe_par=geti(state_new$phe,par_alive);
            }else{
                phe_par=geti(phe,which(phe$pid==inv$pid_par))
            }
           ## print(phe_par);
           ## print(phe);
            ##phe_par=phe;
            
            ##print(inv$phe$pid);
            ##print(state_new$phe$pid);            
            ##print(par_alive);
            
            
            ##print(phe_par);
            a$tree <<- c(a$tree,list(add_phenotype(phe_par,inv$phe)));
                
        }
        else{
            for(i in 1:length(a$tree)){
                buf=a$tree[[i]];
                ##cat("pars:",buf$pid);
                ##cat("inv:",inv$phe$pid,"inv_par:",inv$pid_par,"\n");
                if(sum(buf$pid[length(buf$pid)]==inv$pid_par)>0){
                    ##cat("connect\n");
                    a$tree[[i]]<<-add_phenotype(a$tree[[i]],inv$phe);
                }
            }
        }
        
        for(i in 1:length(a$tree)){
            buf=a$tree[[i]];
            buf=geti(buf,length(buf$pid));
            if(buf$pid>0){
                if(sum(state_new$phe$pid==buf$pid)==0){
                    ##print(buf$pid);
                    ##   a$tree_phe$te[(buf$pid+1)]<<-a$timen;
                    buf$t=a$timen;
                    buf$pid=-1; ## -1 means extinction
                    a$tree[[i]]<<-add_phenotype(a$tree[[i]],buf);
                }
                
            }
        }
        
   

        
        phe=state_new$phe;
        nspe=length(phe[[1]]);
        en=state_new$en;
        n=en[1:nspe];
            
        a$phe<<-phe;
        a$en<<-en;
        a$n<<-n;
        a$e<<- en[(nspe+1):length(en)];

        a$n_all <<- rbind(a$n_all,data.frame(n=n,t=rep(a$timen,length(n)),ninv=rep(a$ninv,length(n)),pid=phe$pid));
        
        level_invf=a$sparam$level_invf;
        if(amp_invf_fix==FALSE){
            a$sparam$amp_invf<<-max(level_invf,level_invf/mean(a$sparam$invf[max((t-100),1):length(a$sparam$invf)]));
        }

        
        ##.simevol_func$output(a$timen,phe,n);
                
        ##cat(a$sparam$outcount,file=a$sparam$file_outcount);
        ##a$sparam$outcount<<-a$sparam$outcount+1;

        
        if((a$sparam$flag_halt==T)||(a$ninv%%a$sparam$out_interval==0)){
          ##  .simevol_func$output(a$timen,phe,n);
                
        ##cat(a$sparam$outcount,file=a$sparam$file_outcount);
        ##a$sparam$outcount<<-a$sparam$outcount+1;

                    
            a$traj$phe<<-add_phenotype(a$traj$phe,phe);
            if(0){
                .simevol_func$output_tree(a$file_data_tree);
            }
            
            if(a$edim>0)a$traj$e<<-rbind(a$traj$e,en[(nspe+1):length(en)]);
            a$traj$n<<-c(a$traj$n,n);
            a$traj$t<<-c(a$traj$t,rep(a$timen,nspe));
            a$traj$te<<-c(a$traj$te,a$timen);
            
        }

        if(a$sparam$show_interval>500){
           ## if(a$ninv%%500==0)cat("ninv: ",a$ninv, "  nspe:", length(a$n),"  ymax: ", max(phe$y), " \n");
        }
        if((a$ninv<=1)||(a$sparam$flag_halt==T)||(a$ninv%%a$sparam$show_interval==0)){
            runname1="";
            if(runname!="")runname1=paste0("\"",runname,"\"");
            cat("runid:",runid,runname1,"time:",a$timen, " residents:",nspe, " invasion:", a$ninv,"amp:",a$sparam$amp_invf,"fit_over:",length(a$sparam$fit_over),"irad:",state_new$irad,"\n");
            
            if(!pparam$omit_plot).simevol_func$plot_func(out=TRUE,analyze=TRUE);
            if(a$sparam$flag_halt){
                cat("\n last:",a$ninv,"\n");
                
            }
            
        }
##        if(a$sparam$flag_halt==T)break;
    }
    
    if((length(a$res$ninv)==0)||(a$res$ninv < a$ninv)){
        
       
        if(!pparam$omit_plot).simevol_func$plot_func(out=TRUE,analyze=TRUE);
       }
}


