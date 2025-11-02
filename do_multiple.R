
## Test functions written by Hiroshi C. Ito (2024)
## Email: hiroshibeetle@gmail.com
##
## copyright (C) 2024 Hiroshi C. Ito
## This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License Version 2 as 
## published by the Free Software Foundation.
## http://www.r-project.org/Licenses/


do_simulation=TRUE; 

njob=10; 
nsim=10; 

data_type=1;

do_unit_simulation_R="ddflw.R";
if(data_type==1)do_unit_simulation_R="ddflwd.R";


runids= 1:nsim; 

seed_offset=10; 
seeds= runids + seed_offset; 

halt_edge_y=10; 

odir0="tmp";
if(data_type==0)dtype="";
if(data_type==1)dtype="_discrete";

odir=paste0(odir0,dtype);

ofile_base="trial"; 
ofiles=paste0(odir,"/",ofile_base,runids,".RData");

command_do_sim=rep("",nsim);
    
if(dir.exists(odir)){

  
        file.remove(list.files(path=odir,pattern=".RData",full.names=TRUE));        
}else{

  
        dir.create(odir);
    }
    

    
    for(i in 1:nsim){

     
        command_do_sim[i]=sprintf("Rscript %s runid=%d seed=%d halt_edge_y=%f file_RData=\\\"%s\\\"",do_unit_simulation_R,i,seeds[i],halt_edge_y,ofiles[i])
    }

print(command_do_sim); 

write(command_do_sim,file=sprintf("%s/tmp_do_multiple_simulation.sh",odir)); 

command_do_multiple_sim=sprintf("cat %s/tmp_do_multiple_simulation.sh | parallel -t -j %d",odir,njob); 

if(do_simulation)system(command_do_multiple_sim); 



