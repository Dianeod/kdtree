\l kdtree.q

\d .kd

imax:{x?max x};
imin:{x?min x};

cure.clust:{[sample;numR;com;diml;kd]
 kd2:select from kd where valid;
 j:first select from kd2 where closDist=min closDist;

 j2:select clust from kd2 where idx=first j`closIdx;

 old:select from kd2 where clust in ((j2`clust),j`clust);
 j0:(exec initi from kd2 where closIdx in old`idx);

 mean:avg pts:sample idxs:(distinct raze old`clustIdx);

 deleteClust:tree.deleteN/[kd;idxs];

 maxFromMean:idxs imax sum each{x*x}mean-/:pts;

 rep:distinct sample sami:numR{[sample;idxs;x] x,maxIdx imax {[sample;maxMean;x] min{[sample;x;maxMean]sum x*x:sample[maxMean]-sample[x]
     }[sample;x] each maxMean}[sample;x]each maxIdx:idxs except x}[sample;idxs]/maxFromMean;

 rep:(rep*1-com)+\:com*mean;
 a:{[rep;x]min (sum each x*x:rep-\:x)}[rep]each (exec rep from deleteClust where valid);
 j3:(exec idx from deleteClust where valid)[n:where a<(exec closDist from deleteClust where valid)];
 b: (exec idx from deleteClust where valid)where a=c:min a;
 insertClust:tree.insertKd[diml]/[deleteClust;rep;sami;(first idxs)];
 insertClust2:$[(count j3)=0;insertClust;{[n;j3;a;insertClust;x] update closDist:a(n x)
  ,closIdx:enlist (max (insertClust`idx)) from
    insertClust where idx=j3(x)}[n;j3;a]/[insertClust;til count j3]];

 clustDist:update closIdx:first b,closDist:c from insertClust2 where clust=first idxs;
 j5:exec idx from insertClust2 where clust=first idxs,valid;
 j4:(exec idx from clustDist where initi in j0,valid) except j5;

 recalc:tree.distC/[clustDist;j4];
 insertIdx:enlist idxs;

  kd:{[insertIdx;kd;x] update clustIdx:(insertIdx) from kd where initi=x,valid}[insertIdx]/[recalc;sami]; /update into kdtree

  kd}


/Create tree, search init nearest neighbours
cure.createTree:{[sample;diml]

 root:flip `idx`initi`rep`left`right`dim`parent`rDim`clust`clustIdx`valid!
   (0;0;enlist sample[0];enlist 0;enlist 0;enlist 0;enlist 0;enlist sample[0;0];0;enlist (0 0);1b); /insert first cluster
 
   
 kds:tree.insertKd1[diml]/[root;(1_sample);1_til count sample;1_til count sample]; /insert the rest of the clusters

 kds:update clustIdx:enlist each til count sample from kds; /insert the cluster indices

 kds:tree.distC/[kds;kds`initi]; /get closest cluster to each cluster
 kds
 }

cure.cure:{[sample;numR;com;numCl]
 diml:(til count[first sample]),0;
 cureTab:{[numCl;kd] (count distinct (select from kd where valid)`clust)>numCl}[numCl] cure.clust[sample;numR;com;diml]/cure.createTree [sample;diml];
 distinct (select from cureTab where valid)`clustIdx}

