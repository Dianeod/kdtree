\l kdtree.q

\d .kd

imax:{x?max x};
imin:{x?min x};

cure.clust:{[sample;numR;com;diml;kd]
 kd2:select from kd where valid; /only count the valid rows
 
 j:first select from kd2 where closDist=min closDist; /get the closest clusters(only want details of one of them)

 j2:select clust from kd2 where idx=first j`closIdx; /get cluster no.

 old:select from kd2 where clust in ((j2`clust),j`clust); /get all points from that cluster
 j0:(exec initi from kd2 where closIdx in old`idx); /get initial indexes from pts in the cluster who had these pts as their closestIdx

 mean:avg pts:sample idxs:(distinct raze old`clustIdx); /get the mean of the clusters

 deleteClust:tree.deleteN/[kd;idxs]; /delete the closest points from the tree

 maxFromMean:idxs imax sum each{x*x}mean-/:pts; /get the maximum points from the mean 

 rep:distinct sample sami:numR{[sample;idxs;x] x,maxIdx imax {[sample;maxMean;x] min{[sample;x;maxMean]sum x*x:sample[maxMean]-sample[x]
     }[sample;x] each maxMean}[sample;x]each maxIdx:idxs except x}[sample;idxs]/maxFromMean; /get the representative points from the cluster by choosing the most spread out points

 rep:(rep*1-com)+\:com*mean; /move it towards the centre
 a:{[rep;x]min (sum each x*x:rep-\:x)}[rep]each (exec rep from deleteClust where valid); /get the distances from the new rep points to every other rep point in the tree 
 j3:(exec idx from deleteClust where valid)[n:where a<(exec closDist from deleteClust where valid)]; /get the index of the pts that are now closer to the new rep pts than their previous closest pts
 b: (exec idx from deleteClust where valid)where a=c:min a; /get the idx of the closest pnt to any of the new rep points
 insertClust:tree.insertKd[diml]/[deleteClust;rep;sami;(first idxs)]; / insert the new rep pts into the tree
 insertClust2:$[0=count[j3];insertClust;{[n;j3;a;insertClust;x] update closDist:a(n x)
  ,closIdx:enlist (max (insertClust`idx)) from
    insertClust where idx=j3(x)}[n;j3;a]/[insertClust;til count j3]]; /update the tree if the new rep pts are closest to any other clusters

 clustDist:update closIdx:first b,closDist:c from insertClust2 where clust=first idxs; /update the closest pt to the new rep pt
 j5:exec idx from insertClust2 where clust=first idxs,valid; /get idx of all new inserted rep pts
 j4:(exec idx from clustDist where initi in j0,valid) except j5; /the index of pts who use to have the rep pts as their closestidx
 
 recalc:tree.distC/[clustDist;j4]; /recalc these new distances
 
 insertIdx:enlist idxs; /get the initial indexes of pts in the new merged cluster

  kd:{[insertIdx;kd;x] update clustIdx:(insertIdx) from kd where initi=x,valid}[insertIdx]/[recalc;sami]; /update these into kdtree
  
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

