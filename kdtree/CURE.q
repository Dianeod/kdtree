\l kdtree.q

\d .kd

imax:{x?max x};
imin:{x?min x};

cure.clust:{[sample;numR;com;diml;kd]
 kd2:select from kd where valid;                                                                    		/only count the valid rows
 j:first select from kd2 where closDist=min closDist; 						 		/get the closest clusters(only want details of one of them)
 j2:first exec clust from kd2 where idx=first j`closIdx; 							/get cluster info.
 old:select from kd2 where clust in (j2,j`clust); 								/get all points from that cluster
 j0:(exec initi from kd2 where closIdx in old`idx)except old`initi;						/get initial indexes from pts in the cluster who had these pts as their closestIdx
 mean:avg pts:sample idxs:distinct raze old`clustIdx;	 							/get the mean of the clusters
 deleteClust:tree.deleteN/[kd;idxs]; 										/delete the closest points from the tree
 maxFromMean:idxs imax sum each{x*x}mean-/:pts; 								/get the maximum points from the mean 
 rep:distinct sample sami:numR{[samp;idxs;nn]nn,maxI imax{[samp;nn;maxI]min{[samp;nn;maxI]sum x*x:
  samp[maxI]-samp[nn]}[samp;maxI]each nn}[samp;nn]each maxI:idxs except nn}[sample;idxs]/maxFromMean;  	 	/get the representative points from the cluster by choosing the most spread out points
 rep:(rep*1-com)+\:com*mean; 											/move it towards the centre
 repdist:{[rep;x]min sum each x*x:rep-\:x}[rep]each exec rep from deleteClust where valid; 			/get the distances from the new rep points to every other rep point in the tree 
 j3:(exec idx from deleteClust where valid)n:where repdist<exec closDist from deleteClust where valid;		/get the index of pts that are now closer to the new rep pts than their previous closest pts
 nclos:(exec idx from deleteClust where valid)where repdist=brep:min repdist; 					/get the idx of the closest pt to any of the new rep points
 insertCl:tree.insertKd[diml]/[deleteClust;rep;sami;first idxs]; 						/ insert the new rep pts into the tree
 $[0<count[j3];insertCl:{[rd;insertCl;n;j3]update closDist:rd[n],closIdx:enlist max[insertCl`idx] 
   from insertCl where idx=j3}[repdist]/[insertCl;n;j3];]; 							/update the tree if the new rep pts are closest to any other clusters 
 clustD:update closIdx:first nclos,closDist:brep from insertCl where clust=first idxs; 				/update the closest pt to the new rep pt
 recalc:tree.distC/[clustD;exec idx from clustD where initi in j0,valid]; 					/recalc the new distances of pts who use to have the rep pts as their closestidx
 {[insertIdx;kd;x] update clustIdx:insertIdx from kd where initi=x,valid}[enlist idxs]/[recalc;sami] 		/get the initial indexes of pts in the new merged clust&update these into kdtree
 }

cure.cure:{[sample;numR;com;numCl]
 cureTab:{[numCl;kd] (count distinct (select from kd where valid)`clust)>numCl}[numCl] cure.clust[sample;numR;com;diml]/tree.createTree[sample;diml:(til count[first sample]),0];
 distinct (select from cureTab where valid)`clustIdx}
