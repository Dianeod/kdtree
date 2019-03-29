\d .kd

tree.createTree:{[sample;diml]
 root:flip `idx`initi`rep`dir`dim`parent`clust`clustIdx`valid!(0;0;enlist @[sample;0];2;0;0;0;enlist (0 0);1b); 			/insert first cluster
 kds:tree.insertKd[diml]/[root;1_sample;cl;cl:1_til count sample]; 								/insert the rest of the clusters
 kds:update clustIdx:enlist each til count[sample] from kds; 									/insert the cluster indices
 tree.distC/[kds;kds`initi]
 }

/insert cluster into tree
tree.insertKd:{[diml;kd;samp;L;cl] 												/check if its to the left or right of initial cluster in tree
 dirn:{0<=first @[x;0]`idx}{[samp;kd;nn]a:@[nn;0];dir1:$[samp[first a`dim]>(raze a`rep)[first a`dim];1;0];
  i:select idx,dim,rep from kd where dir=dir1,parent=first a`idx,valid;                                      			/a is the previous pt, while i is the next pt in the tree to do split on
   (i;a;dir1)}[samp;kd]/(i;i:@[kd;0]); 												/insert cluster into the tree by looking at splitting dimension of each node
 par:@[dirn;1];
 $[0=count select from kd where not valid;kd upsert ([]idx:max[kd`idx]+1;initi:L;clust:L;rep:enlist samp;
  dim:diml par[`dim]+1;valid:1b;parent:par`idx;dir:@[dirn;2]);
  update idx:(max kd`idx)+1,initi:L,clust:cl,rep:enlist samp,valid:1b,dim:diml par[`dim]+1,
  dir:@[dirn;2],parent:par`idx from kd where idx=first exec idx from kd where not valid] 
 } 																/update the info of new node into tree

/Calculating distances between clusters using kd tree
tree.distC:{[kd;pt] 
 distCalc:{[kd;query;cl;bestD]
 newn:nn where{[cl;kd;x](first exec clust from kd where idx=x,valid)<>cl}[cl;kd]each nn:@[bestD;1]; 				/nodes to search that arent in the same cluster
 newD:imins,newn ii?imins:min ii:{[kd;query;x]sum m*m:(first exec rep from kd where idx=x,valid)-query}[kd;query] each newn; 	/get minimum dist of all searched nodes
 $[(@[newD;0]<@[bestD;0])&count[newn]<>0;(bestD[0]:@[newD;0];bestD[2]:@[newD;1]);]; 						/if new dist is less than current best dist, then that becomes new best dist
 axisD:(raze{[kd;bestD;query;nn]ll:exec dim,rep,parent from kd where idx=nn,valid;dirn:$[(qdim:query[dd])<rdim:(first ll[`rep])
  [dd:first ll`dim];0;1];$[@[bestD;0]>=m*m:rdim-qdim;exec idx from kd where parent=nn,valid;exec idx 
  from kd where parent=nn,dir=dirn,valid],ll`parent}[kd;bestD;query]each nn)except bestD[3]:@[bestD;3],nn;			/bestD[3]=nodes already searched
 (@[bestD;0];distinct axisD;@[bestD;2];@[bestD;3];@[bestD;4])}; 								/get dists between node and pts based on split dim.if=<than best Dist than 																    /search the leafs.Go up the tree and to the lft/rght depending on position
 dist:{[distCalc;kd;pts] idpts:select parent,clust,rep from kd where idx=pts,valid;{0<>count @[x;1]}distCalc[kd;first idpts`rep;
 first idpts`clust]/(0W;(raze idpts[`parent],raze exec idx from kd where parent=pts,valid)except pts;pts;pts;pts
  )}[distCalc;kd;pt];update closDist:@[dist;0],closIdx:@[dist;2] from kd where idx=pt 						/update new closDist and idx in tree
 }

/ delete from kd tree
tree.deleteN:{[kd;X]
 delN:{kd:@[x;0];X:@[x;1]; 													/kdtree,pt to be deleted
  delNode:first exec dim from kd where idx=X,valid; 										/details of deleted pt
  dirn:$[0=count ll:exec idx from kd where dir=1,parent=X,valid;first exec idx from kd where parent=X,dir=0,valid;first ll]; 
  mindim:raze{[kd;x]0<>count exec idx from kd where parent=first x,valid}[kd]{[kd;x]raze exec idx from kd where parent in
   x,valid}[kd]\dirn;	 													/ if has no right child then left child replaces  /X=pt to be deleted
  newP:mindim imin raze ({[kd;x]first exec rep from kd where idx=x,valid}[kd]each mindim)[;delNode]; 				/the min value of rep pts from mindim based on splitting dimension
  newNode:select from kd where idx=newP,valid; 											/get info from kdtree of newP
  tree:update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,clust:newNode`clust,clustIdx:newNode`clustIdx,
    closIdx:newNode`closIdx from kd where idx=X,valid; 										/newNode replaces the deleted pt,but dim and idx stay the same as before
  (update closIdx:X from tree where closIdx=first newNode`idx,valid;newP) 							/any pt that had newP as closest idx has to update closIdx to new position
  };
 delCl:{0<>count select from @[x;0] where parent=@[x;1],valid}delN/(kd;first exec idx from kd where initi=X,valid); 		/repeat this until reach a node with no children
  update valid:0b from first delCl where idx=last delCl} 									/delete the node with no children
