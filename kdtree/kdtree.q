/initially insert cluster into tree

\d .kd

tree.insertKd1:{[diml;kd;samp;L;cl]
   /check if its to the left or right of initial cluster in tree

 dir:{(first x[0]`idx)>=0}{[samp;kd;nn] a:nn[0];
  $[first (samp[first nn[0]`dim]>(first nn[0]`rDim));
  i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep,rDim from kd where right=1,parent=first a`idx;
  i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep,rDim from kd where left=1,parent=first a`idx];
  (i;a)}[samp;kd]/(i;i:exec from kd where parent=0,left=0,right=0); /insert cluster into the tree by looking at splitting dimension of each node &going left or right


 ii:dir[0];
 par:dir[1];
 dim:diml (first par`dim)+1; /get its new splitting dimension


 root:kd upsert flip update idx:(max kd`idx)+1,initi:L,clust:cl,rep:enlist samp,
    rDim:enlist samp[dim],dim:enlist dim,valid:1b,
    parent:enlist first par`idx from exec left,right from ii; /update the info of new node into tree
 root}

/insert cluster into tree

tree.insertKd:{[diml;kd;samp;L;cl]
   /check if its to the left or right of initial cluster in tree

 dir:{(first x[0]`idx)>=0}{[samp;kd;nn] a:nn[0];
  $[first (samp[first nn[0]`dim]>(first nn[0]`rDim));
  i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep,rDim from kd where right=1,parent=first a`idx,valid; /a is the previous pt, while i is the next pt in the tree to do split on
  i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep,rDim from kd where left=1,parent=first a`idx,valid];
  (i;a)}[samp;kd]/(i;i:exec from kd where parent=0,left=0,right=0); /insert cluster into the tree by looking at splitting dimension of each node &going left or right

 ii:dir[0];
 par:dir[1];
 dim:diml (first par`dim)+1; /get its new splitting dimension
 nval:exec idx from kd where not valid;
 root:update idx:(max kd`idx)+1, initi:L, clust:cl,rep:enlist samp,
    valid:enlist 1b,dim:dim, rDim:samp[dim],left:ii[`left],right:ii[`right],
    parent:par`idx from kd where idx=first nval;  /update the info of new node into tree

 root}

/Calculating distances between clusters using kd tree.


tree.distC:{[kd;pt]

 distCalc:{[kd;query;bestD]

 prevn:bestD[3]; / nodes that were already searched, not to be searched again

 cl:bestD[5]; / what cluster the node belongs to, as not to search any points in the same cluster
 nn:bestD[1]; /nodes to search
 newn:nn where {[cl;kd;x](first exec clust from kd where idx=x,valid)<>cl}[cl;kd]each nn; /nodes to search that arent in the same cluster

 newD:imins,newn ii?imins:min ii:{[kd;query;x]
    sum m*m:(raze exec rep from kd where idx=x,valid)-query}[kd;query] each newn; /get minimum dist of all searched nodes

  $[(newD[0]<bestD[0])&(count newn)<>0;(bestD[0]:newD[0];bestD[2]:newD[1])
     ;(bestD:bestD)]; /if new dist is less than current best dist, then that becomes new best dist

  axisD:(raze {[kd;bestD;query;nn] $[(m*m:
    (first exec rDim from kd where idx=nn,valid)-
    query(first exec dim from kd where idx=nn,valid))<=bestD[0];
    (exec idx from kd where parent=nn,valid),exec parent from kd where idx=nn,valid;
    (query(first exec dim from kd where idx=nn,valid))<
    first exec rDim from kd where idx=nn,valid;
    (exec idx from kd where parent=nn,left=1,valid),exec parent from kd where idx=nn,valid;
    (exec idx from kd where parent=nn,right=1,valid),exec parent from kd where idx=nn,valid]
    }[kd;bestD;query]each nn)except bestD[3]:prevn,nn; /get dists between node and search pts based on splitting dimension.
    /if =< than best Dist, than search the children of that node &parents.
    /Go up the tree and to the left/right based on whether that pt is <or> search pt

 (bestD[0];distinct axisD;bestD[2];bestD[3];bestD[4];cl)};
 dist:{[distCalc;kd;pts] ({(count x[1])<>0}distCalc[kd;
    raze exec rep from kd where idx=pts,valid]/(0W;
    (raze (exec idx from kd where parent=pts,valid),exec parent from kd where idx=pts,valid)except pts;pts;pts;pts
    ;first exec clust from kd where idx=pts))
    }[distCalc;kd;pt];

    kdC:update closDist:dist[0],closIdx:dist[2] from kd where idx=pt; /update new closDist and idx in tree
    kdC}



/ delete from kd tree

tree.deleteN:{[kd;X]
 n:first exec idx from kd where initi=X,valid;
 delN:{
  kd:x[0]; /kdtree
  X:x[1]; /point to be deleted
  delNode:select from kd where idx=X,valid; /details of deleted pt
  axis:delNode`dim; /splitting dim

  mindim:$[(count exec idx from kd where right=1,parent=X,valid)=0;
    raze {[kd;x](count exec idx from kd where
       parent=first x,valid)<>0}[kd]{[kd;x]raze exec idx from kd where parent in x,valid}[kd]\
    first exec idx from kd where parent=X,left=1,valid; / if has no right child then left child replaces
    raze {[kd;x] (count exec idx from kd where
       parent=first x,valid)<>0}[kd]{[kd;x]raze exec idx from kd where parent in x,valid}[kd]\
    first exec idx from kd where parent=X,right=1,valid]; / get all the right children if there

  newP:mindim repPt?min repPt:(raze{[kd;x] exec rep from kd where idx=x,valid}[kd]each mindim)axis; /the min value of rep pts from mindim based on splitting dimension
  newNode:select from kd where idx=newP,valid; /get info from kdtree of newP
  tree:update rep:newNode`rep,initi:newNode`initi,closDist:newNode`closDist,
    clust:newNode`clust,clustIdx:newNode`clustIdx,rDim:(first newNode`rep)axis,
    closIdx:newNode`closIdx from kd where idx=X,valid; /newNode replaces the deleted pt info in the tree, but dim and idx stay the same as deleted pt


  tree:update closIdx:X from tree where closIdx=first newNode`idx,valid; /any pt that had the replacing pt as its closest idx has to update closIdx to its new value in the tree


    (tree;newP)};

 delCl:{(count select from (first x) where parent=last x,valid)<>0}delN/(kd;n); /repeat this until reach a node with no children


 update valid:0b from first delCl where idx=last delCl} /delete the node with no children

