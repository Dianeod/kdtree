/initially insert cluster into tree

\d .kd

tree.insertKd1:{[diml;kd;sample2;L;clust]
   /check if its to the left or right of initial cluster in tree

 dir:{(first x[0]`idx)>=0}{[sample2;kd;x] a:x[0];
  $[first (sample2[first x[0]`dim]>(first x[0]`rDim));
  i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep,rDim from kd where right=1,parent=first a`idx;
  i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep,rDim from kd where left=1,parent=first a`idx];
  (i;a)}[sample2;kd]/(i;i:exec from kd where parent=0,left=0,right=0); /insert cluster into the tree by looking at splitting dimension of each node &going left or right


 i:dir[0];
 a:dir[1];
 dim:diml (first a`dim)+1; /get its new splitting dimension


 root:kd upsert flip update idx:(max kd`idx)+1,initi:L,clust:clust,rep:enlist sample2,
    rDim:enlist sample2[dim],dim:enlist dim,valid:1b,
    parent:enlist first a`idx from exec left,right from i; /update the info of new node into tree
 root}

/insert cluster into tree

tree.insertKd:{[diml;kd;sample2;L;cl]
   /check if its to the left or right of initial cluster in tree

 dir:{(first x[0]`idx)>=0}{[sample2;kd;x] a:x[0];
  $[first (sample2[first x[0]`dim]>(first x[0]`rDim));
  i:update left:enlist 0,right:enlist 1 from exec idx,dim,rep,rDim from kd where right=1,parent=first a`idx,valid;
  i:update left:enlist 1,right:enlist 0 from exec idx,dim,rep,rDim from kd where left=1,parent=first a`idx,valid];
  (i;a)}[sample2;kd]/(i;i:exec from kd where parent=0,left=0,right=0); /insert cluster into the tree by looking at splitting dimension of each node &going left or right

 i:dir[0];
 left1:i`left;
 right1:i`right;
 a:dir[1];
 dim:diml (first a`dim)+1; /get its new splitting dimension
 b:exec idx from kd where valid=0b;
 root:update idx:(max kd`idx)+1, initi:L, clust:cl,rep:enlist sample2,
    valid:enlist 1b,dim:dim, rDim:sample2[dim],left:left1,right:right1,
    parent:a`idx from kd where idx=first b;  /update the info of new node into tree

 root}

/Calculating distances between clusters using kd tree.


tree.distC:{[kd;pt]

 distCalc:{[kd;query;bestD]

 X:bestD[3]; / nodes that were already searched, not to be searched again

 cl:bestD[5]; / what cluster the node belongs to, as not to search any points in the same cluster
 n:bestD[1]; /nodes to search
 a:n where {[cl;kd;x](first exec clust from kd where idx=x,valid)<>cl}[cl;kd]each n; /nodes to search that arent in the same cluster

 newD:imins,a i?imins:min i:{[kd;query;x]
    sum m*m:(raze exec rep from kd where idx=x,valid)-query}[kd;query] each a; /get minimum dist of all searched nodes

  $[(newD[0]<bestD[0])&(count a)<>0;(bestD[0]:newD[0];bestD[2]:newD[1])
     ;(bestD:bestD)]; /if new dist is less than current best dist, then that becomes new best dist

  axisD:(raze {[kd;bestD;query;x] $[(m*m:
    (first exec rDim from kd where idx=x,valid)-
    query(first exec dim from kd where idx=x,valid))<=bestD[0];
    (exec idx from kd where parent=x,valid),exec parent from kd where idx=x,valid;
    (query(first exec dim from kd where idx=x,valid))<
    first exec rDim from kd where idx=x,valid;
    (exec idx from kd where parent=x,left=1,valid),exec parent from kd where idx=x,valid;
    (exec idx from kd where parent=x,right=1,valid),exec parent from kd where idx=x,valid]
    }[kd;bestD;query]each n)except bestD[3]:X,n; /get dists between node and search pts based on splitting dimension.
    /if =< than best Dist, than search the children of that node &parents.
    /Go up the tree and to the left/right based on whether that pt is <or> search pt

 (bestD[0];distinct axisD;bestD[2];bestD[3];bestD[4];cl)};
 dist:{[distCalc;kd;X] ({(count x[1])<>0}distCalc[kd;
    raze exec rep from kd where idx=X,valid]/(0W;
    (raze (exec idx from kd where parent=X,valid),exec parent from kd where idx=X,valid)except X;X;X;X
    ;first exec clust from kd where idx=X))
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

