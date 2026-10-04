make_pairs <- function(spec) {
  X <- t(vapply(spec,function(s) {x<-numeric(4);x[s]<-c(6,3,1);x},numeric(4)))
  colnames(X)<-LETTERS[1:4];rownames(X)<-paste0('s',seq_len(nrow(X)));X
}

test_that('one/two/three configurations produce zero/one/three edges and one witnessed face', {
  for(n in 1:3) {
    X<-make_pairs(list(c(1,2,3),c(1,3,2),c(2,3,1))[seq_len(n)])
    g<-linf.udcst.graph(X)
    expect_equal(nrow(g$edges),c(0,1,3)[n]);expect_equal(nrow(g$faces),as.integer(n==3))
    expect_true(all(g$edges$length==1));expect_equal(sum(g$incidence$count),n)
  }
})
test_that('full denominators survive vertex filtering and support thresholds preserve isolates', {
  X<-make_pairs(c(rep(list(c(1,2,3)),5),rep(list(c(1,3,2)),3),list(c(2,3,1))))
  g<-linf.udcst.graph(X);h<-linf.udcst.graph(X,states=c('1,2','1,3'))
  expect_equal(h$edges$n.triple,9);expect_equal(h$edges$i.a,5);expect_equal(h$edges$i.b,3)
  expect_equal(h$edges$jaccard.mean,(5/9+3/9)/2)
  expect_equal(h$edges,g$edges[1,,drop=FALSE]);expect_equal(nrow(h$faces),0)
  z<-linf.udcst.graph(X,min.support=4)
  expect_equal(nrow(z$nodes),3);expect_equal(nrow(z$edges),0);expect_equal(z$metadata$components,3)
  expect_equal(z$metadata$isolates,3);expect_equal(nrow(z$faces),0)
  empty<-linf.udcst.graph(X,states=character());expect_equal(nrow(empty$nodes),0)
  expect_equal(nrow(empty$edges),0);expect_equal(empty$metadata$components,0)
})
test_that('triangles from different witnesses have no face', {
  X<-make_pairs(list(c(1,2,3),c(1,3,2),c(1,2,4),c(1,4,2),c(1,3,4),c(1,4,3)))
  g<-linf.udcst.graph(X);expect_equal(nrow(g$edges),3);expect_equal(nrow(g$faces),0)
  expect_equal(length(unique(g$edges$triple)),3)
})
test_that('missing witness configuration distinguishes shared from observed edges', {
  X<-make_pairs(list(c(1,2,3),c(1,3,4)))
  expect_equal(nrow(linf.udcst.graph(X)$edges),0)
  g<-linf.udcst.graph(X,adjacency='shared_feature');expect_equal(nrow(g$edges),1)
  expect_equal(g$edges$min.support,0);expect_equal(g$metadata$components,1)
})
test_that('coverage retains ties and permits unobserved dictionary features', {
  X<-make_pairs(list(c(1,2,3),c(1,3,2),c(2,3,1)))
  g<-linf.udcst.graph(X,coverage=60);expect_equal(nrow(g$nodes),3);expect_equal(g$metadata$coverage,100)
  expect_error(linf.udcst.graph(X,coverage=0),'coverage')
  expect_error(linf.udcst.graph(X,coverage=60,n0=2),'only one')
  expect_error(linf.udcst.graph(X,states='1,4'),'observed')
  expect_error(linf.udcst.graph(X,min.support=0),'positive')
  expect_error(linf.udcst.graph(X,sample.ids=rep('same',3)),'unique')
})
test_that('dictionary, row order and pure ranks are checked for fitted inputs', {
  X<-make_pairs(list(c(1,2,3),c(1,3,2),c(2,3,1)))
  a<-linf.udcsts(X,2,n0=99,low.freq.policy='filter');b<-linf.udcsts(X,3,n0=99)
  g<-linf.udcst.graph(pair.fit=a,triple.fit=b);expect_equal(g$metadata$reference.n,3)
  expect_equal(g$membership$sample.id,rownames(X))
  bad<-b;bad$pure.id[1]<-'1,2,4';expect_error(linf.udcst.graph(pair.fit=a,triple.fit=bad),'disagree')
  bad<-b;bad$feature.ids<-rev(bad$feature.ids);expect_error(linf.udcst.graph(pair.fit=a,triple.fit=bad),'dictionary')
  expect_error(linf.udcst.graph(pair.fit=a,triple.fit=linf.udcsts(X[3:1,],3)),'row order')
  expect_error(linf.udcst.graph(X,pair.fit=a),'not both')
  expect_error(linf.udcst.graph(pair.fit=a),'both fitted')
})
test_that('ties, zero thirds, row permutations and sparse input preserve identities', {
  X<-rbind(c(1,1,0,0),c(0,0,1,0),c(3,1,1,1));colnames(X)<-LETTERS[1:4];rownames(X)<-c('a','b','c')
  g<-linf.udcst.graph(X);expect_equal(g$membership$pair,c('1,2','1,3','1,2'))
  expect_equal(g$membership$triple,rep('1,2,3',3))
  h<-linf.udcst.graph(X[3:1,]);expect_equal(g$edges,h$edges)
  expect_equal(g$membership,h$membership[match(g$membership$sample.id,h$membership$sample.id),],ignore_attr=TRUE)
  s<-linf.udcst.graph(Matrix::Matrix(X,sparse=TRUE));expect_equal(g$edges,s$edges)
  shuffled<-X[,c(4,2,1,3)];realigned<-shuffled[,colnames(X)];expect_equal(linf.udcst.graph(realigned)$membership,g$membership)
})
