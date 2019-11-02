# ANOVA VIA LM OBJECT



sesamesim$site <- as.factor(sesamesim$site)
anov <- lm(sesamesim$postnumb~sesamesim$site-1)
set.seed(100)
z<-bain(anov, "site1=site2=site3=site4=site5;
      site2>site5>site1>site3=site4;
     site1=site2>site3=site4>site5;
     site1<site2>site3<site4>site5;
     site1=site5>site3=site4<site2;
     site2>site3>site4;
     (site1,site2,site5)>(site3,site4);
     site2>(site1,site3,site4,site5)")

z_gor<-gorica(anov, "site1=site2=site3=site4=site5;
      site2>site5>site1>site3=site4;
     site1=site2>site3=site4>site5;
     site1<site2>site3<site4>site5;
     site1=site5>site3=site4<site2;
     site2>site3>site4;
     (site1,site2,site5)>(site3,site4);
     site2>(site1,site3,site4,site5)")

test_that("gorica and bain give similar results for anova", {
   expect_equivalent(z_gor$fit$gorica_weights, z$fit$PMPb, tolerance = .2)
})
