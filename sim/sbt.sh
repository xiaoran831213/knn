## batched

# distribution change
c='N=512, P=10000, Q=10, R=2, efn=EST, eps=1.5, frq=.2, ubz=64, vcs=3';  d="sim/run/bd0"; rm $d/*.rds $d/{log,std,cms}/*
for i in {1..200}; do
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"ref\", lnk=SC, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"exp\", lnk=XP, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"bin\", lnk=BN, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"chi\", lnk=CH, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"poi\", lnk=PS, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
done | hpcwp - -d$d -q10 -m10 --wtm 4 --ln R,data --cp sim/sbt.R --ln sim,R --tag ${d##*/}

# high order response and variant interaction
c='N=512, P=10000, Q=10, R=2, efn=EST, eps=1.5, frq=.2, ubz=64, vcs=3';  d="sim/run/bn0"; rm $d/*.rds $d/{log,std,cms}/*
for i in {1..200}; do
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"ref\", lnk= SC, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"e^2\", lnk= O2, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"e^3\", lnk= O3, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"hy1\", lnk=HY1, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"rc1\", lnk=RC1, oks=~PL, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
done | hpcwp - -d$d -q10 -m10 --wtm 4 --ln R,data --cp sim/sbt.R --ln sim,R --tag ${d##*/}

# non-linear
c='N=512, P=10000, Q=10, R=2, efn=EST, eps=1.5, frq=.2, ubz=64, lnk=SC'; d="sim/run/bi0"; rm $d/*.rds $d/{log,std,cms}/*
for i in {1..200}; do
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"ref\", oks=~PL,    vcs=c(3, 0), ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"g:2\", oks=~PL+I2, vcs=c(2, 1), ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"g:3\", oks=~PL+I3, vcs=c(2, 1), ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"g*2\", oks=~PL+X2, vcs=c(2, 1), ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sbt.R\"); r=main(seed=$[1000+i], tag=\"g*3\", oks=~PL+X3, vcs=c(2, 1), ${c}); saveRDS(r, \"{n:05d}.rds\")'"
done | hpcwp - -d$d -q10 -m10 --wtm 4 --ln R,data --cp sim/sbt.R --ln sim,R --tag ${d##*/}
sim/run/bd0/sub.sh; sim/run/bi0/sub.sh; sim/run/bn0/sub.sh


# impact of sample size and batch size
c='N=1024, P=10000, R=1, eps=1.5, efn=EST, vcs=3, frq=.2, oks=~PL, lnk=O2, bpt=0'; d="sim/run/bz0"; rm $d/*.rds $d/{log,std,cms}/*
for i in {1..200}; do
    echo "time Rscript -e 'source(\"sb1.R\"); r=main(seed=$[1000+i], tag=\"n2k\", Q=2, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sb1.R\"); r=main(seed=$[1000+i], tag=\"n4k\", Q=4, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sb1.R\"); r=main(seed=$[1000+i], tag=\"n6k\", Q=6, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
    echo "time Rscript -e 'source(\"sb1.R\"); r=main(seed=$[1000+i], tag=\"n8k\", Q=8, ${c}); saveRDS(r, \"{n:05d}.rds\")'"
done | hpcwp - -d$d -q8 -m16 --wtm 4 --ln R,data --cp sim/sb1.R --ln sim,R --tag ${d##*/}
