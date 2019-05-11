# haplotyping using illumina data


I need to used SDhaP program but it is not easy to use.


For installation ATLAS, is enough with the following configuration


```

tar xjf atlas3.10.3.tar.bz2
mv ATLAS ATLAS_source
mkdir ATLAS_out; mkdir ATLAS_out/lib; mkdir ATLAS_out/include;
mkdir ATLAS_build
cd ATLAS_build

YOURPATH= .
${YOURPATH}/ATLAS_source/configure   --prefix=/mnt/scratch/majid001/c/ATLAS_out --with-netlib-lapack-tarfile=/mnt/scratch/majid001/c/lapack-3.8.0.tar.gz   -b 64 -D c -DPentiumCPS=2300

#--cripple-atlas-performance 
#-D c -DWALL 

make 
make check    
make ptcheck  
make time  
make install

```

Some notes:

LAPACK (Linear Algebra PACKage) is written in Fortran. There are some LAPACK APIs for C.  SDhaP needs the c-LAPACK version of ATLAS (Automatically Tuned Linear Algebra Software). ATLAS also contains some others libraries.  

So it seems that installing ATLAS with the option (--with-netlib-lapack-tarfile) is enough. And no need to install clapack or cblas or lapack seperately! 



I'm also trying to make the code better!


## References
```

"SDhaP: haplotype assembly for diploids and polyploids via semi-definite programming
Shreepriya Das* and Haris Vikalo"
```

and my recent paper: 
```
"NGS based haplotype assembly using matrix completion", PLOS ONE, 2018
```




