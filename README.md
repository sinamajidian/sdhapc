# haplotyping using illumina data


I need to used SDhaP program but it is not easy to use.


For installation  ATLAS is enough and following configuration
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



I'm also trying to make the code better!


## References
```

"SDhaP: haplotype assembly for diploids and polyploids via semi-definite programming
Shreepriya Das* and Haris Vikalo"

and 

"NGS based haplotype assembly using matrix completion", PLOS ONE, 2018
```




