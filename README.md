# haplotyping with SDhaP



I need to use [sdhap](https://sourceforge.net/projects/sdhap/) program but it is not easy to install. It has some dependencies. LAPACK (Linear Algebra PACKage) is written in Fortran. There are some LAPACK APIs for C.  SDhaP needs the c-LAPACK version of ATLAS (Automatically Tuned Linear Algebra Software). ATLAS also contains some others libraries. So installing ATLAS with the option (--with-netlib-lapack-tarfile) is enough. And it's not needed to install clapack or cblas or lapack seperately! 




# ATLAS

You first need to build [ATLAS](http://math-atlas.sourceforge.net/) (Automatically Tuned Linear Algebra Software) libraries and also need [lapack](http://www.netlib.org/lapack/#_software) using following instruction.


## Download

```
wget https://netix.dl.sourceforge.net/project/math-atlas/Stable/3.10.3/atlas3.10.3.tar.bz2

tar xjf atlas3.10.3.tar.bz2
mv ATLAS ATLAS_source


mkdir lapack_source; cd lapack_source
wget http://www.netlib.org/lapack/lapack-3.8.0.tar.gz 

cd ..
mkdir ATLAS_out; mkdir ATLAS_out/lib; mkdir ATLAS_out/include;
mkdir ATLAS_build
cd ATLAS_build
```


## Configuration


In the `ATLAS_source` folder, there is pdf file as the installation guid.  Step 3 The ATLAS configure  is very important. You may also use 
```
${PATH}='/mnt/scratch/majid001/sdhap_full_inst/'

${PATH}/ATLAS_source/configure --help
```

You need to check your architecture and  your cpu and type it as Mhz in option DPentiumCPS.
```
cat /proc/cpuinfo 
uname -m
```
for x86 use `-DPentiumCPS=2300` and for non x86 use  `-D c -DWALL`. Default option generates static library. You can add `--shared` to generate also shared library.


For generic libraries
```
${PATH}/ATLAS_source/configure --prefix=${PATH}/ATLAS_out --with-netlib-lapack-tarfile=${PATH}/lapack_source/lapack-3.8.0.tar.gz  -b 64 -V -1 -A x86x87 -t 2 -v 2
```


For your case
```
${PATH}/ATLAS_source/configure --prefix=${PATH}/ATLAS_out --with-netlib-lapack-tarfile=${PATH}/lapack_source/lapack-3.8.0.tar.gz  -b 64 -D c -DPentiumCPS=2300  -v 2
```


It takes several minutes.

## build 



```
make 
```
It may takes few days!


```
make check    
make ptcheck  
make time  
make install

```



# SDhaP



## build

The official line by author
```
gcc SDhaP_poly.c -o hap -L/usr/lib64/atlas -llapack -lcblas -I/usr/include/atlas
```

I'm using the followings for my mac and WUR fisher server. 
```
gcc SDhaP_poly.c -g -o out -L/ATLAS_build/lib -lcblas -llapack  -I/ATLAS_out/include -lm -latlas -lf77refblas -L/usr/local/Cellar/gcc/8.3.0_2/lib/gcc/8/ -lgfortran

gcc SDhaP_poly.c -g -o hap_polya -I/usr/include/x86_64-linux-gnu/ -llapack -lcblas -lm
```


## Run

```
k=3

./hap_poly fragment_file out.hap $k
```







## References

"SDhaP: haplotype assembly for diploids and polyploids via semi-definite programming", 
Shreepriya Das* and Haris Vikalo





