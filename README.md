fastGLOBETROTTER is an updated version of the same GLOBETROTTER model, using the same input, but that is ~4-20 times faster than GLOBETROTTER without sacrificing accuracy. 

The bioRvix paper can be found [here](https://www.biorxiv.org/content/10.1101/2021.08.12.455263v1).

To download and compile:

```
git clone https://github.com/sahwa/fastGLOBETROTTER.git && cd fastGLOBETROTTER
Rscript -e 'install.packages(c("nnls", "data.table", "optparse"))'
R CMD SHLIB -o fastGLOBETROTTERCompanion.so fastGLOBETROTTERCompanion.c -lz -O3
```  
