#!/bin/bash
set -e

ENV_NAME="py39_r42_CR2"

mamba create -n $ENV_NAME -c conda-forge -c bioconda -c defaults \
  python=3.9 r-base=4.2.3 perl=5.32 libnsl=2.0.0 icu=73.2 \
  libgcc-ng libstdcxx-ng libgfortran-ng openblas \
  r-stringi=1.8.3 bioconductor-biocversion=3.16.0 r-biocmanager=1.30.22 bioconductor-topgo \
  bioconductor-genomicfeatures r-mass r-matrix=1.6_5 r-ggplot2 r-seriation bioconductor-edgeR \
  biopython macs2 numpy pandas pyyaml pysam samtools threadpoolctl scipy idr \
  scikit-learn imbalanced-learn picard bedtools natsort \
  ucsc-bedgraphtobigwig ucsc-fatotwobit \
  perl-sort-naturally meme -y || {
  echo "Conda environment setup failed! "
  exit 1
}

source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME || {
  echo "Failed to activate environment!"
  exit 1
}

export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$CONDA_PREFIX/lib64"
export PKG_CONFIG_PATH="$CONDA_PREFIX/lib/pkgconfig:$CONDA_PREFIX/share/pkgconfig"

LIB_PATH="${CONDA_PREFIX}/lib"
if [[ -f "${LIB_PATH}/libnsl.so.3.0.0" && ! -e "${LIB_PATH}/libnsl.so.1" ]]; then
    echo "Creating libnsl.so.1 symlink..."
    ln -s "${LIB_PATH}/libnsl.so.3.0.0" "${LIB_PATH}/libnsl.so.1"
else
    echo "libnsl.so.1 already exists or base file missing. Skipping link creation."
fi


mamba install -n $ENV_NAME -c conda-forge -c bioconda -c defaults \
  r-base=4.2.3 icu=73.2 nano orthofinder=2.5.5 mcscanx \
  git make gxx_linux-64 deeptools \
  r-mclust r-domc r-kknn r-scpred r-ggpubr \
  r-dbscan r-r.utils bioconductor-biostrings \
  r-tidyverse bioconductor-complexheatmap r-rsample r-furrr r-argparse \
  r-spatstat.core r-spatstat.geom r-spatstat.linnet r-spatstat.random r-spatstat.explore \
  r-seurat=5 r-seuratobject=5 r-remotes \
  r-rcppml r-dosnow r-itertools r-glmnet r-viridis r-qlcmatrix r-cmf \
  bioconductor-gviz bioconductor-genomicranges bioconductor-rtracklayer \
  bioconductor-summarizedexperiment r-data.table r-vgam r-igraph r-dplyr \
  r-devtools r-harmony r-phytools r-speedglm \
  bioconductor-preprocesscore r-sm r-ggvenn r-vioplot \
  bioconductor-chromvar bioconductor-motifmatchr \
  bioconductor-jaspar2016 bioconductor-jaspar2018 bioconductor-jaspar2020 \
  r-pheatmap r-ggalluvial -y || {
  echo "Additional packages installation failed! "
  exit 1
}

echo "Reinstalling BSgenome from source to avoid UCSC tool compatibility issues..."

Rscript -e '
  options(repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
  ))

  if ("BSgenome" %in% installed.packages()) {
    remove.packages("BSgenome")
    message("Removed conda-installed BSgenome package")
  }

  BiocManager::install("BSgenome", update = FALSE, ask = FALSE, force = TRUE, type = "source")
  
  if (requireNamespace("BSgenome", quietly = TRUE)) {
    message("✓ BSgenome reinstalled successfully! Version: ", packageVersion("BSgenome"))
  } else {
    stop("✗ BSgenome reinstallation failed!")
  }
' || {
  echo "Failed to reinstall BSgenome!"
  exit 1
}

Rscript -e '
  options(repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
  ))

  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  if (!require("varistran", quietly = TRUE)) {
    remotes::install_github("MonashBioinformaticsPlatform/varistran")
  }
' || {
  echo "varistran installation failed!"
  exit 1
}


echo "Installing GENESPACE..."
Rscript -e '
  options(repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
  ))

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", quiet = TRUE)
  }

  message("Installing GENESPACE...")
  remotes::install_github(
    "jtlovell/GENESPACE", 
    force = TRUE, 
    quiet = TRUE,
    INSTALL_opts = c("--no-multiarch", "--with-keep.source")
  )

  if (requireNamespace("GENESPACE", quietly = TRUE)) {
    message("✓ GENESPACE installed successfully! Version: ", packageVersion("GENESPACE"))
    
    suppressWarnings(library(GENESPACE))
    
    message("Testing basic functionality...")
    essential_funcs <- c("init_genespace", "run_orthofinder", "plot_riparian")
    
    test_results <- sapply(essential_funcs, function(func) {
      tryCatch({
        exists(func, where = as.environment("package:GENESPACE"))
      }, error = function(e) {
        message("Error checking function ", func, ": ", e$message)
        FALSE
      })
    })
    
    if (all(test_results)) {
      message("✓ All core functions available")
    } else {
      warning("✗ Missing functions: ", paste(names(test_results)[!test_results], collapse = ", "))
    }
    
  } else {
    stop("✗ GENESPACE installation failed!")
  }

  message("GENESPACE installation completed!")
' || {
  echo "GENESPACE installation failed!"
  exit 1
}

echo "Environment ${ENV_NAME} is ready!"
