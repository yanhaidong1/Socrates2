#!/bin/bash
set -e

# ============================================================
# 一键环境构建脚本(通用版,不含任何镜像硬编码)
# 用法: bash Socrates2_one_step_envset_XXXX.sh
#
# 国内用户加速 Bioconductor 数据包下载:
#   此脚本不处理镜像(避免污染用户 condarc)。
#   如需加速,请按需运行独立工具:
#     bash patch_biocdata.sh [镜像URL, 默认清华]
# ============================================================

ENV_NAME="Socrates2"

mamba create -n $ENV_NAME --override-channels -c conda-forge -c bioconda \
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

# 自动定位 conda base(多级兜底: CONDA_EXE → command -v → 常见路径)
if [ -n "$CONDA_EXE" ] && [ -f "$(dirname "$(dirname "$CONDA_EXE")")/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="$(dirname "$(dirname "$CONDA_EXE")")"
elif command -v conda >/dev/null 2>&1 && [ -f "$(dirname "$(dirname "$(command -v conda)")")/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="$(dirname "$(dirname "$(command -v conda)")")"
elif [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="$HOME/miniconda3"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="$HOME/anaconda3"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
  CONDA_BASE="/opt/conda"
else
  echo "❌ 无法自动定位 conda, 请手动设置 CONDA_BASE"
  exit 1
fi
echo ">>> 使用 conda base: $CONDA_BASE"
source "$CONDA_BASE/etc/profile.d/conda.sh"
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


# --- 拆分第二批为小批,避免 libmamba 大求解卡死 ---
echo "第二批 2a/3: R 核心 + Seurat 类..."
mamba install -n $ENV_NAME --override-channels -c conda-forge -c bioconda \
  python=3.9 r-base=4.2.3 perl=5.32 libnsl=2.0.0 icu=73.2 r-stringi=1.8.3 bioconductor-biocversion=3.16.0 r-biocmanager=1.30.22 nano orthofinder=2.5.5 mcscanx \
  git make cmake gxx_linux-64 deeptools \
  r-mclust r-domc r-kknn r-scpred r-ggpubr \
  r-dbscan r-r.utils r-remotes \
  r-tidyverse r-rsample r-furrr r-argparse \
  r-seurat=5 r-seuratobject=5 \
  r-rcppml r-dosnow r-itertools r-glmnet r-viridis r-qlcmatrix \
  r-data.table r-vgam r-igraph r-dplyr \
  r-devtools r-harmony r-phytools r-speedglm \
  r-sm r-ggvenn r-vioplot \
  r-pheatmap r-ggalluvial -y || {
  echo "Additional packages 2a installation failed! "
  exit 1
}

echo "第二批 2b/3: Bioconductor 类..."
mamba install -n $ENV_NAME --override-channels -c conda-forge -c bioconda \
  python=3.9 r-base=4.2.3 perl=5.32 libnsl=2.0.0 icu=73.2 r-stringi=1.8.3 bioconductor-biocversion=3.16.0 r-biocmanager=1.30.22 bioconductor-biostrings bioconductor-complexheatmap \
  bioconductor-gviz bioconductor-genomicranges bioconductor-rtracklayer \
  bioconductor-summarizedexperiment bioconductor-preprocesscore \
  bioconductor-chromvar bioconductor-motifmatchr \
  bioconductor-jaspar2016 bioconductor-jaspar2018 bioconductor-jaspar2020 -y || {
  echo "Additional packages 2b installation failed! "
  exit 1
}

echo "第二批 2c/3: spatstat 类..."
mamba install -n $ENV_NAME --override-channels -c conda-forge -c bioconda \
  python=3.9 r-base=4.2.3 perl=5.32 libnsl=2.0.0 icu=73.2 r-stringi=1.8.3 bioconductor-biocversion=3.16.0 r-biocmanager=1.30.22 r-spatstat.core r-spatstat.geom r-spatstat.linnet r-spatstat.random r-spatstat.explore -y || {
  echo "Additional packages 2c installation failed! "
  exit 1
}

echo "第二批 2d: r-cmf(conda 版, 仅 pkgs/r 通道) + 显式 pin r-base=4.2.3 防 R 被顶升..."
mamba install -n $ENV_NAME -c conda-forge -c bioconda -c defaults r-cmf python=3.9 r-base=4.2.3 perl=5.32 libnsl=2.0.0 icu=73.2 r-stringi=1.8.3 bioconductor-biocversion=3.16.0 r-biocmanager=1.30.22 -y || {
  echo "r-cmf installation failed! "
  exit 1
}

echo "Reinstalling BSgenome from source to avoid UCSC tool compatibility issues..."

Rscript -e '
  options(repos = c(
    CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
    Bioc = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
  ))
  options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

  if ("BSgenome" %in% installed.packages()) {
    remove.packages("BSgenome")
    message("Removed conda-installed BSgenome package")
  }

    # 确保数据包存在(BSgenome 依赖 GenomeInfoDbData, 缺了 lazy loading 会失败)
  if (!requireNamespace("GenomeInfoDbData", quietly = TRUE)) {
    message("Installing GenomeInfoDbData (Bioc data package, required by BSgenome)...")
    BiocManager::install("GenomeInfoDbData", update = FALSE, ask = FALSE, force = TRUE, type = "source")
  } else {
    message("GenomeInfoDbData already installed: ", as.character(packageVersion("GenomeInfoDbData")))
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
  options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  if (!require("varistran", quietly = TRUE)) {
    # 优先 GitHub API, 失败则直接下载 tarball(绕过 api.github.com, 兼容国内服务器)
    ok <- tryCatch({
      remotes::install_github("MonashBioinformaticsPlatform/varistran", upgrade = "never")
      TRUE
    }, error = function(e) {
      message("GitHub API 不可达, 改用 tarball 直下: ", conditionMessage(e))
      FALSE
    })
    if (!ok) {
      remotes::install_url("https://github.com/MonashBioinformaticsPlatform/varistran/archive/refs/heads/master.tar.gz")
    }
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
  options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", quiet = TRUE)
  }

  message("Installing GENESPACE...")
  ok <- tryCatch({
    remotes::install_github(
      "jtlovell/GENESPACE",
      force = TRUE,
      quiet = TRUE,
      upgrade = "never",
      INSTALL_opts = c("--no-multiarch", "--with-keep.source")
    )
    TRUE
  }, error = function(e) {
    message("GitHub API 不可达, 改用 tarball 直下: ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    remotes::install_url("https://github.com/jtlovell/GENESPACE/archive/refs/heads/master.tar.gz", INSTALL_opts = c("--no-multiarch", "--with-keep.source"))
  }

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

echo "Installing qs2 from CRAN (conda r-qs2 requires R>=4.4, incompatible with R 4.2.3)..."
Rscript -e '
  options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
  if (!requireNamespace("qs2", quietly = TRUE)) {
    install.packages("qs2", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  }
  if (requireNamespace("qs2", quietly = TRUE)) {
    message("qs2 installed! Version: ", packageVersion("qs2"))
  } else {
    stop("qs2 installation failed!")
  }
' || {
  echo "qs2 installation failed!"
  exit 1
}

echo "编译 log_finite 垫片(兼容 glibc>=2.34,解决 macs2 __log_finite undefined)..."
CC="${CONDA_PREFIX}/bin/x86_64-conda-linux-gnu-cc"
if [ -x "$CC" ]; then
  cat > "${CONDA_PREFIX}/lib/log_finite_shim.c" <<'C_EOF'
/* Full __*_finite shim for glibc >= 2.34 (Ubuntu 22.04+, Rocky 9, ...)
 * Restores symbols removed from libm so old compiled extensions (e.g. MACS2)
 * keep loading. Wrappers call the plain functions (identical math results). */
#include <math.h>
#define SHIM_FN(ret, name, arglist, call) ret name arglist { return call; }
SHIM_FN(double, __acos_finite, (double x), acos(x))
SHIM_FN(double, __acosh_finite, (double x), acosh(x))
SHIM_FN(double, __asin_finite, (double x), asin(x))
SHIM_FN(double, __asinh_finite, (double x), asinh(x))
SHIM_FN(double, __atan2_finite, (double y, double x), atan2(y, x))
SHIM_FN(double, __atan_finite, (double x), atan(x))
SHIM_FN(double, __atanh_finite, (double x), atanh(x))
SHIM_FN(double, __cos_finite, (double x), cos(x))
SHIM_FN(double, __cosh_finite, (double x), cosh(x))
SHIM_FN(double, __exp2_finite, (double x), exp2(x))
SHIM_FN(double, __exp_finite, (double x), exp(x))
SHIM_FN(double, __expm1_finite, (double x), expm1(x))
SHIM_FN(double, __fmod_finite, (double x, double y), fmod(x, y))
SHIM_FN(double, __hypot_finite, (double x, double y), hypot(x, y))
SHIM_FN(double, __log10_finite, (double x), log10(x))
SHIM_FN(double, __log1p_finite, (double x), log1p(x))
SHIM_FN(double, __log2_finite, (double x), log2(x))
SHIM_FN(double, __log_finite, (double x), log(x))
SHIM_FN(double, __pow_finite, (double x, double y), pow(x, y))
SHIM_FN(double, __remainder_finite, (double x, double y), remainder(x, y))
SHIM_FN(double, __sin_finite, (double x), sin(x))
SHIM_FN(double, __sinh_finite, (double x), sinh(x))
SHIM_FN(double, __sqrt_finite, (double x), sqrt(x))
SHIM_FN(double, __tan_finite, (double x), tan(x))
SHIM_FN(double, __tanh_finite, (double x), tanh(x))
SHIM_FN(float, __acosf_finite, (float x), acosf(x))
SHIM_FN(float, __asinf_finite, (float x), asinf(x))
SHIM_FN(float, __atan2f_finite, (float y, float x), atan2f(y, x))
SHIM_FN(float, __atanf_finite, (float x), atanf(x))
SHIM_FN(float, __cosf_finite, (float x), cosf(x))
SHIM_FN(float, __exp2f_finite, (float x), exp2f(x))
SHIM_FN(float, __expf_finite, (float x), expf(x))
SHIM_FN(float, __fmodf_finite, (float x, float y), fmodf(x, y))
SHIM_FN(float, __hypotf_finite, (float x, float y), hypotf(x, y))
SHIM_FN(float, __log10f_finite, (float x), log10f(x))
SHIM_FN(float, __log2f_finite, (float x), log2f(x))
SHIM_FN(float, __logf_finite, (float x), logf(x))
SHIM_FN(float, __powf_finite, (float x, float y), powf(x, y))
SHIM_FN(float, __sinf_finite, (float x), sinf(x))
SHIM_FN(float, __sqrtf_finite, (float x), sqrtf(x))
SHIM_FN(float, __tanf_finite, (float x), tanf(x))
C_EOF
  "$CC" -shared -fPIC -O2 -o "${CONDA_PREFIX}/lib/liblog_finite_shim.so" "${CONDA_PREFIX}/lib/log_finite_shim.c" -lm
  # conda 激活钩子: 激活此环境时自动设置 LD_PRELOAD
  mkdir -p "${CONDA_PREFIX}/etc/conda/activate.d" "${CONDA_PREFIX}/etc/conda/deactivate.d"
  cat > "${CONDA_PREFIX}/etc/conda/activate.d/zz_log_finite.sh" <<EOF
export LD_PRELOAD="${CONDA_PREFIX}/lib/liblog_finite_shim.so:\${LD_PRELOAD}"
EOF
  cat > "${CONDA_PREFIX}/etc/conda/deactivate.d/zz_log_finite.sh" <<EOF
unset LD_PRELOAD
EOF
  echo "✓ log_finite 垫片已编译,LD_PRELOAD 激活钩子已设置"
else
  echo "⚠ 未找到 conda gcc,跳过垫片编译(如遇 macs2 __log_finite 问题需手动处理)"
fi

echo "Environment ${ENV_NAME} is ready!"
