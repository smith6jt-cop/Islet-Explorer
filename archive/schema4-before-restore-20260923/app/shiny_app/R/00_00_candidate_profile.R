Sys.setenv(ISLET_DATA_ROOT = "/data/panc_CODEX/Islet-Explorer/data/development/all-donor-20260918",
LOCAL_IMAGE_ROOT = "/data/panc_CODEX/Islet-Explorer/data/development/all-donor-20260918/images",
ISLET_RUNTIME_MODE = "development",
ISLET_DATA_SCHEMA = "3",
RETICULATE_PYTHON = "/home/smith6jt/miniconda3/envs/scvi-env/bin/python",
ISLET_ENABLE_AI = "FALSE",
SHOW_STATISTICS_TAB = "0",
SHOW_INFLAMMATION_TAB = "0",
DEBUG_CREDENTIALS = "0",
VIEWER_DEBUG = "0",
OPENBLAS_NUM_THREADS = "1",
OMP_NUM_THREADS = "1")

# Route direct workstation app links to nginx for byte-range TIFF serving.
options(islet.network.frontend_port = "8080", islet.network.backend_port = "3838")
