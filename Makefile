ENV_NAME=panc-codex-streamlit

.PHONY: env create update run

env:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

create:
	conda env create -f environment.yml

update:
	conda env update -f environment.yml --prune

run:
	conda run -n $(ENV_NAME) streamlit run app_v2.py --server.address 0.0.0.0 --server.port $${PORT:-8501}

