ENV_NAME=panc-codex-shiny

.PHONY: env create update run

env:
	conda env create -f environment.yml || conda env update -f environment.yml --prune

create:
	conda env create -f environment.yml

update:
	conda env update -f environment.yml --prune

run:
	conda run -n $(ENV_NAME) ./run_shiny.sh
