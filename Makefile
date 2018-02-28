flake:
	@if command -v flake8 > /dev/null; then \
		echo "Running flake8"; \
		flake8 flake8 --ignore N802,N806 `find . -name \*.py | grep -v setup.py | grep -v /doc/`; \
	else \
		echo "flake8 not found, please install it!"; \
		exit 1; \
	fi;
	@echo "flake8 passed"

test:
	py.test --pyargs celloutline --cov-report term-missing --cov=celloutline

## DOCKER TASKS

# Import config.
# Can change default config with `make cnf="config_special.env" build`
cnf ?= config.env
include $(cnf)
export $(shell sed 's/=.*//' $(cnf))

# Build the container
build: ## Build the container
	docker build -t $(APP_NAME) .

build_dl: ## Build the dl container, hosted remotely
	docker build -t $(DL_APP_NAME) $(DL_REMOTE_REPO)

build-nc: ## Build the container without caching
	docker build --no-cache -t $(APP_NAME) .

run: ## Run container 
	# Default to port configured in `config.env`
	# Override with `make run PORT=????` 
	docker run -i -t --rm --env-file=./config.env -p=$(PORT):8888 -v $(LOCAL_MOUNT):/home/$(USER)/nfs --ipc=host $(APP_NAME)

run_dl: ## Run container 
	docker run -i -t --rm --env-file=./config.env -p=$(PORT):8888 -v $(LOCAL_MOUNT):/home/$(DL_USER)/nfs --ipc=host $(DL_APP_NAME)

up: build run ## Run container on port configured in `config.env` (Alias to run)

stop: ## Stop and remove a running container
	docker stop $(APP_NAME); docker rm $(APP_NAME)

