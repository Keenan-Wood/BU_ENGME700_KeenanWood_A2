- Navigate to the project directory

- With a fresh python environment active, run the command: pip install -e .

- To install dev dependencies, run the command: pip install -r requirements-dev.txt
  (this will install dependencies listed in requirements.txt as well)

- To test, run the command: pytest -v --cov=bisect_weighted --cov-report term-missing