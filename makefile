test:
	@echo Running Tests
	@python -m unittest discover

docs:
	doxygen doxygen_config

clean:
	rm -r docs/