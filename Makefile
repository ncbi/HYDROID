
clean-build:
	rm --force --recursive build/
	rm --force --recursive dist/
	rm --force --recursive *.egg-info

sdist:
	rm -rf dist
	rm -rf *.egg-info
	python setup.py sdist

pypi-push:
	twine upload dist/*

gh-push:
	git add .
	git commit -m "Update"
	git push

upd: gh-push sdist pypi-push