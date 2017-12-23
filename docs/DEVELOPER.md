# HYDROID, Developer Docs

## What to do to ship a release:
- increment version file hydroid/VERSION
- tag release
- Modify anaconda recipe https://github.com/intbio/hydroid-conda
- Build conda package and deploy to anaconda cloud


### Deploy HYDROID to PyPi (no FreeSASA module)

```
make sdist
make pypi-push
```

.pipyrc should look like this
```
index-servers =
  pypi

[pypi]
repository=https://upload.pypi.org/legacy/
username=....
password=....
```


