# scitemplate

If you don't already have twine installed run the following:
`python -m pip install --user --upgrade twine`

## The following will create the package
```
python setup.py sdist bdist_wheel
twine check dist/PROJECT_NAME.tar.gz
```

## Install to python environment localling

`pip install PATH_TO_PROJECT/dist/PROJECT_NAME.tar.gz`
You should run this before uploading it and check all works as expected.

## The following will push the package to pip 
**Note you need to set up a pip account first**

```
twine upload dist/*
```

## Have a look at your projects page on pip

`https://pypi.org/project/PROJECT_NAME/`


## Installing semantic stuff

- install using package.json (npm install)
- update names in .releaserc
- Add personal access token to repo --> go to github settings --> developer (at the bottom) (https://github.com/settings/tokens/new)
- add a new personal access token (no exiry) click repo and click workflow & write packages. (then egnerate token)
- Go back to your repo on github and select settings, go to secrets, actions, then add new repository secret. Add in the secret with the name: PERSONAL_ACCESS_TOKEN
- create new branch gh-pages remove all files and commit 
- know how to do this: https://www.conventionalcommits.org/en/v1.0.0/#summary
-  Change into the docs/source directory and run : ln -s ../../../CHANGELOG.md index.md 