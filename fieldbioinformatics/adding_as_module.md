#Adding this repository as a submodule of another

Within the parent repo add the submodule:

```
git submodule add https://github.com/artic-network/<repo_name>.git
```

Commit the change and push:

```
git commit -m "adding submodule"
git push origin master
```

To update all submodules:

```
git submodule update --remote
```
