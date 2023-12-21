# How to migrate existing merge requests to auto-formatted code


1. Fetch the latest changes.

```shell
git fetch
```

2. Merge or rebase the changes before the auto-formatter commit. There might be conflicts at this stage, carefully review those changes and resolve the conflicts.

```shell
git merge f4b6e09884bb2ea2f9af1129bd5e1bf8f9ef2ffc
```

3. Install ruff.
```shell
pip install ruff==0.1.8
```

4. Auto-format and commit your code.
```shell
ruff format . && git commit -am "autoformat"
```

5. Merge the formatted code. This shouldn't have any conflicts.
```shell
git merge 5d5879a4b4f12fec1ebd0d133684190810b6a838 -X ours
```

6. Merge the develop branch. There might be conflicts at this stage, carefully review those changes and resolve the conflicts.
```shell
git merge develop
```

## Optional step:
1. Configure your git blame to ignore the formatting changes.
```shell
git config blame.ignoreRevsFile .git-blame-ignore-revs
```